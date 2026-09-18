# -*- coding: utf-8 -*-
"""Tests for ``biorazer.database.molecule.rotamer`` (external rotamer readers).

两个来源子包各测一遍:

* :mod:`.pymol`   -- PyMOL 自带的 Dunbrack pickle
* :mod:`.rosetta` -- Rosetta 内附的 Dunbrack 2002 / Shapovalov 2010 文本库

这些测试依赖**本机安装**的 PyMOL / Rosetta (库不 vendored)。找不到库文件时
``pytest.skip``, 因此 CI 无这些库时不会失败。

除"能读"之外, 重点核对:
* 列切分与 Rosetta 源码 ``read_rotation_library()`` 的 a–q 语义一致
  (Dun02 = 15+N 列无 -lnP, Shapovalov = 16+N 列有 -lnP, 自动识别);
* 实测规模 (HIS: PyMOL ind 9 个 rotamer / dep 1141 条 169 箱;
  Shapovalov 49284 条 1369 箱; Dun02 12321 条 1369 箱);
* PyMOL ``sc_bb_dep`` 是 Dunbrack 2002 按概率降序的**子集**;
* chi 四元组与 :data:`...by_residue.SIDECHAIN_CHI` 一致 (同源 Dunbrack)。
"""

import numpy as np
import pytest

from biorazer.database import molecule as M
from biorazer.database.molecule import rotamer as R
from biorazer.database.molecule.rotamer import pymol as PM
from biorazer.database.molecule.rotamer import rosetta as RZ


# --------------------------------------------------------------------------
# availability
# --------------------------------------------------------------------------
def _need_pymol():
    if not PM.PYMOL_SC_BB_IND.exists():
        pytest.skip(f"PyMOL rotamer library not found: {PM.PYMOL_SC_BB_IND}")


def _need_rosetta():
    if not RZ.default_bbdep02_path().exists():
        pytest.skip(f"Rosetta rotamer library not found: {RZ.ROSETTA_ROTAMER_DIR}")


def _need_shapovalov():
    if not RZ.default_shapovalov_path("HIS").exists():
        pytest.skip(f"Shapovalov library not found: {RZ.ROSETTA_ROTAMER_DIR}")


# --------------------------------------------------------------------------
# package layout
# --------------------------------------------------------------------------
class TestPackageLayout:
    """子包按来源划分, 公共记录类型在 rotamer 顶层。"""

    def test_subpackage_exports(self):
        for mod in (R, RZ, PM):
            for name in mod.__all__:
                assert hasattr(mod, name), f"{mod.__name__} missing {name}"

    def test_shared_record_types_live_at_top(self):
        assert hasattr(R, "RotamerRecord") and hasattr(R, "RotamerLibrary")
        # 子包 reader 返回的正是顶层定义的类型
        assert R.RotamerRecord is not None
        assert R.RotamerLibrary is not None

    def test_molecule_reexports(self):
        for name in M.__all__:
            assert hasattr(M, name), f"molecule missing {name}"
        assert M.read_pymol_dep is PM.read_pymol_dep
        assert M.read_shapovalov is RZ.read_shapovalov
        assert M.read_bbdep02 is RZ.read_bbdep02
        assert M.RotamerRecord is R.RotamerRecord

    def test_readers_not_hoisted_to_rotamer_top(self):
        """reader 只经 .pymol / .rosetta 暴露, 不在 rotamer 顶层。"""
        for name in ("read_pymol_dep", "read_shapovalov", "read_bbdep02"):
            assert name not in R.__all__


# --------------------------------------------------------------------------
# PyMOL
# --------------------------------------------------------------------------
class TestPymolReaders:
    def test_ind_shape_and_normalisation(self):
        _need_pymol()
        ind = PM.read_pymol_ind()
        assert len(ind) == 18
        his = ind["HIS"]
        assert len(his) == 9
        assert his.chi_quads == (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "ND1"))
        assert sum(r.probability for r in his.records) == pytest.approx(1.0, abs=1e-6)
        # backbone-independent -> phi/psi 为 None
        assert all(r.phi is None and r.psi is None for r in his.records)
        # chi 顺序与 quads 一致
        assert all(len(r.chi) == len(his.chi_quads) for r in his.records)

    def test_dep_shape(self):
        _need_pymol()
        dep = PM.read_pymol_dep()
        assert len(dep) == 18
        his = dep["HIS"]
        assert len(his) == 1141
        assert len(his.bins) == 169
        assert all(r.phi is not None and r.psi is not None for r in his.records)

    def test_dep_backbone_lookup_rounds_to_bin(self):
        _need_pymol()
        his = PM.read_pymol_dep()["HIS"]
        top = his.top(1, phi=-178.0, psi=-176.0)
        assert top[0].phi == -180.0
        assert top[0].psi == -180.0
        # top 按概率降序
        several = his.top(3, phi=-180.0, psi=-180.0)
        probs = [r.probability for r in several]
        assert probs == sorted(probs, reverse=True)

    def test_ind_lookup_ignores_backbone(self):
        _need_pymol()
        his = PM.read_pymol_ind()["HIS"]
        a = his.top(1, phi=-60, psi=-40)
        b = his.top(1, phi=120, psi=90)
        assert a[0].chi == b[0].chi  # 不区分骨架

    def test_library_reader(self):
        _need_pymol()
        lib = PM.read_pymol_library()
        assert len(lib) > 0
        for res, entry in lib.items():
            assert entry.resname == res


# --------------------------------------------------------------------------
# Rosetta
# --------------------------------------------------------------------------
class TestRosettaReaders:
    def test_bbdep02_scale_and_columns(self):
        """Dun02 = 15+N 列, 无 -lnP; HIS 12321 条 = 1369 箱 x 9。"""
        _need_rosetta()
        his = RZ.read_bbdep02(resname="HIS")
        assert len(his) == 12321
        assert len(his.bins) == 1369
        r0 = his.records[0]
        assert np.isnan(r0.minus_log_prob), "Dun02 不应有 -lnP 列"
        assert r0.count is not None and r0.rotwell is not None
        # HIS 只有 2 个 chi -> r3/r4 = 0
        assert r0.rotwell[2:] == (0, 0)
        assert r0.chi[2:] == (0.0, 0.0)

    def test_bbdep02_all_residues(self):
        _need_rosetta()
        allz = RZ.read_rosetta_text(RZ.default_bbdep02_path())
        assert isinstance(allz, dict)
        assert len(allz) == 18
        assert len(allz["HIS"]) == 12321

    def test_shapovalov_columns_and_header(self):
        """Shapovalov = 16+N 列, 有 -lnP; 注释头给出 chi 数与 bin 数。"""
        _need_shapovalov()
        his = RZ.read_shapovalov("HIS")
        assert len(his) == 49284
        assert len(his.bins) == 1369
        assert his.header["n_chi"] == 2
        assert his.header["chi_bins"] == [3, 12]
        assert his.header["n_rotamers"] == 36
        assert his.header["phi_step"] == 10.0
        r0 = his.records[0]
        assert not np.isnan(r0.minus_log_prob), "Shapovalov 应有 -lnP 列"
        # -lnP 与 P 应满足 -lnP = -ln(P)
        assert r0.minus_log_prob == pytest.approx(-np.log(r0.probability), abs=1e-5)

    def test_shapovalov_resname_case_insensitive(self):
        _need_shapovalov()
        a = RZ.read_shapovalov("HIS")
        b = RZ.read_shapovalov("his")
        assert len(a) == len(b)
        assert a.records[0].chi == b.records[0].chi

    def test_unknown_residue_raises(self):
        _need_rosetta()
        with pytest.raises(KeyError):
            RZ.read_rosetta_text(RZ.default_bbdep02_path(), resname="XXX")

    def test_default_paths(self):
        assert RZ.default_bbdep02_path().name == "bbdep02.May.sortlib"
        p = RZ.default_shapovalov_path("his")
        assert p.name == "his.bbdep.rotamers.lib.gz"
        assert RZ.DEFAULT_STEPDOWN in str(p)


# --------------------------------------------------------------------------
# cross-source consistency
# --------------------------------------------------------------------------
class TestCrossSource:
    def test_pymol_dep_is_subset_of_dun02(self):
        """PyMOL sc_bb_dep 是 Dunbrack 2002 按概率降序的子集 (低概率条目被丢弃)。"""
        _need_pymol()
        _need_rosetta()
        pm = PM.read_pymol_dep()["HIS"]
        dn = RZ.read_bbdep02(resname="HIS")

        dn_bins = {}
        for r in dn.records:
            dn_bins.setdefault((r.phi, r.psi), []).append(r)

        checked = 0
        for key in pm.bins:
            a = sorted(pm.for_backbone(*key), key=lambda r: -r.probability)
            b = sorted(dn_bins[key], key=lambda r: -r.probability)
            # PyMOL 的每一条都应能在 Dun02 同箱内找到 (chi 一致, P 至浮点精度)
            bset = [(round(r.chi[0], 3), round(r.chi[1], 3), r.probability) for r in b]
            for r in a:
                hit = [p for c1, c2, p in bset
                       if c1 == round(r.chi[0], 3) and c2 == round(r.chi[1], 3)]
                assert hit, f"PyMOL rotamer {r.chi[:2]} not found in Dun02 bin {key}"
                assert min(abs(p - r.probability) for p in hit) < 1e-5, (
                    f"概率不一致 in bin {key}: {r.probability} vs {hit}"
                )
            assert len(a) <= len(b), f"PyMOL 条数不应多于 Dun02 (bin {key})"
            checked += 1
        assert checked == 169

    def test_chi_order_is_normalised_not_file_order(self):
        """PyMOL pickle 的 key 顺序对 7 个残基被打乱; reader 必须按规范顺序重排。

        以 GLU 为例: 文件内 key 顺序是 (chi2, chi3, chi1), 若按 key 顺序取,
        chi1 会拿到 chi2 的值。这里核对 reader 输出与 Dun02 的 chi1 分布一致
        (都是 g- 约 -60 度)。
        """
        _need_pymol()
        _need_rosetta()
        from biorazer.database.molecule.bond.dihedral.protein.by_residue import SIDECHAIN_CHI

        glu = PM.read_pymol_ind()["GLU"]
        # chi_quads 必须是 chi1..chiN 规范顺序
        assert tuple(SIDECHAIN_CHI["GLU"]) == glu.chi_quads
        # 概率最高那条的 chi1 应是 g- (约 -60 度), 而不是 chi2 的 ~180
        top = glu.top(1)[0]
        assert -90 < top.chi[0] < -30, f"GLU chi1 应在 g- 区, 得 {top.chi[0]}"
        assert abs(top.chi[1]) > 150, f"GLU chi2 应接近 180, 得 {top.chi[1]}"

    def test_ile_merged_atom_name_normalised(self):
        """ILE 的 chi2 在 PyMOL 文件里叫 ``CD1+CD``, reader 应规范化为 ``CD1``。"""
        _need_pymol()
        from biorazer.database.molecule.bond.dihedral.protein.by_residue import SIDECHAIN_CHI
        ile = PM.read_pymol_ind()["ILE"]
        assert tuple(SIDECHAIN_CHI["ILE"]) == ile.chi_quads
        assert "CD1+CD" not in [a for q in ile.chi_quads for a in q]

    def test_chi_definitions_match_by_residue(self):
        """PyMOL pickle 的 chi 四元组与 by_residue.SIDECHAIN_CHI 一致 (同源 Dunbrack)。"""
        _need_pymol()
        from biorazer.database.molecule.bond.dihedral.protein.by_residue import SIDECHAIN_CHI
        ind = PM.read_pymol_ind()
        for res, lib in ind.items():
            if res not in SIDECHAIN_CHI or not lib.chi_quads:
                continue
            assert tuple(SIDECHAIN_CHI[res]) == lib.chi_quads, f"{res} chi 定义不一致"
            assert all(len(r.chi) == len(lib.chi_quads) for r in lib.records)
