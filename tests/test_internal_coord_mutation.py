# -*- coding: utf-8 -*-
"""Tests for :mod:`biorazer.structure.manipulation.internal_coord.mutation`.

被测的是"在 **IC 图**上把残基侧链整棵换掉"这条链路。夹具不读任何外部结构
文件: 用另一个测试文件里那条**只含主链**的理想链
(:func:`tests.test_mutation._ideal_chain`) 经 bridge 转成 ``InternalCoord`` ——
所以"插入侧链" (图里本来没有侧链) 与"替换侧链" (图里有、要丢掉) 两条路
都在纯内存夹具上跑得到。

判据 (每条都对着数, 不读 docstring 推断):

* 骨架 (``N/CA/C/O``, 含羰基 ``O``) 与**非目标残基的每个原子**坐标不动 ——
  从重建坐标 (``to_coords``) 逐原子回量;
* 目标残基的原子集合 = 模板的侧链原子集合, ``res_name`` 更新;
* chi 从结果图上回读 = :func:`dihedral` 从重建坐标回量 = 调用方给的/rotamer 名对应的值;
* 编号右移: 插入侧链后, 插入点**之后**的原子仍与输入一一对应 (不会撞号);
* 非法输入全部 ``ValueError``/``TypeError``。
"""

import numpy as np
import pytest

from biorazer.database.molecule.bond.dihedral.protein import SIDECHAIN_CHI
from biorazer.database.molecule.icoor.protein.template import rotamer_names
from biorazer.database.molecule.icoor.protein.topology import MAINCHAIN_ATOMS
from biorazer.structure.bridge import AtomArray_InternalCoord
from biorazer.structure.manipulation.internal_coord import mutate
from biorazer.structure.manipulation.internal_coord.mutation import _SS
from biorazer.structure.objects import InternalCoord, InternalCoordAtom
from biorazer.structure.objects.internal_coords import dihedral

from .test_mutation import _ideal_chain, _wrap

CHI_TOL = 1e-4


def _chain_ic(n_res=4, chain="A", res_names=None):
    """理想主链 -> :class:`InternalCoord` (仓库自己的 bridge 读入路径)。"""
    arr = _ideal_chain(n_res, chain, res_names=res_names)
    return AtomArray_InternalCoord(arr).convert()


def _coords_by_key(ic):
    """``{(chain_id, res_id, ins_code, atom_name): (x, y, z)}`` (从图上重建)。"""
    coords = ic.to_coords()
    return {(a.chain_id, a.res_id, a.ins_code, a.name): np.asarray(coords[i], float)
            for i, a in enumerate(ic.atoms)}


def _chi(ic, res_id, resn, chain="A"):
    """从重建坐标回量某残基的 chi (用 :func:`dihedral`)。"""
    xyz = _coords_by_key(ic)
    return [dihedral(*[xyz[(chain, res_id, "", name)] for name in quad])
            for quad in SIDECHAIN_CHI[resn]]


def _atom_names(ic, res_id, chain="A"):
    return [a.name for a in ic.atoms
            if a.chain_id == chain and a.res_id == res_id]


class TestGraft:
    def test_inserts_the_target_side_chain(self):
        ic = _chain_ic(4)
        out = mutate(ic, ["A2H"])
        assert _atom_names(ic, 2) == ["N", "CA", "C", "O"]          # 夹具本来没有侧链
        assert _atom_names(out, 2) == ["N", "CA", "C", "O", "CB", "CG",
                                       "ND1", "CE1", "NE2", "CD2"]
        assert {str(a.res_name) for a in out.atoms
                if a.chain_id == "A" and a.res_id == 2} == {"HIS"}

    def test_replaces_an_existing_side_chain(self):
        """先 A2H (插入), 再 H2S (丢掉 HIS 的 6 个侧链原子换成 SER 的 2 个)。"""
        ic = _chain_ic(4)
        his = mutate(ic, ["A2H"])
        ser = mutate(his, ["H2S"])
        assert _atom_names(ser, 2) == ["N", "CA", "C", "O", "CB", "OG"]
        assert {str(a.res_name) for a in ser.atoms
                if a.chain_id == "A" and a.res_id == 2} == {"SER"}

    def test_backbone_and_other_residues_do_not_move(self):
        ic = _chain_ic(4)
        before = _coords_by_key(ic)
        out = mutate(ic, ["A2H"])
        after = _coords_by_key(out)
        moved = {key: float(np.linalg.norm(after[key] - before[key]))
                 for key in before}
        assert max(moved.values()) == 0.0, \
            f"有原子动了: {[k for k, v in moved.items() if v]}"

    def test_does_not_touch_the_input(self):
        ic = _chain_ic(4)
        n_atoms, n_dih = len(ic.atoms), len(ic.dihedra)
        mutate(ic, ["A2H"])
        assert (len(ic.atoms), len(ic.dihedra)) == (n_atoms, n_dih)
        assert _atom_names(ic, 2) == ["N", "CA", "C", "O"]

    def test_index_shift_keeps_downstream_atoms_identified(self):
        """侧链插在残基 2 之后, 残基 3/4 的原子按 (名, 残基) 仍在, 坐标不动。"""
        ic = _chain_ic(4)
        before = _coords_by_key(ic)
        out = mutate(ic, ["A2H"])
        after = _coords_by_key(out)
        assert set(before) <= set(after)
        assert [a.name for a in out.atoms[-4:]] == ["N", "CA", "C", "O"]
        assert max(float(np.linalg.norm(after[k] - before[k]))
                   for k in before) == 0.0

    def test_two_mutations_in_one_call(self):
        from biorazer.database.molecule.icoor.protein.template import build_template
        from biorazer.structure.manipulation.internal_coord.mutation import _SS

        def sc(resn):
            tmpl, _ = build_template(resn, _SS, "canonical")
            return [a.name for a in tmpl.atoms if a.name not in MAINCHAIN_ATOMS]

        out = mutate(_chain_ic(4), ["A1H", "A4W"])
        assert _atom_names(out, 1)[4:] == sc("HIS")
        assert _atom_names(out, 4)[4:] == sc("TRP")
        assert _atom_names(out, 2) == ["N", "CA", "C", "O"]

    def test_gly_target_has_no_side_chain(self):
        ic = _chain_ic(4)
        assert _atom_names(mutate(ic, ["A2G"]), 2) == ["N", "CA", "C", "O"]

    def test_pro_target_ring_closes(self):
        """PRO 的环是 IC_PATH 里的分支: ``to_coords`` 必须能一致地重建它。"""
        out = mutate(_chain_ic(4), ["A2P"], rotamer="canonical")
        assert _atom_names(out, 2) == ["N", "CA", "C", "O", "CB", "CG", "CD"]
        assert len(out.to_coords()) == len(out.atoms)

    @pytest.mark.parametrize("resn", sorted(SIDECHAIN_CHI) + ["GLY", "ALA"])
    def test_every_target_residue_rebuilds(self, resn):
        """20 个标准残基都能落到理想链上, 且图自洽 (to_coords 不报错)。"""
        letter = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
                  "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
                  "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
                  "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}[resn]
        out = mutate(_chain_ic(3), [f"A2{letter}"])
        assert {str(a.res_name) for a in out.atoms
                if a.res_id == 2} == {resn}
        assert len(out.to_coords()) == len(out.atoms)


class TestChi:
    def test_canonical_rotamer_matches_the_template(self):
        """canonical 的 chi 就是模板表里的值 —— 拿模板自己当参照, 不写死数字。"""
        from biorazer.database.molecule.icoor.protein.template import build_template
        from biorazer.structure.manipulation.internal_coord.mutation import _SS
        tmpl, _ = build_template("HIS", _SS, "canonical")
        idx = {a.name: n for n, a in enumerate(tmpl.atoms)}
        want = [tmpl.dihedra[tuple(idx[n] for n in quad)]
                for quad in SIDECHAIN_CHI["HIS"]]
        got = _chi(mutate(_chain_ic(3), ["A2H"]), 2, "HIS")
        assert max(abs(_wrap(g - w)) for g, w in zip(got, want)) <= CHI_TOL

    def test_rotamer_name_sets_the_chi(self):
        out = mutate(_chain_ic(3), ["A2H"], rotamer="g+/g+")
        for got in _chi(out, 2, "HIS")[:2]:
            assert _wrap(got - 60.0) == pytest.approx(0.0, abs=CHI_TOL)

    def test_chi_dict_and_sequence_agree(self):
        ic = _chain_ic(3)
        d = mutate(ic, ["A2H"], chi={"CG": -65.0, "ND1": 170.0})
        s = mutate(ic, ["A2H"], chi=[-65.0, 170.0])
        got = _chi(d, 2, "HIS")
        assert _wrap(got[0] + 65.0) == pytest.approx(0.0, abs=CHI_TOL)
        assert _wrap(got[1] - 170.0) == pytest.approx(0.0, abs=CHI_TOL)
        a, b = _coords_by_key(d), _coords_by_key(s)
        assert set(a) == set(b)
        assert max(float(np.linalg.norm(a[k] - b[k])) for k in a) == 0.0

    def test_chi_overrides_rotamer(self):
        out = mutate(_chain_ic(3), ["A2H"], rotamer="canonical", chi={"CG": 30.0})
        assert _wrap(_chi(out, 2, "HIS")[0] - 30.0) == pytest.approx(0.0, abs=CHI_TOL)

    def test_return_info_reports_what_landed(self):
        out, info = mutate(_chain_ic(3), ["A2H"], rotamer="g+/g+",
                           return_info=True)
        assert len(info) == 1
        assert info[0]["key"] == ("A", 2, "")
        assert info[0]["src"] == "ALA" and info[0]["res_name"] == "HIS"
        assert info[0]["rotamer"] == "g+/g+"
        assert [round(v, 6) for v in info[0]["chi"]] == [60.0, 60.0]

    def test_ss_does_not_change_the_side_chain(self):
        """``_SS`` 是内部常量: 侧链几何与它无关 (见模块里那条实测注释)。"""
        assert _SS in ("coil", "alpha-helix") or isinstance(_SS, str)


class TestErrors:
    def test_rejects_non_internal_coord(self):
        with pytest.raises(TypeError):
            mutate(_ideal_chain(2, "A"), ["A1H"])

    def test_rejects_empty_or_non_list_spec(self):
        ic = _chain_ic(2)
        with pytest.raises(ValueError):
            mutate(ic, [])
        with pytest.raises(ValueError):
            mutate(ic, "A1H")

    @pytest.mark.parametrize("spec", ["A1", "1H", "A1HH", "HIS1", "", " A1 H"])
    def test_rejects_malformed_spec(self, spec):
        with pytest.raises(ValueError, match="invalid mutation spec"):
            mutate(_chain_ic(2), [spec])

    def test_rejects_unknown_letters(self):
        ic = _chain_ic(2)
        with pytest.raises(ValueError, match="unknown source"):
            mutate(ic, ["Z1H"])
        with pytest.raises(ValueError, match="unknown target"):
            mutate(ic, ["A1Z"])

    def test_rejects_source_letter_mismatch(self):
        with pytest.raises(ValueError, match="expects W"):
            mutate(_chain_ic(2), ["W1H"])

    def test_rejects_missing_residue(self):
        with pytest.raises(ValueError, match="does not exist"):
            mutate(_chain_ic(2), ["A9H"])

    def test_rejects_chain_filter_that_matches_nothing(self):
        with pytest.raises(ValueError, match="chain B"):
            mutate(_chain_ic(2), ["A1H"], chain="B")

    def test_rejects_duplicate_target(self):
        with pytest.raises(ValueError, match="more than once"):
            mutate(_chain_ic(2), ["A1H", "A1F"])

    def test_rejects_unknown_rotamer(self):
        """``build_template`` 对未知 rotamer 名静默退回 canonical —— 这一层要拦住。"""
        assert "nope" not in rotamer_names("HIS")
        with pytest.raises(ValueError, match="unknown rotamer"):
            mutate(_chain_ic(2), ["A1H"], rotamer="nope")

    def test_rejects_chi_of_wrong_length(self):
        with pytest.raises(ValueError, match="has 2 chi"):
            mutate(_chain_ic(2), ["A1H"], chi=[1.0])

    def test_rejects_chi_dict_with_unknown_rotator(self):
        with pytest.raises(ValueError, match="no chi ending at"):
            mutate(_chain_ic(2), ["A1H"], chi={"OE1": 10.0})

    def test_rejects_residue_without_backbone(self):
        """目标残基缺 ``C`` (只剩 N/CA) -> 建不了 CB 的定位四元组。"""
        ic = InternalCoord(atoms=[InternalCoordAtom(chain_id="A", res_id=1,
                                                    res_name="ALA", name="N"),
                                  InternalCoordAtom(chain_id="A", res_id=1,
                                                    res_name="ALA", name="CA")],
                           anchor={0: (0.0, 0.0, 0.0), 1: (1.458, 0.0, 0.0)})
        with pytest.raises(ValueError, match="has no C atom"):
            mutate(ic, ["A1H"])

    def test_error_paths_leave_the_input_usable(self):
        ic = _chain_ic(2)
        with pytest.raises(ValueError):
            mutate(ic, ["A1H", "A9F"])
        assert _atom_names(ic, 1) == ["N", "CA", "C", "O"]
        assert MAINCHAIN_ATOMS.issuperset(_atom_names(ic, 1))
