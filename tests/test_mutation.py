# -*- coding: utf-8 -*-
"""Tests for :mod:`biorazer.structure.manipulation.mutation`.

被测的是"在真实骨架上按 rotamer 库**重建**侧链"这条链路:

* :func:`dihedral` 必须是仓库**唯一**的二面角实现 (re-export 自
  ``structure.objects.internal_coords``), 且与 :func:`_place` 严格互逆;
* :func:`ss_from_phi_psi` 的粗判;
* :func:`rotamer_candidates` 的概率降序与 (phi, psi) 箱;
* :func:`build_side_chain` 的原子集合与 chi 回量, 以及羰基 ``O`` 的定位
  (sp2 共面性: 真实 ``N_{i+1}`` (``next_n``) > ``psi`` > ``ss`` 类均值);
* :func:`_backbone_neighbours` 的肽键判据 (号相邻**且** ``C-N`` 距离成立);
* :func:`mutate` 的骨架不动、链作用域、各种非法输入。

夹具不依赖任何外部结构文件 (另一个仓库的 ``*.cif`` 不进测试), 也不用被测代码
造输入: 用仓库自己的 :func:`_place` 加上数据库的理想键长/键角
(:data:`...bond.length.protein.AMINO_ACID_BOND_LENGTH` /
:data:`...bond.angle.protein.AMINO_ACID_BACKBONE_BOND_ANGLE`) 长一条理想
alpha 螺旋骨架, 只含 ``N/CA/C/O`` —— 侧链留给被测的 :func:`mutate` 自己建。
"""

import numpy as np
import pytest
from biotite.structure import AtomArray

from biorazer.database.molecule.bond.angle.protein import (
    AMINO_ACID_BACKBONE_BOND_ANGLE as BANG,
)
from biorazer.database.molecule.bond.dihedral.protein import SIDECHAIN_CHI
from biorazer.database.molecule.bond.length.protein import (
    AMINO_ACID_BOND_LENGTH as BLEN,
)
from biorazer.database.molecule.icoor.protein import template
from biorazer.database.molecule.rotamer import pymol as PM
from biorazer.structure.manipulation import mutation
from biorazer.structure.manipulation.mutation import (
    build_side_chain,
    dihedral,
    mutate,
    rotamer_candidates,
    ss_from_phi_psi,
)
from biorazer.structure.objects import internal_coords as _ic
from biorazer.structure.selection.index.annotation import group_atoms_by_residue


# --------------------------------------------------------------------------
# fixtures / helpers
# --------------------------------------------------------------------------
PHI, PSI, OMEGA = -60.0, -45.0, 180.0

#: 从 ``AtomArray`` 回量 chi 时的容差 (度)。
#: ``biotite`` 的 ``coord`` 是 **float32** (相对精度 1e-7, 坐标量级 ~10 A),
#: 折成二面角后约 1e-5 度 —— 这是存储精度下限, 不是几何误差。实测
#: 6VY1 A31 HIS chi 回量偏差 1.3e-5 度。要求 1e-6 度等价于要求 float64 存储。
CHI_TOL = 1e-4


def _wrap(x):
    return (x + 180.0) % 360.0 - 180.0


def _mean(value):
    """数据库记录是 ``{mean, std, lb, up, source}`` (或裸 float) —— 取点值。"""
    return float(value["mean"]) if isinstance(value, dict) else float(value)


def _ideal_chain(n_res, chain_id, phi=PHI, psi=PSI, omega=OMEGA,
                 res_names=None):
    """长一条理想 alpha 螺旋骨架 (只含 ``N/CA/C/O``), 返回 ``AtomArray``。

    第 1 个残基的 ``N/CA/C`` 直接用模板的 anchor 帧 (模板的 anchor 就是理想
    键长/键角), 之后每个原子都由 :func:`_place` 从已放好的三个父原子生长 ——
    即仓库里 ``InternalCoord.to_coords`` 用的同一个原语。所以夹具的几何和
    数据库表严格一致, 且不 hardcode 任何数字。
    """
    if res_names is None:
        res_names = ["ALA"] * n_res
    assert len(res_names) == n_res

    ic = template.build_template("ALA", "alpha-helix")
    pos = {a.name: np.asarray(ic.anchor[i], float)
           for i, a in enumerate(ic.atoms) if i in ic.anchor}
    n_coord, ca_coord, c_coord = pos["N"], pos["CA"], pos["C"]

    l_cn, l_nca, l_cac = (_mean(BLEN[("C", "N")]), _mean(BLEN[("N", "CA")]),
                          _mean(BLEN[("CA", "C")]))
    l_co = _mean(BLEN[("C", "O")])
    a_ncac = _mean(BANG[("N", "CA", "C")])
    a_cacn = _mean(BANG[("CA", "C", "N")])
    a_cnca = _mean(BANG[("C", "N", "CA")])
    a_caco = _mean(BANG[("CA", "C", "O")])

    atoms = []          # (atom_name, res_id, xyz)
    names = ("N", "CA", "C", "O")
    for i in range(n_res):
        if i == 0:
            first = (n_coord, ca_coord, c_coord)
        else:
            pn, pca, pc = [a[2] for a in atoms[-4:-1]]     # 前一个残基的 N/CA/C
            first = (
                _ic._place(pn, pca, pc, l_cn, a_cacn, psi),
                None,
                None,
            )
            n_i = first[0]
            ca_i = _ic._place(pca, pc, n_i, l_nca, a_cnca, omega)
            c_i = _ic._place(pc, n_i, ca_i, l_cac, a_ncac, phi)
            first = (n_i, ca_i, c_i)
        n_i, ca_i, c_i = first
        o_i = _ic._place(n_i, ca_i, c_i, l_co, a_caco, psi - 180.0)
        for nm, xyz in zip(names, (n_i, ca_i, c_i, o_i)):
            atoms.append((nm, i + 1, np.asarray(xyz, float)))

    arr = AtomArray(len(atoms))
    arr.coord = np.array([a[2] for a in atoms], float)
    arr.atom_name = np.array([a[0] for a in atoms], dtype="U4")
    arr.res_id = np.array([a[1] for a in atoms], dtype=np.int32)
    arr.res_name = np.array([res_names[a[1] - 1] for a in atoms], dtype="U3")
    arr.chain_id = np.array([chain_id] * len(atoms), dtype="U4")
    arr.element = np.array([a[0][0] for a in atoms], dtype="U2")
    arr.ins_code = np.array([""] * len(atoms), dtype="U4")
    arr.hetero = np.array([False] * len(atoms))
    return arr


def _next_n_from(bb, psi):
    """在 ``bb`` 的 ``C`` 上按理想肽键几何放一个 ``N_{i+1}`` (绕 ``CA-C`` 到 ``psi``)。

    用被测代码外部的原语 :func:`_place` + 数据库的键长/键角, 所以不循环依赖
    ``build_side_chain``: 夹具给的是一个**真实存在**的下一个 N。
    """
    return _ic._place(bb["N"], bb["CA"], bb["C"], _mean(BLEN[("C", "N")]),
                      _mean(BANG[("CA", "C", "N")]), psi)


def _bb_dict(arr, chain, res_id):
    """把某残基的 ``N/CA/C`` 取成 ``build_side_chain`` 要的字典。"""
    m = (arr.chain_id == chain) & (arr.res_id == res_id)
    return {str(n): np.asarray(arr.coord[m & (arr.atom_name == n)][0], float)
            for n in ("N", "CA", "C")}


def _xyz(arr, name):
    """按原子名取坐标 (arr 只含一个残基时用)。"""
    return np.asarray(arr.coord[[str(n) for n in arr.atom_name].index(name)], float)


@pytest.fixture
def two_chain_array():
    """两条链 (A/B), 各 6 个 ALA 残基, 只有主链原子。"""
    a = _ideal_chain(6, "A")
    b = _ideal_chain(6, "B")
    from biotite.structure import concatenate
    return concatenate([a, b])


@pytest.fixture
def backbone_31():
    """6VY1 A31 (真实 alpha 螺旋位置) 的 ``{N, CA, C, O}`` —— 用于 build_side_chain。

    这个坐标是 6VY1 链 A 31 号残基的实测值; 用理想骨架自己的理想值替代也行,
    但真实坐标能顺带覆盖"骨架不完美"的情况。
    """
    return {
        "N": np.array([-1.48913, 8.35542, 1.02062]),
        "CA": np.array([-0.31083, 7.50055, 0.90623]),
        "C": np.array([0.74021, 8.13639, 0.00566]),
        "O": np.array([1.72553, 7.51708, -0.39438]),
    }


def _residue_atoms(arr, chain, res_id):
    m = (arr.chain_id == chain) & (arr.res_id == res_id)
    return m, [str(n) for n in arr.atom_name[m]]


def _measure_chi(arr, chain, res_id, res_name):
    """从 ``arr`` 里回量某残基的 chi 列表 (用被测的 :func:`dihedral`)。"""
    c = {}
    m = (arr.chain_id == chain) & (arr.res_id == res_id)
    for i in np.nonzero(m)[0]:
        c[str(arr.atom_name[i])] = np.asarray(arr.coord[i], float)
    return [dihedral(*[c[a] for a in quad]) for quad in SIDECHAIN_CHI[res_name]]


def _assert_chi(got, want, tol=CHI_TOL):
    """比 chi 要按**圆周**比: ``180.0`` 与 ``-180.0`` 是同一个角。"""
    assert len(got) == len(want), f"chi 个数不符: {got} vs {want}"
    worst = max(abs(_wrap(g - w)) for g, w in zip(got, want))
    assert worst <= tol, f"chi 偏差 {worst:.3e} deg 超过 {tol:.0e}: {got} vs {want}"


def _need_pymol():
    if not PM.PYMOL_SC_BB_DEP.exists():
        pytest.skip(f"PyMOL rotamer library not found: {PM.PYMOL_SC_BB_DEP}")


# --------------------------------------------------------------------------
# 1. dihedral -- 单一数据源
# --------------------------------------------------------------------------
class TestDihedral:
    def test_matches_internal_coords(self):
        rng = np.random.default_rng(11)
        for _ in range(200):
            p = rng.normal(size=(4, 3))
            assert dihedral(*p) == pytest.approx(_ic._dihedral(*p), abs=1e-12)

    def test_is_the_same_object_as_internal_coords(self):
        """不能有两套约定并存: mutation.dihedral 就是 internal_coords 那个。"""
        assert dihedral is _ic.dihedral
        assert _ic._dihedral is _ic.dihedral      # 旧私有名保留为别名

    def test_iupac_sign_against_ideal_helix(self):
        """alpha 螺旋的 phi/psi 必须是负值 (-60/-45 附近), 不是 +60/+45。"""
        arr = _ideal_chain(4, "A")
        phis = _measure_phi(arr)
        assert max(phis) < -30.0, f"phi 符号约定错了: {phis}"

    def test_round_trip_with_place(self):
        """dihedral 与 _place 必须严格互逆 (这是 build_side_chain 正确的前提)。"""
        p0 = np.array([1.0, 0.3, -0.7])
        p1 = np.array([0.0, 0.0, 0.0])
        p2 = np.array([0.0, 1.5, 0.2])
        for ang in (-180.0, -120.0, -60.0, 0.0, 60.0, 179.0):
            p3 = _ic._place(p0, p1, p2, 1.5, 110.0, ang)
            assert _wrap(dihedral(p0, p1, p2, p3) - ang) == pytest.approx(
                0.0, abs=1e-9)


def _measure_phi(arr, chain="A"):
    """整条链的 phi (按原子名取邻居, 与 mutate 内部同一套逻辑)。"""
    out = []
    ids = sorted({int(r) for r in arr.res_id[arr.chain_id == chain]})
    for k, rid in enumerate(ids):
        if k == 0 or ids[k - 1] != rid - 1:
            continue
        c = {}

        def get(name, r):
            m = ((arr.chain_id == chain) & (arr.res_id == r)
                 & (arr.atom_name == name))
            return np.asarray(arr.coord[m][0], float)
        out.append(dihedral(get("C", ids[k - 1]), get("N", rid),
                            get("CA", rid), get("C", rid)))
    return np.array(out)


# --------------------------------------------------------------------------
# 2. ss_from_phi_psi
# --------------------------------------------------------------------------
class TestSecondaryStructureGuess:
    @pytest.mark.parametrize("phi,psi,expected", [
        (-60.0, -45.0, "alpha-helix"),
        (-120.0, 130.0, "beta-strand"),
        (0.0, 0.0, "coil"),
    ])
    def test_representative_values(self, phi, psi, expected):
        assert ss_from_phi_psi(phi, psi) == expected

    def test_returns_a_template_ss_class(self):
        for phi, psi in ((-60, -45), (-120, 130), (0, 0), (-75, 145)):
            assert ss_from_phi_psi(phi, psi) in template.SS_CLASSES


# --------------------------------------------------------------------------
# 3. rotamer_candidates
# --------------------------------------------------------------------------
class TestRotamerCandidates:
    def test_probability_descending(self):
        _need_pymol()
        cands = rotamer_candidates("HIS", -65.0, -35.0, library="pymol")
        assert len(cands) > 1
        p = [r.probability for r in cands]
        assert all(p[i] >= p[i + 1] for i in range(len(p) - 1))

    def test_top_limits_and_chi_length(self):
        _need_pymol()
        cands = rotamer_candidates("HIS", -65.0, -35.0, library="pymol", top=2)
        assert len(cands) == 2
        assert all(len(r.chi) == len(SIDECHAIN_CHI["HIS"]) for r in cands)

    def test_unknown_residue_raises(self):
        _need_pymol()
        with pytest.raises(ValueError):
            rotamer_candidates("XXX", -60.0, -45.0)

    def test_unknown_library_raises(self):
        with pytest.raises(ValueError):
            rotamer_candidates("HIS", -60.0, -45.0, library="nope")


# --------------------------------------------------------------------------
# 4/5. build_side_chain
# --------------------------------------------------------------------------
class TestBuildSideChain:
    def test_ala_has_exactly_the_backbone_plus_cb(self):
        aa = build_side_chain("ALA", {
            "N": [0, 0, 0], "CA": [1.458, 0, 0], "C": [2.0095, 1.4218, 0]},
            ss="alpha-helix")
        assert [str(n) for n in aa.atom_name] == ["N", "CA", "C", "O", "CB"]

    def test_his_chi_round_trip(self, backbone_31):
        want = (-60.0, -60.0)
        aa = build_side_chain("HIS", backbone_31, chi=want, ss="alpha-helix")
        _assert_chi(_measure_chi(aa, "A", 1, "HIS"), want)

    @pytest.mark.parametrize("want", [(-60.0, 180.0, 0.0), (65.0, -70.0, 120.0)])
    def test_glu_chi_round_trip(self, backbone_31, want):
        aa = build_side_chain("GLU", backbone_31, chi=want, ss="alpha-helix")
        _assert_chi(_measure_chi(aa, "A", 1, "GLU"), want)

    def test_backbone_lands_on_the_given_coordinates(self, backbone_31):
        """``N/CA/C`` 必须落在传入的真实坐标上。

        ``O`` 不在此列: 它由模板按 ``psi`` 重新生长 (见 ``build_side_chain``
        的 Parameters), 传入的 ``O`` 会被忽略 —— ``mutate`` 最后是从
        ``remove_side_chains`` 的输出里保留原来的 ``O``, 不靠这里生长的那个。
        """
        aa = build_side_chain("HIS", backbone_31, chi=(-60, -60),
                              ss="alpha-helix")
        for i, nm in enumerate(aa.atom_name):
            nm = str(nm)
            if nm in ("N", "CA", "C"):
                # 输入是 float64, 存进 AtomArray 后被降为 float32
                assert np.allclose(aa.coord[i], backbone_31[nm], atol=1e-6)

    def test_unknown_residue_raises(self, backbone_31):
        with pytest.raises(ValueError):
            build_side_chain("XXX", backbone_31)

    def test_unknown_ss_raises(self, backbone_31):
        with pytest.raises(ValueError):
            build_side_chain("HIS", backbone_31, ss="not-a-class")

    def test_dict_chi_keyed_by_terminal_atom(self, backbone_31):
        aa = build_side_chain("HIS", backbone_31, chi={"CG": -70.0, "ND1": 90.0},
                              ss="alpha-helix")
        got = _measure_chi(aa, "A", 1, "HIS")
        assert got == pytest.approx([-70.0, 90.0], abs=CHI_TOL)

    def test_too_many_chi_values_raises(self, backbone_31):
        with pytest.raises(ValueError):
            build_side_chain("HIS", backbone_31, chi=(-60, -60, -60))

    # ---- 羰基 O: sp2 共面性, 三档信息 next_n > psi > ss 类均值 ----

    def test_carbonyl_o_is_anti_to_the_real_next_n(self, backbone_31):
        """(1) 给真实 ``N_{i+1}``: O 与它绕 ``CA-C`` 反平行, 且落在同一平面内。

        sp2 共面性只保证"反平行 + 共面", 不含任何键角数值 —— 所以这里同时查
        圆周关系 (差 180) 和出平面距离 (0)。角度容差用 ``CHI_TOL``、距离容差用
        ``1e-6`` A: 都是 ``AtomArray`` 的 float32 存储精度下限 (见 ``CHI_TOL``
        注释), 不是几何误差。
        """
        psi = -49.1
        n_next = _next_n_from(backbone_31, psi)
        aa = build_side_chain("ALA", backbone_31, ss="coil", next_n=n_next)
        got = dihedral(_xyz(aa, "N"), _xyz(aa, "CA"), _xyz(aa, "C"), _xyz(aa, "O"))
        assert _wrap(got - (psi - 180.0)) == pytest.approx(0.0, abs=CHI_TOL)
        nrm = np.cross(_xyz(aa, "CA") - _xyz(aa, "C"),
                       n_next - _xyz(aa, "C"))
        oob = abs(np.dot(_xyz(aa, "O") - _xyz(aa, "C"),
                         nrm / np.linalg.norm(nrm)))
        assert oob < 1e-6, f"O 不在 (CA, C, N_next) 平面内: 出平面 {oob:.2e} A"

    def test_next_n_overrides_the_ss_class(self, backbone_31):
        """(2) ``ss="coil"`` (均值 psi = 0 -> O 180) 也要被真实 ``N_{i+1}`` 覆盖。"""
        aa = build_side_chain("ALA", backbone_31, ss="coil",
                              next_n=_next_n_from(backbone_31, -45.0))
        got = dihedral(_xyz(aa, "N"), _xyz(aa, "CA"), _xyz(aa, "C"), _xyz(aa, "O"))
        assert _wrap(got - 135.0) == pytest.approx(0.0, abs=CHI_TOL)

    def test_explicit_psi_drives_the_carbonyl_o(self, backbone_31):
        """(3) 没有 ``N_{i+1}`` 时 ``psi`` 直接用; 两者都没有才退回 ``ss`` 均值。"""
        aa = build_side_chain("ALA", backbone_31, ss="coil", psi=-100.0)
        got = dihedral(_xyz(aa, "N"), _xyz(aa, "CA"), _xyz(aa, "C"), _xyz(aa, "O"))
        assert _wrap(got - 80.0) == pytest.approx(0.0, abs=CHI_TOL)
        fallback = build_side_chain("ALA", backbone_31, ss="coil")
        assert _wrap(dihedral(_xyz(fallback, "N"), _xyz(fallback, "CA"),
                              _xyz(fallback, "C"), _xyz(fallback, "O"))
                     - 180.0) == pytest.approx(0.0, abs=CHI_TOL)

    def test_o_reproduces_the_chain_own_carbonyl(self):
        """在夹具链上重建: 传真实 ``N_{i+1}`` 应复现夹具自己生长的那个 ``O``。"""
        arr = _ideal_chain(4, "A")
        bb = _bb_dict(arr, "A", 2)
        n_next = np.asarray(arr.coord[(arr.chain_id == "A") & (arr.res_id == 3)
                                      & (arr.atom_name == "N")][0], float)
        aa = build_side_chain("ALA", bb, ss="coil", next_n=n_next)
        o_fixture = np.asarray(arr.coord[(arr.chain_id == "A") & (arr.res_id == 2)
                                         & (arr.atom_name == "O")][0], float)
        assert np.allclose(_xyz(aa, "O"), o_fixture, atol=1e-4)

    def test_next_n_that_is_not_the_peptide_partner_raises(self, backbone_31):
        with pytest.raises(ValueError, match="peptide bond"):
            build_side_chain("ALA", backbone_31,
                             next_n=np.array([50.0, 50.0, 50.0]))

    def test_next_n_wrong_shape_raises(self, backbone_31):
        with pytest.raises(ValueError, match="3-vector"):
            build_side_chain("ALA", backbone_31, next_n=np.zeros((2, 3)))


# --------------------------------------------------------------------------
# 5b. 肽键邻居判据 (号相邻但肽键不成 -> 必须判为"没有邻居")
# --------------------------------------------------------------------------
class TestPeptideNeighbours:
    """``_backbone_neighbours`` 除了 ``res_id`` 相邻, 还要量 ``C-N`` 距离。

    这条判据不严的话, ``next_n`` 会把一个假肽平面直接喂给羰基 ``O`` 的定位
    (假 ``psi`` -> O 偏几 A), 比退回 ``ss`` 均值更糟。
    """

    @staticmethod
    def _broken_gap():
        """3 个 ALA 的链, 中间那个残基整体平移 20 A (号仍相邻)。"""
        arr = _ideal_chain(3, "A")
        arr.coord[arr.res_id == 2] += np.array([20.0, 0.0, 0.0])
        return arr

    def test_intact_chain_gives_both_neighbours(self):
        arr = _ideal_chain(3, "A")
        groups = group_atoms_by_residue(arr)
        prev_c, next_n = mutation._backbone_neighbours(arr, groups, ("A", 2, ""))
        assert prev_c is not None and next_n is not None

    def test_broken_peptide_bond_gives_none(self):
        arr = self._broken_gap()
        groups = group_atoms_by_residue(arr)
        assert mutation._backbone_neighbours(arr, groups, ("A", 2, "")) == (None, None)
        assert mutation._backbone_neighbours(arr, groups, ("A", 1, ""))[1] is None
        assert mutation._backbone_neighbours(arr, groups, ("A", 3, ""))[0] is None


# --------------------------------------------------------------------------
# 6-15. mutate
# --------------------------------------------------------------------------
class TestMutate:
    def test_backbone_unchanged(self, two_chain_array):
        out = mutate(two_chain_array, ["A3H"], chain="A")
        for nm in ("N", "CA", "C", "O"):
            for r in (3, 4):
                a = two_chain_array.coord[
                    (two_chain_array.chain_id == "A")
                    & (two_chain_array.res_id == r)
                    & (two_chain_array.atom_name == nm)]
                b = out.coord[(out.chain_id == "A") & (out.res_id == r)
                              & (out.atom_name == nm)]
                assert np.allclose(a, b, atol=1e-6), f"A{r} {nm} 动了"

    def test_ala_to_his_has_ten_atoms(self, two_chain_array):
        out = mutate(two_chain_array, ["A3H"], chain="A")
        m, names = _residue_atoms(out, "A", 3)
        assert str(out.res_name[m][0]) == "HIS"
        assert sorted(names) == sorted(
            ["N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"])

    def test_chain_scope(self, two_chain_array):
        out = mutate(two_chain_array, ["A3H"], chain="A")
        assert str(out.res_name[(out.chain_id == "A") & (out.res_id == 3)][0]) \
            == "HIS"
        assert str(out.res_name[(out.chain_id == "B") & (out.res_id == 3)][0]) \
            == "ALA"
        # 夹具的 ALA 是 4 原子 (主链), HIS 是 10 原子
        assert len(out) == len(two_chain_array) + (10 - 4)

    def test_no_chain_argument_hits_every_chain(self, two_chain_array):
        out = mutate(two_chain_array, ["A3H"])
        for ch in ("A", "B"):
            assert str(out.res_name[(out.chain_id == ch)
                                    & (out.res_id == 3)][0]) == "HIS"

    def test_source_letter_mismatch_raises(self, two_chain_array):
        with pytest.raises(ValueError, match="expects"):
            mutate(two_chain_array, ["K3H"], chain="A")

    @pytest.mark.parametrize("spec", ["31H", "H", "3", "A31", "A31HZ",
                                      "A-3H", "A31*", ""])
    def test_bad_spec_format_raises(self, two_chain_array, spec):
        with pytest.raises((ValueError, TypeError)):
            mutate(two_chain_array, [spec], chain="A")

    def test_missing_residue_raises(self, two_chain_array):
        with pytest.raises(ValueError, match="does not exist"):
            mutate(two_chain_array, ["A99H"], chain="A")

    def test_empty_spec_raises(self, two_chain_array):
        with pytest.raises(ValueError):
            mutate(two_chain_array, [])

    def test_duplicate_residue_raises(self, two_chain_array):
        with pytest.raises(ValueError, match="more than once"):
            mutate(two_chain_array, ["A3H", "A3F"], chain="A")

    def test_gly_to_x_builds_side_chain(self):
        """GLY 没有 CB, 也必须能建出完整侧链。"""
        arr = _ideal_chain(6, "A", res_names=["ALA", "ALA", "GLY"] + ["ALA"] * 3)
        m, names = _residue_atoms(arr, "A", 3)
        assert "CB" not in names
        out = mutate(arr, ["G3H"], chain="A")
        m, names = _residue_atoms(out, "A", 3)
        assert sorted(names) == sorted(
            ["N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"])

    def test_x_to_gly_drops_the_side_chain(self, two_chain_array):
        out = mutate(two_chain_array, ["A3G"], chain="A")
        m, names = _residue_atoms(out, "A", 3)
        assert str(out.res_name[m][0]) == "GLY"
        assert sorted(names) == ["C", "CA", "N", "O"]

    def test_rotamer_index_out_of_range_raises(self, two_chain_array):
        with pytest.raises(ValueError, match="out of range"):
            mutate(two_chain_array, ["A3H"], chain="A", rotamer=999)

    def test_explicit_chi_values(self, two_chain_array):
        want = (-70.0, 100.0)
        out = mutate(two_chain_array, ["A3H"], chain="A", rotamer=want)
        got = _measure_chi(out, "A", 3, "HIS")
        assert got == pytest.approx(list(want), abs=CHI_TOL)

    def test_return_info_matches_the_chosen_rotamer(self, two_chain_array):
        _need_pymol()
        out, info = mutate(two_chain_array, ["A3H"], chain="A", rotamer=0,
                           return_info=True)
        assert len(info) == 1
        d = info[0]
        assert d["key"] == ("A", 3, "")
        assert d["src"] == "ALA" and d["res_name"] == "HIS"
        cands = rotamer_candidates("HIS", d["phi"], d["psi"], library="pymol")
        assert d["chi"] == pytest.approx(tuple(cands[0].chi), abs=1e-9)
        assert d["probability"] == pytest.approx(cands[0].probability)
        assert d["n_candidates"] == len(cands)
        # 报出来的 chi 就是结构里量到的 chi
        assert _measure_chi(out, "A", 3, "HIS") == pytest.approx(
            list(d["chi"]), abs=CHI_TOL)

    def test_info_default_chi_is_the_best_rotamer(self, two_chain_array):
        _need_pymol()
        _, info = mutate(two_chain_array, ["A3H"], chain="A",
                         return_info=True)
        cands = rotamer_candidates("HIS", info[0]["phi"], info[0]["psi"],
                                   library="pymol")
        assert info[0]["chi"] == pytest.approx(tuple(cands[0].chi), abs=1e-9)

    def test_rotamer_index_selects_that_candidate(self, two_chain_array):
        _need_pymol()
        _, i0 = mutate(two_chain_array, ["A3H"], chain="A", rotamer=0,
                       return_info=True)
        _, i1 = mutate(two_chain_array, ["A3H"], chain="A", rotamer=1,
                       return_info=True)
        assert i0[0]["chi"] != i1[0]["chi"]

    def test_phi_psi_are_the_real_backbone_torsions(self, two_chain_array):
        """phi/psi 必须按原子名取邻居 C/N —— 早期版本取了上一个残基的最后一个
        原子 (侧链末端), phi 会静默错掉。这里用理想骨架核对。"""
        _, info = mutate(two_chain_array, ["A3H"], chain="A",
                         return_info=True)
        assert info[0]["phi"] == pytest.approx(PHI, abs=0.05)
        assert info[0]["psi"] == pytest.approx(PSI, abs=0.05)

    def test_multiple_mutations_at_once(self, two_chain_array):
        out = mutate(two_chain_array, ["A2H", "A5F"], chain="A")
        assert str(out.res_name[(out.chain_id == "A")
                                & (out.res_id == 2)][0]) == "HIS"
        assert str(out.res_name[(out.chain_id == "A")
                                & (out.res_id == 5)][0]) == "PHE"

    def test_ala_to_his_does_not_touch_other_residues(
            self, two_chain_array):
        out = mutate(two_chain_array, ["A3H"], chain="A")
        keep = ~((out.chain_id == "A") & (out.res_id == 3))
        src = ~((two_chain_array.chain_id == "A")
                & (two_chain_array.res_id == 3))
        assert np.allclose(out.coord[keep], two_chain_array.coord[src],
                           atol=1e-9)
        assert [str(n) for n in out.atom_name[keep]] == \
            [str(n) for n in two_chain_array.atom_name[src]]

    def test_proline_keeps_its_ring(self, two_chain_array):
        """PRO 是环状残基, 原子集合必须完整 (CD 不能丢)。"""
        _need_pymol()
        out = mutate(two_chain_array, ["A3P"], chain="A")
        m, names = _residue_atoms(out, "A", 3)
        assert sorted(names) == ["C", "CA", "CB", "CD", "CG", "N", "O"]
