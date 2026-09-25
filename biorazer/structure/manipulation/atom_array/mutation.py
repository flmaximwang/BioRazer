# -*- coding: utf-8 -*-
"""按 rotamer 库**重建**侧链的点突变 (不是从别的结构搬运侧链)。

操作对象: biotite ``AtomArray`` —— 读入它的一段**真实骨架坐标**, 输出一个目标残基
已带上重建侧链的 ``AtomArray`` (``InternalCoord`` 只在 :func:`build_side_chain`
里作为一次性模板出现)。

与同目录 :mod:`.modification` 的分工
─────────────────────────────────────
``modification`` 提供的是零件:

* :func:`~biorazer.structure.manipulation.atom_array.modification.remove_side_chains` --
  把某些残基削成只剩 ``N/CA/C/O`` 的 GLY 骨架;
* :func:`~biorazer.structure.manipulation.atom_array.modification.mutate_without_side_chains` --
  同上, 但把残基名改成目标类型;
* :func:`~biorazer.structure.manipulation.atom_array.modification.replace_side_chains` --
  **移植**另一个结构里已有侧链的坐标 (graft)。

本模块补上缺的那一环: **在真实骨架上, 按 rotamer 库给出的连续 chi 值把侧链重建出来**。
区别很关键 —— graft 要求手头已经有一个构象合适的同种侧链, 而设计场景下往往没有;
重建只需要一个骨架 + 一个 rotamer 库。

算法
────
对每个要突变的残基 ``r``:

1. 取 ``r`` 的真实骨架 ``N/CA/C`` 坐标;
2. 从 rotamer 库的 (phi, psi) 箱里取候选 rotamer (连续 chi, 按概率降序);
3. 用 :func:`~biorazer.database.molecule.icoor.protein.template.build_template`
   取该残基的**理想 icoor 模板** (键长/键角/非 chi 二面角全部来自数据库表);
4. **把模板的 anchor 换成 ``r`` 的真实 ``N/CA/C``**, 并同步 anchor 的键长与键角 ——
   这一步是全部关键: 模板自带的 anchor 是"N 在原点、CA 沿 +x"的理想坐标系,
   直接用会把侧链丢在错误的地方;
5. 把 chi 四元组的二面角覆盖成候选 rotamer 的连续值;
6. :meth:`InternalCoord.to_coords` 从 anchor 生长出全部原子 —— 侧链就长在真实骨架上了。

因为 rotamer 是一次**刚性子树旋转**, 键长/键角/非 chi 二面角在旋转下不变,
所以第 3 步取理想值、第 5 步只改 chi 是**严格正确**的, 不是近似。

羰基 ``O`` (与侧链无关, 但同一函数要顺手给出)
───────────────────────────────────────────
侧链**不**依赖 ``O``: :data:`...topology.IC_PATH` 里没有任何含 ``O`` 的 quad,
``CB`` 的参考帧是本残基的 ``C/N/CA``。但 :func:`build_side_chain` 返回的残基要带
``O``。``O`` 绕 ``CA-C`` 的那个二面角**不是自由度**: 羰基碳是 sp2, 三个取代基
共面且该平面包含旋转轴 ``CA-C``, 因此
``dihedral(N, CA, C, O) = psi - 180`` (与键角数值无关, 见
:func:`~biorazer.database.molecule.icoor.protein.topology.carbonyl_o_dihedral`)。
放 ``O`` 的信息有三档, 从精确到近似:

1. 真实的下一个 N (``next_n``) -> 直接量出 ``psi``, **精确**;
2. 该残基的 ``psi`` -> 用它 (``psi`` 就是未知 ``N_{i+1}`` 的方位);
3. 都没有 -> 所属 ``ss`` 类的均值 ``psi`` (模板的默认值, 精确到该类的散布)。

下游再交给 :func:`~biorazer.structure.manipulation.atom_array.modification.replace_side_chains`
装配回整结构 (它保留原结构的 ``N/CA/C/O/OXT``, 所以 ``mutate`` 的输出用的是**原来
的** ``O``; 这套定位服务于直接调用 :func:`build_side_chain` 的场景)。

二面角约定 (单一数据源)
────────────────────────
本模块的 :func:`dihedral` **不是自己实现的**, 而是直接 re-export
``biorazer.structure.objects.internal_coords.dihedral`` —— 定义在
``biorazer/structure/objects/internal_coords.py:101``, 该文件把原先私有的
``_dihedral`` 提升为公开名并在同文件 ``:125`` 保留 ``_dihedral`` 别名。
那里的实现是
``b0 = -(p1-p0)`` / ``b1 = p2-p1`` / ``b2 = p3-p2`` 后取
``v = b0 - (b0·b1̂)b1̂``、``w = b2 - (b2·b1̂)b1̂``、
``atan2((b1×v)·w, v·w)`` —— **IUPAC 符号约定**, 与 Dunbrack/Rosetta 的 chi
值同号, 也与 :func:`biorazer.structure.objects.internal_coords._place` 的框架
严格互逆。二面角是全仓共用的基础量, 任何模块都不得再带一份私有副本。

一个实测过的坑 (2026-09, 5 分钟就能踩): 另一种常见写法
``n1 = b0×b1``、``m = n1×b1̂``、``atan2((m·n2), (n1·n2))`` (``n2 = b1×b2``)
看着也"合理", 但它给出的是 ``θ_false = -(θ_true + 180)``。在 6VY1 上量一条
真实的 LEU chi1 得到 ``-116.62`` 而真值是 ``-63.38``; 用它回量自己刚建好的
HIS 侧链 (要求 ``-60``) 会得到 ``-120``, 于是被误判成"建侧链的代码有 bug"。
判据: 同一残基内 chi1/chi2 的偏差**完全相等**、不同残基的偏差不同, 就是这种
公式错而不是随机几何错。改用 :func:`dihedral` 后回量偏差 0.0000 deg。

已知取舍 (实测值见 :func:`build_side_chain` 的 Notes)
──────────────────────────────────────────────────────
``CB`` 由骨架帧 ``(C, N, CA)`` + 理想二面角生长得到。真实结构里的 ``CB`` 会略有
差异, 因此重建的 ``CB`` 与野生型 ``CB`` 不严格重合 (量级 <0.1 A)。这对侧链构象
无影响, 但**不要把重建结果当作野生型的精确复现**。
"""

from __future__ import annotations

import numpy as np
import biotite.structure as bio_struct

from biorazer.database.molecule.bond.dihedral.protein import SIDECHAIN_CHI
from biorazer.database.molecule.bond.length.protein import AMINO_ACID_BOND_LENGTH
from biorazer.database.molecule.icoor.protein import template
from biorazer.database.molecule.icoor.protein.topology import carbonyl_o_dihedral
from biorazer.database.alphabet import AMINO_ACIDS_1TO3_UPPER

from biorazer.structure.objects import AtomArray
from ...objects.internal_coords import dihedral
from ...selection.index.annotation import group_atoms_by_residue
from .modification import remove_side_chains, replace_side_chains

__all__ = [
    "dihedral",
    "ss_from_phi_psi",
    "rotamer_candidates",
    "build_side_chain",
    "mutate",
]

#: 肽键 ``C-N`` 理想键长 (**A**)。判"这真的是相邻残基的 N/C 吗"用的期望值。
_PEPTIDE_CN = float(AMINO_ACID_BOND_LENGTH[("C", "N")]["mean"])
#: 肽键 ``C-N`` 的键长上界 (**A**)。超过它就不是肽键 —— 链断或传错残基。
#: 与 ``AtomArray_InternalCoord`` (``structure.bridge``) 的 peptide 环节用**同一个**判据, 不另立阈值。
_MAX_PEPTIDE_CN = float(AMINO_ACID_BOND_LENGTH[("C", "N")]["up"])

#: rotamer 库来源 -> 读取函数名 (惰性 import, 免得只为读一个库就把重依赖拉起来)
_LIBRARY_READERS = {
    "pymol-ind": ("biorazer.database.molecule.rotamer.pymol", "read_pymol_ind"),
    "pymol": ("biorazer.database.molecule.rotamer.pymol", "read_pymol_dep"),
    "pymol-dep": ("biorazer.database.molecule.rotamer.pymol", "read_pymol_dep"),
    "rosetta": ("biorazer.database.molecule.rotamer.rosetta", "read_bbdep02"),
    "shapovalov": ("biorazer.database.molecule.rotamer.rosetta", "read_shapovalov"),
}

_LIB_CACHE: dict[str, dict] = {}


def _load_library(library: str) -> dict:
    """按名字加载 rotamer 库, 结果缓存 (同一个库只读一次, pickle 不小)。"""
    key = library.lower()
    if key in _LIB_CACHE:
        return _LIB_CACHE[key]
    if key not in _LIBRARY_READERS:
        raise ValueError(
            f"unknown rotamer library {library!r}; "
            f"available: {sorted(_LIBRARY_READERS)}")
    mod_name, fn_name = _LIBRARY_READERS[key]
    import importlib
    fn = getattr(importlib.import_module(mod_name), fn_name)
    _LIB_CACHE[key] = fn()
    return _LIB_CACHE[key]


# --------------------------------------------------------------------------- #
# 小工具
# --------------------------------------------------------------------------- #
# ``dihedral`` 不是本模块自己实现的, 而是从
# :mod:`biorazer.structure.objects.internal_coords` re-export —— 见模块 docstring
# 的"二面角约定"一节。


def ss_from_phi_psi(phi: float, psi: float, tol: float = 30.0) -> str:
    """由 (phi, psi) 粗判二级结构类别, 返回 :data:`template.SS_CLASSES` 里的键。

    只用于给模板选一个 ``ss``; 侧链几何本身**不依赖** ``ss`` (见
    :func:`build_side_chain` 的 Notes), 所以这里只需要粗判。

    Examples
    --------
    >>> ss_from_phi_psi(-60, -45)
    'alpha-helix'
    >>> ss_from_phi_psi(-120, 130)
    'beta-strand'
    >>> ss_from_phi_psi(80, -170)
    'polyproline-II'
    >>> ss_from_phi_psi(0, 0)
    'coil'
    """
    def near(a, b):
        d = abs((a - b + 180.0) % 360.0 - 180.0)
        return d <= tol

    if near(phi, -60) and near(psi, -45):
        return "alpha-helix"
    if near(phi, -120) and (near(psi, 130) or near(psi, 120)):
        return "beta-strand"
    if near(phi, -75) and near(psi, 145):
        return "polyproline-II"
    if near(phi, -60) and near(psi, -30):
        return "3-10-helix"
    return "coil"


def _residue_backbone_xyz(arr, idxs, key):
    """从残基原子下标里取出 ``{atom_name: xyz}``; 缺 N/CA/C 时报错。"""
    got = {}
    names = arr.atom_name
    for i in idxs:
        got[str(names[i])] = np.asarray(arr.coord[i], float)
    missing = [n for n in ("N", "CA", "C") if n not in got]
    if missing:
        raise ValueError(
            f"residue {key} lacks backbone atom(s) {missing}; "
            "cannot rebuild a side chain on it")
    return got


def _atom_xyz(arr, idxs, name):
    """在 ``idxs`` 指定的残基里按**原子名**找坐标; 找不到返回 ``None``。"""
    for i in idxs:
        if str(arr.atom_name[i]) == name:
            return np.asarray(arr.coord[i], float)
    return None


def _backbone_neighbours(arr, groups, key):
    """取 ``key`` 的肽键邻居坐标, 返回 ``(prev_C, next_N)``, 缺则 ``None``。

    邻接按**残基在数组里的顺序**判定 (与 :func:`group_atoms_by_residue` 一致),
    并要求 ``res_id`` 不跳号 (``0 <= res_id - prev_res_id <= 1``; 差 0 即插入码
    残基)。缺残基 / 断链时返回 ``None`` —— 跨断口算出来的 phi/psi 是假的,
    比返回 NaN 更危险。

    号相邻**还不够**: 重编号过的链、``TER`` 断口、或传错的残基都可能号相邻而
    肽键不成立。所以这里再量一次 ``C-N`` 距离, 超过键长表上界
    (``_MAX_PEPTIDE_CN``, 与 ``AtomArray_InternalCoord`` 同一判据) 就当作
    没有邻居。这条现在是硬要求: ``next_n`` 会被 :func:`build_side_chain` 直接
    用来**精确**定位羰基 ``O``, 假的邻居会给出假的肽平面, 比退回 ``ss`` 均值更糟。

    Notes
    -----
    这里**必须按原子名取** ``C`` 与 ``N``。早期版本用 ``groups[prev][-1]`` 当作
    "前一个残基的 C", 对 6VY1 A30 实测取到的是 ``OD2`` (侧链末端), 于是 phi
    从真值 ``-65.4`` 变成 ``-120.3`` —— 而且这个错误 phi 只是拿去查 rotamer 的
    (phi, psi) 箱, 不会报错, 只会悄悄选错 rotamer。
    """
    order = [k for k in groups if k[0] == key[0]]
    i = order.index(key)
    own = groups[key]
    n_self = _atom_xyz(arr, own, "N")
    c_self = _atom_xyz(arr, own, "C")
    prev_c = next_n = None
    if i > 0 and 0 <= key[1] - order[i - 1][1] <= 1:
        prev_c = _atom_xyz(arr, groups[order[i - 1]], "C")
        if (prev_c is not None and n_self is not None
                and float(np.linalg.norm(n_self - prev_c)) > _MAX_PEPTIDE_CN):
            prev_c = None
    if i + 1 < len(order) and 0 <= order[i + 1][1] - key[1] <= 1:
        next_n = _atom_xyz(arr, groups[order[i + 1]], "N")
        if (next_n is not None and c_self is not None
                and float(np.linalg.norm(next_n - c_self)) > _MAX_PEPTIDE_CN):
            next_n = None
    return prev_c, next_n


# --------------------------------------------------------------------------- #
# 1. rotamer 候选
# --------------------------------------------------------------------------- #
def rotamer_candidates(
    res_name: str,
    phi: float | None = None,
    psi: float | None = None,
    library: str = "pymol",
    top: int | None = None,
):
    """取某残基在某 (phi, psi) 箱内按概率降序的候选 rotamer (连续 chi 值)。

    Parameters
    ----------
    res_name : 三字母残基名 (大写, 如 ``"HIS"``)
    phi, psi : 主链二面角 (度)。给 ``None`` 则用 backbone-independent 行为 ——
        库是 backbone-dependent 时会退回"全局按概率降序"。
    library : ``"pymol"`` (默认, backbone-dependent) / ``"pymol-ind"`` /
        ``"rosetta"`` (Dunbrack 2002) / ``"shapovalov"`` (2010)
    top : 只取前 ``n`` 条; ``None`` = 全部

    Returns
    -------
    list[RotamerRecord]
        每条记录带 ``chi`` (chi1..chiN 顺序, 度) 与 ``probability``。

    Examples
    --------
    >>> recs = rotamer_candidates("HIS", -60, -45, top=1)   # doctest: +SKIP
    >>> recs[0].chi                                          # doctest: +SKIP
    (-64.5, -73.3)
    """
    res_name = res_name.upper()
    libs = _load_library(library)
    lib = libs.get(res_name)
    if lib is None:
        raise ValueError(
            f"{res_name} not present in rotamer library {library!r}; "
            f"available: {sorted(libs)}")
    if phi is not None and psi is not None:
        recs = lib.for_backbone(float(phi), float(psi))
    else:
        recs = sorted(lib.records, key=lambda r: -r.probability)
    return recs[:top] if top else recs


# --------------------------------------------------------------------------- #
# 2. 在真实骨架上重建一个残基
# --------------------------------------------------------------------------- #
def build_side_chain(
    res_name: str,
    backbone_xyz: dict[str, np.ndarray],
    chi=None,
    ss: str | None = None,
    phi: float | None = None,
    psi: float | None = None,
    next_n=None,
    tol: float = 1e-6,
) -> AtomArray:
    """按理想 icoor 模板 + 指定 chi, 在**真实骨架**上重建一个残基。

    Parameters
    ----------
    res_name : 目标三字母残基名 (大写)
    backbone_xyz : ``{"N": xyz, "CA": xyz, "C": xyz}`` —— 真实骨架坐标。
        额外的 ``O`` 会被忽略: ``O`` 一律按下面的三档信息重新定位。
    chi : chi 值, 单位度。可以是 ``None`` (用模板 canonical 值)、序列
        ``(chi1, chi2, ...)``、或 ``{原子名: 角度}`` (键是 chi 四元组的**末端**
        原子名, 如 ``{"CG": -60.0, "ND1": 90.0}``)。
    ss : 二级结构类别; ``None`` 时由 ``phi``/``psi`` 粗判 (见
        :func:`ss_from_phi_psi`)。**侧链几何不依赖 ss**, 它只决定骨架 ``O``
        的位置 —— 而且只在既没有 ``next_n`` 也没有 ``psi`` 时才轮得到它。
    phi, psi : 三级结构二面角 (度)。``psi`` 决定羰基 ``O`` 的二面角
        (``psi - 180``, 见 Notes); ``phi`` 只参与粗判 ``ss``。
    next_n : array_like (3,), optional
        下一个残基 (同链 ``i+1``) 的酰胺 ``N`` 的**真实坐标**。给了它, 羰基
        ``O`` 就按 C 的 sp2 共面性**精确**定位 (O 与 ``N_{i+1}`` 绕 ``CA-C``
        反平行), 不再需要 ``psi``/``ss``。若它与本残基 ``C`` 的距离超出肽键键长
        上界 (说明不是下一个残基的 N), 报 :exc:`ValueError` —— 而不是悄悄给出
        一个假肽平面。

    Returns
    -------
    biotite.structure.AtomArray
        该残基的全部重原子 (``N, CA, C, O, CB, <侧链>``); ``N/CA/C`` 落在传入的
        真实骨架上, 其余由模板生长。``chain_id``/``res_id`` 为占位值, 由调用方改写。

    Notes
    -----
    **为什么这样是对的不是近似的**: rotamer 是一次刚性子树旋转 —— 侧链绕着一条
    chi 键整体转动。键长、键角、以及任何"不在被转动子树内部"的二面角在旋转下
    都不变。所以模板里除 chi 以外的量可以直接用数据库的理想值, 只有 chi 需要
    换成目标值。

    **羰基 ``O`` 为什么需要 ``psi`` (或 ``next_n``)**: ``O`` 绕 ``CA-C`` 的
    二面角**不是**自由参数 —— 羰基碳 sp2 使 ``CA/O/N_{i+1}`` 共面, 而该平面包含
    旋转轴 ``CA-C``, 于是 ``dihedral(N, CA, C, O) = psi - 180`` 与键角数值无关
    (推导与实测残差见
    :func:`~biorazer.database.molecule.icoor.protein.topology.carbonyl_o_dihedral`)。
    但 ``psi = (N, CA, C, N_{i+1})`` 需要的 ``N_{i+1}`` 不在本残基的
    ``{N, CA, C}`` 里, 所以这里按 ``next_n`` > ``psi`` > ``ss`` 类均值三级取值
    (第一级精确, 后两级是退路)。实测 (6VY1 A10, 真 ``psi = -49.1``, 真 O 二面角
    ``131.4``):

    * 传 ``next_n`` -> ``130.90`` (真 O 是 ``131.44``, 差 ``-0.54`` = 该残基实际的
      出平面量; O 位移 ``0.025`` A) —— 与直接传 ``psi`` 同值;
    * 只给 ``ss="alpha-helix"`` -> ``135.00`` (差 ``+3.56``, O 偏 ``0.070`` A);
    * 只给 ``ss="coil"`` (均值 psi = 0) -> ``180.00`` (差 ``+48.56``, O 偏 ``0.876`` A);
    * 误判成 ``ss="beta-strand"`` -> ``-50.00`` (差 ``178.6``, O 偏 ``2.13`` A)。

    **已知取舍**: ``CB`` 是从骨架帧 ``(C, N, CA)`` 用理想二面角生长的, 与真实
    结构中观测到的 ``CB`` 会有细微差别。若要精确复现野生型, 应当把野生型的
    ``CB`` 直接传进来 (本函数不接受 ``CB``, 因为 GLY 没有 ``CB``, 统一用理想值
    才能一致处理"带 CB 的"和"不带 CB 的"突变起点)。

    Examples
    --------
    >>> bb = {"N": [0, 0, 0], "CA": [1.458, 0, 0], "C": [2.0095, 1.4218, 0]}
    >>> aa = build_side_chain("ALA", bb, ss="alpha-helix")
    >>> [str(n) for n in aa.atom_name]
    ['N', 'CA', 'C', 'O', 'CB']
    """
    res_name = res_name.upper()
    if res_name not in template.RESIDUES:
        raise ValueError(
            f"unknown residue {res_name!r}; expected one of {template.RESIDUES}")

    if ss is None:
        if phi is not None and psi is not None:
            ss = ss_from_phi_psi(phi, psi)
        else:
            ss = "coil"
    if ss not in template.SS_CLASSES:
        raise ValueError(
            f"unknown secondary-structure class {ss!r}; "
            f"expected one of {template.SS_CLASSES}")

    ic, _ = template.build_template(res_name, ss, "canonical")
    idx = {a.name: n for n, a in enumerate(ic.atoms)}

    # ---- 关键一步: 把 anchor 换成真实骨架, 并同步 anchor 内部的键长/键角 ----
    n_x = np.asarray(backbone_xyz["N"], float)
    ca_x = np.asarray(backbone_xyz["CA"], float)
    c_x = np.asarray(backbone_xyz["C"], float)
    ic.anchor[idx["N"]] = tuple(n_x)
    ic.anchor[idx["CA"]] = tuple(ca_x)
    ic.anchor[idx["C"]] = tuple(c_x)
    ic.bond_distances[(idx["N"], idx["CA"])] = float(np.linalg.norm(n_x - ca_x))
    ic.bond_distances[(idx["CA"], idx["C"])] = float(np.linalg.norm(ca_x - c_x))
    ang_acc = ca_x - n_x
    ang_acc = ang_acc / np.linalg.norm(ang_acc)
    ang_ccc = c_x - ca_x
    ang_ccc = ang_ccc / np.linalg.norm(ang_ccc)
    ic.bond_angles[(idx["N"], idx["CA"], idx["C"])] = float(
        np.degrees(np.arccos(np.clip(np.dot(ang_acc, ang_ccc), -1.0, 1.0))))

    # ---- 覆盖 chi ----
    if chi is not None:
        chis = SIDECHAIN_CHI.get(res_name, [])
        if isinstance(chi, dict):
            targets = {str(k): float(v) for k, v in chi.items()}
            for quad in chis:
                last = quad[3]
                if last in targets:
                    ic.dihedra[(idx[quad[0]], idx[quad[1]], idx[quad[2]],
                                idx[quad[3]])] = targets[last]
        else:
            vals = list(chi)
            if len(vals) > len(chis):
                raise ValueError(
                    f"{res_name} has {len(chis)} chi but {len(vals)} values given")
            for quad, val in zip(chis, vals):
                ic.dihedra[(idx[quad[0]], idx[quad[1]], idx[quad[2]],
                            idx[quad[3]])] = float(val)

    # ---- 羰基 O: 由 C 的 sp2 共面性定位 (三档信息, 从精确到近似) ----
    # 见 topology.carbonyl_o_dihedral: C 是 sp2 中心, CA/O/N_{i+1} 共面且该平面
    # 含旋转轴 CA-C, 所以 dihedral(N, CA, C, O) = psi - 180, 与键角数值无关。
    # 模板已按 ss 类的均值 psi 放好了 O (build_template), 这里用更硬的信息覆盖:
    #   1. 真实 N_{i+1} -> 量出 psi, 精确;
    #   2. 传入的 psi   -> 用它 (psi 就是未知 N_{i+1} 的方位);
    #   3. 两者都没有   -> 保留模板的值 (ss 类均值)。
    if next_n is not None:
        n2 = np.asarray(next_n, float).reshape(-1)
        if n2.shape != (3,):
            raise ValueError(f"next_n must be a 3-vector, got shape {n2.shape}")
        d_cn = float(np.linalg.norm(n2 - c_x))
        if d_cn > _MAX_PEPTIDE_CN:
            raise ValueError(
                f"next_n is {d_cn:.2f} A from this residue's C; a peptide bond "
                f"C-N is {_PEPTIDE_CN:.3f} A (max {_MAX_PEPTIDE_CN:.3f}), so this "
                "is not the next residue's amide N")
        psi = dihedral(n_x, ca_x, c_x, n2)
    if psi is not None:
        ic.dihedra[(idx["N"], idx["CA"], idx["C"], idx["O"])] = \
            carbonyl_o_dihedral(psi)

    coords = ic.to_coords(tol=tol)
    n_atoms = len(ic.atoms)
    aa = AtomArray(n_atoms)
    aa.coord = np.array([coords[i] for i in range(n_atoms)], float)
    aa.atom_name = np.array([a.name for a in ic.atoms], dtype="U4")
    aa.res_name = np.array([res_name] * n_atoms, dtype="U3")
    aa.element = np.array([a.element for a in ic.atoms])
    aa.chain_id = np.array(["A"] * n_atoms, dtype="U4")
    aa.res_id = np.array([1] * n_atoms, dtype=np.int32)
    aa.ins_code = np.array([""] * n_atoms, dtype="U4")
    aa.hetero = np.array([False] * n_atoms)
    return aa


# --------------------------------------------------------------------------- #
# 3. 完整突变
# --------------------------------------------------------------------------- #
def mutate(
    atom_array: AtomArray,
    mutation_spec,
    *,
    chain: str | None = None,
    library: str = "pymol",
    rotamer="best",
    return_info: bool = False,
):
    """把残基突变到目标类型, 并用 rotamer 库**重建**其侧链。

    与 :func:`~biorazer.structure.manipulation.atom_array.modification.mutate_without_side_chains`
    的区别: 那个只把侧链砍掉、改个名字, 得到的是**没有侧链**的骨架;
    本函数接着按 rotamer 库把侧链**建出来**。

    Parameters
    ----------
    atom_array : 输入结构 (biotite ``AtomArray``)
    mutation_spec : 突变列表, 每项形如 ``"A31H"`` ——
        目标残基的单字母代码 + 残基号 + 原残基的单字母代码? **不是**。
        本项目沿用 :func:`mutate_without_side_chains` 的既有约定:
        ``"<原残基字母><残基号><目标残基字母>"``, 例如 ``"A31H"`` 表示
        "31 位上的 ALA 变成 HIS"。原残基字母会与结构里的实际类型核对。
    chain : 只对该链做。``None`` = 所有链 (与既有函数一致)。
        注意既有函数按 ``res_id`` 匹配**全部链**, 多链同号时会一起改;
        要只改一条链就必须显式给 ``chain``。
    library : rotamer 库名字, 见 :func:`rotamer_candidates`
    rotamer : 用哪个 rotamer。``"best"`` = 该 (phi, psi) 箱里概率最高的;
        整数 = :func:`rotamer_candidates` 返回列表里的下标;
        序列 = 直接给 chi 值; 字典 = ``{chi 末端原子名: 角度}``。
    return_info : 是否连带返回每个残基实际用了哪个 rotamer

    Returns
    -------
    AtomArray, 或 ``(AtomArray, list[dict])`` (``return_info=True``)。
    ``info`` 每项含 ``key``/``res_name``/``phi``/``psi``/``ss``/``chi``/
    ``probability``/``n_candidates``。

    Raises
    ------
    ValueError
        规格格式错误、残基不存在、原残基字母对不上、缺少骨架原子等。

    Examples
    --------
    >>> out = mutate(arr, ["A31H"], chain="A")           # doctest: +SKIP
    >>> out, info = mutate(arr, ["A20Y"], return_info=True)   # doctest: +SKIP
    >>> info[0]["chi"]                                    # doctest: +SKIP
    (-64.5, -73.3)
    """
    import re

    if not isinstance(mutation_spec, (list, tuple)) or len(mutation_spec) == 0:
        raise ValueError("mutation_spec must be a non-empty list of strings")

    # ---- 解析规格 ----
    pattern = re.compile(r"^([A-Za-z])(\d+)([A-Za-z])$")
    wanted: dict[tuple[str, int], tuple[str, str]] = {}
    for entry in mutation_spec:
        if not isinstance(entry, str):
            raise TypeError(f"mutation_spec entries must be str, got {type(entry).__name__}")
        m = pattern.fullmatch(entry.strip())
        if m is None:
            raise ValueError(
                f"invalid mutation spec {entry!r}: expected 'A31H' "
                "(source letter, residue id, target letter)")
        src, res_id, tgt = m.group(1).upper(), int(m.group(2)), m.group(3).upper()
        for letter, which in ((src, "source"), (tgt, "target")):
            if letter not in AMINO_ACIDS_1TO3_UPPER:
                raise ValueError(f"unknown {which} residue letter {letter!r} in {entry!r}")
        key = (chain, res_id) if chain is not None else (None, res_id)
        if key in wanted:
            raise ValueError(f"residue {res_id} targeted more than once in mutation_spec")
        wanted[key] = (src, AMINO_ACIDS_1TO3_UPPER[tgt])

    groups = group_atoms_by_residue(atom_array)

    # ---- 定位每个目标残基 ----
    targets: dict[tuple, dict] = {}
    for (want_chain, res_id), (src, tgt_res_name) in wanted.items():
        keys = [k for k in groups
                if k[1] == res_id and (want_chain is None or k[0] == want_chain)]
        if not keys:
            raise ValueError(
                f"mutation_spec references residue {res_id}"
                + (f" of chain {want_chain}" if want_chain else "")
                + ", which does not exist in the input structure")
        for k in keys:
            idxs = groups[k]
            cur = str(atom_array.res_name[idxs[0]])
            expected = AMINO_ACIDS_1TO3_UPPER[src]
            if cur != expected:
                raise ValueError(
                    f"residue {k} is {cur}, but the mutation spec expects "
                    f"{src} ({expected})")
            bb = _residue_backbone_xyz(atom_array, idxs, k)
            prev_c, next_n = _backbone_neighbours(atom_array, groups, k)
            phi = psi = None
            if prev_c is not None:
                phi = dihedral(prev_c, bb["N"], bb["CA"], bb["C"])
            if next_n is not None:
                psi = dihedral(bb["N"], bb["CA"], bb["C"], next_n)
            targets[k] = dict(res_name=tgt_res_name, backbone=bb,
                              phi=phi, psi=psi, next_n=next_n, idxs=idxs)

    # ---- 逐残基选 rotamer 并重建 ----
    implants, info = [], []
    for k, t in targets.items():
        res_name = t["res_name"]
        phi, psi = t["phi"], t["psi"]
        ss = ss_from_phi_psi(phi, psi) if (phi is not None and psi is not None) else "coil"

        chis, proba, n_cand = None, float("nan"), 0
        if res_name not in ("GLY", "ALA"):
            cands = rotamer_candidates(res_name, phi, psi, library=library)
            n_cand = len(cands)
            if not cands:
                raise ValueError(
                    f"no rotamer candidates for {res_name} at phi={phi}, psi={psi} "
                    f"in library {library!r}")
            if rotamer == "best":
                rec = cands[0]
                chis, proba = rec.chi, rec.probability
            elif isinstance(rotamer, int):
                if not 0 <= rotamer < len(cands):
                    raise ValueError(
                        f"rotamer index {rotamer} out of range "
                        f"(0..{len(cands) - 1})")
                rec = cands[rotamer]
                chis, proba = rec.chi, rec.probability
            elif isinstance(rotamer, dict):
                chis = rotamer
            else:
                chis = rotamer
            if chis is None:
                chis = ()

        aa = build_side_chain(res_name, t["backbone"], chi=chis, ss=ss,
                              phi=phi, psi=psi, next_n=t["next_n"])
        aa.chain_id[:] = k[0]
        aa.res_id[:] = k[1]
        aa.ins_code[:] = k[2]
        aa.hetero[:] = bool(atom_array.hetero[t["idxs"][0]])
        implants.append(aa)
        info.append(dict(key=k, res_name=res_name, src=str(atom_array.res_name[t["idxs"][0]]),
                         phi=phi, psi=psi, ss=ss,
                         chi=tuple(float(c) for c in chis) if chis else (),
                         probability=float(proba), n_candidates=n_cand,
                         library=library))

    if not implants:
        raise ValueError("no residues were rebuilt; nothing to mutate")

    # ---- 削骨架 ----
    target_mask = np.zeros(len(atom_array), dtype=bool)
    for t in targets.values():
        target_mask[t["idxs"]] = True
    backbone_only = remove_side_chains(atom_array, target_mask)
    # remove_side_chains 会把目标残基改成 GLY, 这里改回目标名
    for k, t in targets.items():
        sel = ((backbone_only.chain_id == k[0])
               & (backbone_only.res_id == k[1])
               & (backbone_only.ins_code == k[2]))
        backbone_only.res_name[sel] = t["res_name"]

    # ---- 装配 (复用既有 graft 装配器; 传入的 mask 是残基级的) ----
    implant_array = bio_struct.concatenate(implants)
    mask_map = []
    for aa, k in zip(implants, [i["key"] for i in info]):
        b_mask = ((backbone_only.chain_id == k[0])
                  & (backbone_only.res_id == k[1])
                  & (backbone_only.ins_code == k[2]))
        i_mask = ((implant_array.chain_id == k[0])
                  & (implant_array.res_id == k[1])
                  & (implant_array.ins_code == k[2]))
        mask_map.extend([b_mask, i_mask])

    out = replace_side_chains(backbone_only, implant_array, mask_map)
    return (out, info) if return_info else out
