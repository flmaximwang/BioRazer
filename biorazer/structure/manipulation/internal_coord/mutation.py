# -*- coding: utf-8 -*-
"""Rebuild a residue's side chain inside an :class:`InternalCoord`.

操作对象是 :class:`~biorazer.structure.objects.internal_coords.InternalCoord` ——
侧链重建发生在**图上**, 不是在一份坐标上: 目标残基的骨架原子 (连同它在图里的
``phi``/``psi``/羰基 ``O`` 与 :attr:`InternalCoord.anchor` 条目) 原样留着, 该残基的
侧链原子**连同它们的整棵生成树**换掉。

与 :func:`~biorazer.structure.manipulation.atom_array.mutation.mutate` 的关系: 那个的
输入/输出是 biotite ``AtomArray`` (骨架取自输入、侧链按 rotamer 库重建后按掩码拼装);
这里做的是同一件事在 IC 图上的版本。两边的侧链生成树来自**同一张表**
(:data:`~biorazer.database.molecule.icoor.protein.topology.IC_PATH`) ——
模板构建与
:class:`~biorazer.structure.bridge.atom_array.AtomArray_InternalCoord`
的读入走的是它, 所以换进来的子树与换出去的子树同构, 只是原子数与 chi 值不同。

侧链几何来自
:func:`~biorazer.database.molecule.icoor.protein.template.build_template` 的单残基模板。
本模块**不**查 rotamer 库 (那是 :mod:`~biorazer.structure.manipulation.atom_array.mutation`
的事, 且需要 PyMOL / Rosetta 数据): ``rotamer`` 用模板自己的词表
(:func:`~biorazer.database.molecule.icoor.protein.template.rotamer_names`), 库外的
chi 用 ``chi`` 直接给 —— 两边都能接。

用法::

    from biorazer.structure.manipulation.internal_coord import mutate

    out = mutate(ic, ["A31H"])                      # 31 位 ALA -> HIS, 理想 canonical 侧链
    out = mutate(ic, ["A31H"], rotamer="g+/g+")     # 两个 chi 都用 g+ 箱均值 (60 度)
    out, info = mutate(ic, ["A31H"], chi={"CG": -65.0, "ND1": 170.0})
    out, info = mutate(ic, ["A31H"], return_info=True)
"""

from __future__ import annotations

import re

from biorazer.database.alphabet import AMINO_ACIDS_1TO3_UPPER
from biorazer.database.molecule.bond.dihedral.protein import SIDECHAIN_CHI
from biorazer.database.molecule.icoor.protein import template
from biorazer.database.molecule.icoor.protein.template import build_template
from biorazer.database.molecule.icoor.protein.topology import MAINCHAIN_ATOMS

from ...objects import InternalCoord, InternalCoordAtom

__all__ = [
    "mutate",
]

#: ``build_template`` 要一个二级结构类键, 但**侧链几何与它无关**: 实测 20 个标准残基
#: (``coil`` vs ``alpha-helix``, ``rotamer="canonical"``) 的 ``bond_distances`` /
#: ``bond_angles`` 与除羰基 ``O`` 的 ``(N, CA, C, O)`` 之外的每个 ``dihedra`` 条目
#: **逐键相同**; 而羰基 ``O`` 本来就不进结果 (骨架沿用输入)。所以这里取固定值,
#: 也不给调用方一个不起作用的 ``ss`` 参数。
_SS = "coil"

#: 骨架原子名 (留在输入那一份里, 不被侧链换掉)。
_BACKBONE = frozenset(MAINCHAIN_ATOMS)

_SPEC = re.compile(r"^([A-Za-z])(\d+)([A-Za-z])$")


def _residue_indices(ic):
    """``{(chain_id, res_id, ins_code): [原子下标, 按顺序]}``。"""
    groups: dict[tuple, list[int]] = {}
    for i, a in enumerate(ic.atoms):
        groups.setdefault((a.chain_id, a.res_id, a.ins_code), []).append(i)
    return groups


def _rotate_in(chi, resn, quads):
    """把调用方给的 ``chi`` 归一成 ``{rotator 原子名: 度}``。

    ``chi`` 可以是 ``{rotator: 度}``, 也可以是按 :data:`SIDECHAIN_CHI` 顺序给出的序列
    (与 :func:`~biorazer.structure.manipulation.atom_array.mutation.mutate` 返回的
    ``info["chi"]`` 同一顺序, 便于直接串起来用)。
    """
    if chi is None:
        return {}
    if isinstance(chi, dict):
        unknown = set(chi) - {q[3] for q in quads}
        if unknown:
            raise ValueError(
                f"{resn} has no chi ending at {sorted(unknown)}; "
                f"rotators are {[q[3] for q in quads]}")
        return {k: float(v) for k, v in chi.items()}
    values = [float(v) for v in chi]
    if len(values) != len(quads):
        raise ValueError(
            f"{resn} has {len(quads)} chi, got {len(values)} values "
            f"({[q[3] for q in quads]})")
    return {q[3]: v for q, v in zip(quads, values)}


def _graft(ic, key, tgt_resn, rotamer, chi):
    """把残基 ``key`` 的侧链换成 ``tgt_resn`` 的模板, 返回 ``(新 IC, info)``。"""
    idxs = _residue_indices(ic)[key]
    backbone, drop = {}, set()
    for i in idxs:
        name = ic.atoms[i].name
        if name in _BACKBONE:
            backbone.setdefault(name, i)
        else:
            drop.add(i)
    missing = [n for n in ("N", "CA", "C") if n not in backbone]
    if missing:
        raise ValueError(
            f"residue {key} has no {'/'.join(missing)} atom; "
            "a side chain cannot be grafted onto it")

    # build_template 对未知 rotamer 名**不报错** (rotamer_targets 把不认识的 bin
    # 直接跳过 -> 静默退回 canonical), 所以这一层自己校验, 免得调用方写错名字还拿到结果
    if rotamer not in template.rotamer_names(tgt_resn):
        raise ValueError(
            f"unknown rotamer {rotamer!r} for {tgt_resn}; "
            f"available: {template.rotamer_names(tgt_resn)}")
    tmpl, _ = build_template(tgt_resn, _SS, rotamer)
    tmpl_sc = [(n, a) for n, a in enumerate(tmpl.atoms) if a.name not in _BACKBONE]

    # ---- 新原子表: 丢掉目标残基侧链, 新侧链插在该残基骨架块之后 (编号随之右移) ----
    kept = [i for i in range(len(ic.atoms)) if i not in drop]
    cut = kept.index(max(backbone.values())) + 1    # 插入点: 目标残基骨架块的末尾
    head, tail = kept[:cut], kept[cut:]
    new_records = [InternalCoordAtom(ins_code=ic.atoms[backbone["N"]].ins_code,
                                     chain_id=key[0],
                                     res_name=tgt_resn,
                                     res_id=key[1],
                                     name=a.name,
                                     element=a.element)
                   for _, a in tmpl_sc]
    added = {rec.name: cut + n for n, rec in enumerate(new_records)}
    # 目标残基的**骨架**原子也要换成新记录: res_name/res_id 改在结果上, 且输入那份
    # 记录不被就地改写 (InternalCoordAtom 是可变的 —— 直接改就改了输入)
    backbone_new = {i: InternalCoordAtom(ins_code=ic.atoms[i].ins_code,
                                         chain_id=key[0],
                                         res_name=tgt_resn,
                                         res_id=key[1],
                                         name=ic.atoms[i].name,
                                         element=ic.atoms[i].element)
                    for i in backbone.values()}

    def pick(i):
        return backbone_new.get(i) or ic.atoms[i]

    # 插入点之后的原子整体右移 len(new_records) —— 不做这一步, 新侧链的编号会与
    # 被挤开的原子撞号 (撞号表现为 to_coords 报 "Inconsistent cycle")
    remap = {old: n for n, old in enumerate(head)}
    remap.update({old: cut + len(new_records) + n for n, old in enumerate(tail)})
    out = InternalCoord(atoms=[pick(i) for i in head]
                              + new_records
                              + [pick(i) for i in tail])

    # ---- 原图的几何原样搬, 只丢触及被删原子的条目 ----
    for (i, j), dist in ic.bond_distances.items():
        if i in remap and j in remap:
            out.bond_distances[(remap[i], remap[j])] = dist
    for (i, j, k), ang in ic.bond_angles.items():
        if i in remap and j in remap and k in remap:
            out.bond_angles[(remap[i], remap[j], remap[k])] = ang
    for quad, dih in ic.dihedra.items():
        if all(x in remap for x in quad):
            i, j, k, l = quad
            out.dihedra[(remap[i], remap[j], remap[k], remap[l])] = dih
    out.anchor = {remap[i]: c for i, c in ic.anchor.items() if i in remap}

    # ---- 接上模板的侧链子树: 模板里指向它自己骨架 (N/CA/C) 的父原子按**名字**
    #      改指输入的骨架, 所以 CB 的定位四元组与 chi1 都落在输入的真实骨架上 ----
    def here(name):
        if name in added:
            return added[name]
        if name in backbone:
            return remap[backbone[name]]
        raise ValueError(
            f"template for {tgt_resn} references atom {name!r}, which residue "
            f"{key} does not have")

    for quad, dih in tmpl.dihedra.items():
        names = tuple(tmpl.atoms[i].name for i in quad)
        if names[3] in _BACKBONE:
            # 模板的羰基 O / OXT 分支: 骨架 (含 O) 沿用输入, 不能被模板值覆盖
            continue
        out.dihedra[(here(names[0]), here(names[1]), here(names[2]), here(names[3]))] = dih
    for (k, l), dist in tmpl.bond_distances.items():
        if tmpl.atoms[k].name in _BACKBONE and tmpl.atoms[l].name in _BACKBONE:
            continue        # 模板 anchor 三元组自身的键: 输入的骨架里已有
        out.bond_distances[(here(tmpl.atoms[k].name), here(tmpl.atoms[l].name))] = dist
    for (j, k, l), ang in tmpl.bond_angles.items():
        if all(tmpl.atoms[i].name in _BACKBONE for i in (j, k, l)):
            continue
        out.bond_angles[(here(tmpl.atoms[j].name), here(tmpl.atoms[k].name),
                         here(tmpl.atoms[l].name))] = ang

    # ---- chi 覆盖 (可选): 模板已给了理想值, 这里按调用方给的改 ----
    chi_quads = SIDECHAIN_CHI.get(tgt_resn, [])
    for rotator, value in _rotate_in(chi, tgt_resn, chi_quads).items():
        for quad in chi_quads:
            if quad[3] == rotator:
                out.dihedra[(here(quad[0]), here(quad[1]), here(quad[2]), here(quad[3]))] = value

    info = dict(key=key,
                src=str(ic.atoms[idxs[0]].res_name),
                res_name=tgt_resn,
                rotamer=rotamer,
                chi=tuple(float(out.dihedra[(here(q[0]), here(q[1]), here(q[2]), here(q[3]))])
                          for q in chi_quads))
    return out, info


def mutate(ic, mutation_spec, *, chain=None, rotamer="canonical", chi=None,
           return_info=False):
    """把残基换成目标类型, 并在 IC 图上**重建**它的侧链。

    与 :func:`~biorazer.structure.manipulation.atom_array.mutation.mutate` 的区别只在
    操作对象与几何来源: 那个吃/吐 ``AtomArray``、侧链按 rotamer 库选;
    本函数吃/吐 :class:`InternalCoord`, 侧链来自模板的**理想**几何 + 调用方给的 chi。
    共同点: 骨架 (``N/CA/C/O``, 含羰基 ``O``) 与其余残基**一个原子都不动**。

    规格沿用既有约定 ``"<原残基字母><残基号><目标残基字母>"`` —— ``"A31H"`` 表示
    "31 位的 ALA 变成 HIS", 原残基字母会与图里的实际类型核对。

    Parameters
    ----------
    ic : InternalCoord
        输入结构 (不被改动; 返回的是新对象)。
    mutation_spec : list of str
        突变列表, 每项形如 ``"A31H"``。
    chain : str, optional
        只对该链做。``None`` = 所有链 (与既有函数一致: 按 ``res_id`` 匹配**全部链**,
        多链同号时会一起改; 要只改一条链就必须显式给 ``chain``)。
    rotamer : str
        模板词表里的 rotamer 名 (:func:`~biorazer.database.molecule.icoor.protein.
        template.rotamer_names`), 例如 ``"canonical"`` / ``"g+/g+"`` / ``"g-/t"``;
        写错的**名字**会报错 (可选值随残基不同: 1 个 chi 轴的是 ``g-/t/g+``,
        2 个轴的是 ``g-/g-`` 这种 ``<chi1>/<chi2>``)。它决定**理想** chi;
        与 ``chi`` 同时给时 ``chi`` 赢。
    chi : dict or sequence, optional
        ``{chi 末端原子名: 度}`` (如 ``{"CG": 60.0}``) 或按
        :data:`~biorazer.database.molecule.bond.dihedral.protein.SIDECHAIN_CHI`
        顺序给出的序列。给 rotamer 库选出来的 chi 走这里 —— 本模块不查库。
    return_info : bool
        是否连带返回每个残基实际落了什么 rotamer / chi。

    Returns
    -------
    InternalCoord, 或 ``(InternalCoord, list[dict])`` (``return_info=True``)。
    ``info`` 每项含 ``key``/``src``/``res_name``/``rotamer``/``chi`` —— ``chi`` 是
    从**结果图上回读**的值 (度), 不是入参。

    Raises
    ------
    ValueError
        规格格式错误、残基不存在、原残基字母对不上、目标残基缺 ``N/CA/C``、
        ``chi`` 的键/个数与目标残基不符。

    Examples
    --------
    >>> out = mutate(ic, ["A31H"])                  # doctest: +SKIP
    >>> out.atoms[3].res_name                       # doctest: +SKIP
    'HIS'
    >>> _, info = mutate(ic, ["A31H"], rotamer="g+", return_info=True)   # doctest: +SKIP
    >>> info[0]["chi"][0]                           # doctest: +SKIP
    60.0
    """
    if not isinstance(ic, InternalCoord):
        raise TypeError(f"ic must be an InternalCoord, got {type(ic).__name__}")
    if not isinstance(mutation_spec, (list, tuple)) or len(mutation_spec) == 0:
        raise ValueError("mutation_spec must be a non-empty list of strings")

    wanted: dict[tuple, tuple[str, str]] = {}
    for entry in mutation_spec:
        if not isinstance(entry, str):
            raise TypeError(
                f"mutation_spec entries must be str, got {type(entry).__name__}")
        m = _SPEC.fullmatch(entry.strip())
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

    out, info = ic, []
    for (want_chain, res_id), (src, tgt_resn) in wanted.items():
        groups = _residue_indices(out)
        keys = [k for k in groups
                if k[1] == res_id and (want_chain is None or k[0] == want_chain)]
        if not keys:
            raise ValueError(
                f"mutation_spec references residue {res_id}"
                + (f" of chain {want_chain}" if want_chain else "")
                + ", which does not exist in the input structure")
        expected = AMINO_ACIDS_1TO3_UPPER[src]
        for k in keys:
            cur = str(out.atoms[groups[k][0]].res_name)
            if cur != expected:
                raise ValueError(
                    f"residue {k} is {cur}, but the mutation spec expects "
                    f"{src} ({expected})")
            out, rec = _graft(out, k, tgt_resn, rotamer, chi)
            info.append(rec)

    return (out, info) if return_info else out
