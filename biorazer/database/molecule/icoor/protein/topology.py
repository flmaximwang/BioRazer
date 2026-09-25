# -*- coding: utf-8 -*-
"""Protein internal-coordinate (grow-path) topology.

This is the **canonical home** of the per-residue side-chain IC topology, shared
by the read path
(:class:`~biorazer.structure.bridge.atom_array.AtomArray_InternalCoord`) and the
write path (``biorazer.database.molecule.icoor.protein.template``).  (It was
originally authored as ``biorazer/database/amino_acid_internal_coords.py``, then
``internal_coord_template/_topology.py``.)

The main chain (backbone) of a protein is a **uniform** polymer walk
``N -> CA -> C -> O`` with a peptide bond ``C - N`` to the next residue;
every residue follows the same path, so unlike the side chain it is not a
per-residue table.  Its grow specs are recorded once in
:data:`BACKBONE_IC_PATH` (grouped into per-residue ``intra`` carbonyl
branches and cross-residue ``peptide`` links), which both the read path
(:class:`~biorazer.structure.bridge.atom_array.AtomArray_InternalCoord`) and
the write path (the template builder) consume.

The side chain, by contrast, is **different for each amino acid**: a tree
grafted at ``CA``, whose branching is described by the official chi
rotamers (``chi1`` rotates about ``CA-CB``, ``chi2`` about ``CB-CG``, ...).
This module stores, per residue, the ordered list of grow specs that
traverse that tree (:data:`IC_PATH`), plus the **uniform backbone** grow
specs every residue shares (:data:`BACKBONE_IC_PATH`).

Each grow spec is a quad ``(i, j, k, l)`` of **atom names**, matching the
``InternalCoord`` convention: atom ``l`` is grown from parents ``(i, j, k)``,
where ``k`` is the *bonded* parent of ``l`` (``bond (k, l)``), the bond
angle is at ``k`` (``bond_angles[(j, k, l)]``), and the dihedral is
``dihedra[(i, j, k, l)]``.  The quads are listed in topological order:
every ``i/j/k`` is either a backbone atom (``N``/``CA``/``C``/``O``,
already placed by the main-chain pass) or an earlier side-chain atom.

The table is translated 1:1 from Rosetta's ``fa_standard`` residue params
(``ICOOR_INTERNAL`` atom-tree).  Source files, release 408:

    main/database/chemical/residue_type_sets/fa_standard/residue_types/l-caa/<AA>.params

Each ``ICOOR_INTERNAL <child> <dih> <ang> <len> <P1> <P2> <P3>`` row places
the child from the stub frame ``(P3, P2, P1)``: ``P1`` is the *bonded*
parent (bond ``(P1, child)``), the bond angle is at ``P1`` between ``P2``
and ``child``, and the dihedral is ``dihedral(P3, P2, P1, child)`` (Rosetta
measures it around the ``P2-P1`` axis in the ``(P3, P2, P1)`` frame).  The
IC quad therefore keeps **the official order** ``(P3, P2, P1, child)``:
the bonded parent ``P1`` stays in the third slot ``k``, and the stored
dihedral equals the official torsion definition (e.g. ``chi1 =
(N, CA, CB, CG)``).

Backbone ``O``/``OXT`` are placed by the main-chain pass (Biopython
convention ``(N, CA, C, O)``); hydrogens and Pro's virtual ring-closure
atom ``NV`` are not represented (heavy-atom structure only).  ``CH2`` (Trp
indole carbon) is a heavy atom despite the ``H`` in its name.
"""

from __future__ import annotations

#: Backbone atoms placed by the uniform main-chain pass (not in the
#: side-chain table).  ``O``/``OXT`` are branches off ``C``; ``H_n``
#: backbone protons exist only in explicit-H structures.
MAINCHAIN_ATOMS = frozenset(("N", "CA", "C", "O", "OXT"))

#: Backbone (main-chain) grow-path: the **uniform** IC grow quads that
#: every residue shares (the backbone is not per-residue, unlike
#: :data:`IC_PATH`).  Two groups:
#:
#: * ``"intra"`` -- per-residue quads: the carbonyl ``O`` (and C-terminal
#:   ``OXT``) branch off ``C``, grown from the ``(N, CA, C)`` anchor frame.
#:   This is the only backbone growth a single-residue builder can record
#:   (``N``/``CA``/``C`` themselves are the anchor frame, placed
#:   analytically, not grown).  ``OXT`` is present only in C-terminal
#:   residues and in explicit-H / capped structures; the bond tables carry
#:   no ``OXT`` geometry, so template builders record just the ``O`` quad.
#:   The ``O`` quad's dihedral is **not** a fixed 180: a trans peptide plane
#:   places O anti to the next residue's amide N across the C-N bond, so
#:   ``dihedral(N, CA, C, O) = psi - 180`` with ``psi`` the carbonyl
#:   residue's own psi (``dihedral(N, CA, C, N_{i+1})``).  For a terminal
#:   residue a trans plane is assumed (``psi = 180`` -> ``O`` dihedral 0).
#: * ``"peptide"`` -- cross-residue quads linking residue ``i`` to ``i+1``,
#:   with atom names subscripted ``_i`` (current residue) / ``_{i+1}``
#:   (next residue).  These are the main-chain pass of
#:   :class:`~biorazer.structure.bridge.atom_array.AtomArray_InternalCoord`:
#:   each quad grows one atom of residue
#:   ``i+1`` from the ``(i, j, k)`` frame spanning the peptide bond.  The
#:   quads are listed in the **official torsion order** of
#:   :data:`~biorazer.database.molecule.bond.dihedral.protein.ALIAS_QUAD`
#:   -- the dihedral each quad stores is the official definition:
#:   ``(N_i, CA_i, C_i, N_{i+1})`` = ``psi`` of residue ``i``,
#:   ``(CA_i, C_i, N_{i+1}, CA_{i+1})`` = ``omega`` of residue ``i``,
#:   ``(C_i, N_{i+1}, CA_{i+1}, C_{i+1})`` = ``phi`` of residue ``i+1``.
BACKBONE_IC_PATH = {
    "intra": (
        ("N", "CA", "C", "O"),
        ("N", "CA", "C", "OXT"),
    ),
    "peptide": (
        ("N_i", "CA_i", "C_i", "N_{i+1}"),
        ("CA_i", "C_i", "N_{i+1}", "CA_{i+1}"),
        ("C_i", "N_{i+1}", "CA_{i+1}", "C_{i+1}"),
    ),
}


def carbonyl_o_dihedral(psi: float) -> float:
    """羰基 ``O`` 的放置二面角 ``(N, CA, C, O)`` —— 由羰基碳 C 的 sp2 共面性给出。

    **单一定义处**: 模板构建 (:func:`~biorazer.database.molecule.icoor.protein.
    template.build_template`)、真实结构读入 (:class:`~biorazer.structure.bridge.
    atom_array.AtomArray_InternalCoord`) 与真实骨架上的重建
    (:func:`~biorazer.structure.manipulation.atom_array.mutation.build_side_chain`) 都调用
    这里, 不要各自再写一份 ``- 180``。

    推导 (为什么 ``- 180`` 与键角数值无关)
    ─────────────────────────────────────
    羰基碳 C 是 sp2 中心, 三个取代基 ``CA``/``O``/``N_{i+1}`` 共面; 而绕轴
    ``CA-C`` 的这根轴**本身就在该平面内**, 所以 ``C->O`` 与 ``C->N_{i+1}`` 在
    垂直于轴的平面里的投影必然共线; 又由 ``angle(O-C-N) = 360 - angle(CA-C-O)
    - angle(CA-C-N)`` 知 O 与 ``N_{i+1}`` 分居 ``C->CA`` 射线两侧, 两个投影恰好
    **反平行**。于是

        dihedral(N, CA, C, O) = dihedral(N, CA, C, N_{i+1}) - 180 = psi - 180

    只用到共面性, **不含任何键角数值** (键角取 130/100 或 150/80 结果一样)。
    反过来: 这个二面角**不是**自由度, 它是 ``psi`` 的确定性函数 —— 真正自由的
    是 ``psi`` (绕 ``CA-C`` 的旋转)。所以只要知道 ``N_{i+1}``, ``O`` 就被唯一
    确定; ``psi`` 或所属二级结构类的均值只是"没有 ``N_{i+1}`` 时"的退路。

    实测残差 (6VY1 链 A, n=120, 用仓库自己的 dihedral):
    ``dihedral(N,CA,C,O) - (psi - 180)`` 均值 ``-0.50`` 度 / 最大 ``5.79`` 度 ——
    这是真实结构的**出平面量**, 不是自由度: O 到 ``(CA, C, N_{i+1})`` 平面距离
    均值 ``0.020`` A / 最大 ``0.108`` A, 而 ``asin(0.108 / 1.057) = 5.9`` 度,
    两者自洽 (1.057 A = O 绕轴的半径 ``1.231 * sin(120.8)``)。C 上三键角之和
    实为 ``360.0 - 0.02`` 度 (max 偏差 0.27 度), sp2 平面性成立。

    Parameters
    ----------
    psi : float
        该残基的 ``psi = (N_i, CA_i, C_i, N_{i+1})`` (**度**)。没有真实
        ``N_{i+1}`` 时可传该残基所属二级结构类的均值 ``psi`` (近似)。

    Returns
    -------
    float
        ``(N, CA, C, O)`` 二面角 (**度**, 落在 ``[0, 360)``)。

    Examples
    --------
    >>> carbonyl_o_dihedral(-45.0)          # alpha 螺旋 -> O 在 135 度
    135.0
    >>> carbonyl_o_dihedral(130.0)          # beta 折叠 -> O 在 -50 度
    310.0
    """
    return (float(psi) - 180.0) % 360.0


#: Side-chain grow-path per residue: ``{res_name: ((i, j, k, l), ...)}``,
#: one ``(i, j, k, l)`` atom-name quad per side-chain heavy atom ``l``
#: (atom ``k`` is its bonded parent).  Quads are in the **official dihedral
#: order** ``(P3, P2, P1, child)`` from Rosetta ICOOR, so the first
#: side-chain quads equal the official chi definitions (``chi1 =
#: (N, CA, CB, CG)``, ``chi2 = (CA, CB, CG, CD)``, ...).  ``GLY`` has no
#: side chain so its tuple is empty.
IC_PATH = {
    "ALA": (
        ("C", "N", "CA", "CB"),
    ),
    "ARG": (
        ("C", "N", "CA", "CB"),
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "CD"),
        ("CB", "CG", "CD", "NE"),
        ("CG", "CD", "NE", "CZ"),
        ("CD", "NE", "CZ", "NH1"),
        ("NH1", "NE", "CZ", "NH2"),
    ),
    "ASN": (
        ("C", "N", "CA", "CB"),
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "OD1"),
        ("OD1", "CB", "CG", "ND2"),
    ),
    "ASP": (
        ("C", "N", "CA", "CB"),
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "OD1"),
        ("OD1", "CB", "CG", "OD2"),
    ),
    "CYS": (
        ("C", "N", "CA", "CB"),
        ("N", "CA", "CB", "SG"),
    ),
    "GLN": (
        ("C", "N", "CA", "CB"),
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "CD"),
        ("CB", "CG", "CD", "OE1"),
        ("OE1", "CG", "CD", "NE2"),
    ),
    "GLU": (
        ("C", "N", "CA", "CB"),
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "CD"),
        ("CB", "CG", "CD", "OE1"),
        ("OE1", "CG", "CD", "OE2"),
    ),
    "GLY": (),
    "HIS": (
        ("C", "N", "CA", "CB"),
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "ND1"),
        ("CB", "CG", "ND1", "CE1"),
        ("CG", "ND1", "CE1", "NE2"),
        ("ND1", "CE1", "NE2", "CD2"),
    ),
    "ILE": (
        ("C", "N", "CA", "CB"),
        ("N", "CA", "CB", "CG1"),
        ("CA", "CB", "CG1", "CD1"),
        ("CG1", "CA", "CB", "CG2"),
    ),
    "LEU": (
        ("C", "N", "CA", "CB"),
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "CD1"),
        ("CD1", "CB", "CG", "CD2"),
    ),
    "LYS": (
        ("C", "N", "CA", "CB"),
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "CD"),
        ("CB", "CG", "CD", "CE"),
        ("CG", "CD", "CE", "NZ"),
    ),
    "MET": (
        ("C", "N", "CA", "CB"),
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "SD"),
        ("CB", "CG", "SD", "CE"),
    ),
    "PHE": (
        ("C", "N", "CA", "CB"),
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "CD1"),
        ("CB", "CG", "CD1", "CE1"),
        ("CG", "CD1", "CE1", "CZ"),
        ("CD1", "CE1", "CZ", "CE2"),
        ("CE1", "CZ", "CE2", "CD2"),
    ),
    "PRO": (
        ("C", "N", "CA", "CB"),
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "CD"),
    ),
    "SER": (
        ("C", "N", "CA", "CB"),
        ("N", "CA", "CB", "OG"),
    ),
    "THR": (
        ("C", "N", "CA", "CB"),
        ("N", "CA", "CB", "OG1"),
        ("OG1", "CA", "CB", "CG2"),
    ),
    "TRP": (
        ("C", "N", "CA", "CB"),
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "CD1"),
        ("CB", "CG", "CD1", "NE1"),
        ("CG", "CD1", "NE1", "CE2"),
        ("CD1", "NE1", "CE2", "CZ2"),
        ("NE1", "CE2", "CZ2", "CH2"),
        ("CE2", "CZ2", "CH2", "CZ3"),
        ("CZ2", "CH2", "CZ3", "CE3"),
        ("CH2", "CZ3", "CE3", "CD2"),
    ),
    "TYR": (
        ("C", "N", "CA", "CB"),
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "CD1"),
        ("CB", "CG", "CD1", "CE1"),
        ("CD1", "CB", "CG", "CD2"),
        ("CB", "CG", "CD2", "CE2"),
        ("CG", "CD2", "CE2", "CZ"),
        ("CD2", "CE2", "CZ", "OH"),
    ),
    "VAL": (
        ("C", "N", "CA", "CB"),
        ("N", "CA", "CB", "CG1"),
        ("CG1", "CA", "CB", "CG2"),
    ),
}
