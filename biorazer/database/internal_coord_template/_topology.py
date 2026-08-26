# -*- coding: utf-8 -*-
"""Per-amino-acid side-chain internal-coordinate (grow-path) topology.

This is the **canonical home** of the per-residue side-chain IC topology, shared
by the read path (``InternalCoord.from_atomarray`` in
``biorazer/structure/objects/internal_coords.py``) and the write path
(``biorazer.database.internal_coord_template``).  (It was originally authored as
``biorazer/database/amino_acid_internal_coords.py``.)

The main chain (backbone) of a protein is a **uniform** polymer walk
``N -> CA -> C -> O`` with a peptide bond ``C - N`` to the next residue; every
residue follows the same path, so it is built algorithmically (see
``InternalCoord.from_atomarray``) and is *not* stored here.

The side chain, by contrast, is **different for each amino acid**: a tree
grafted at ``CA``, whose branching is described by the official chi rotamers
(``chi1`` rotates about ``CA-CB``, ``chi2`` about ``CB-CG``, ...).  This module
stores, per residue, the ordered list of grow specs that traverse that tree.

Each grow spec is a quad ``(i, j, k, l)`` of **atom names**, matching the
``InternalCoord`` convention: atom ``l`` is grown from parents ``(i, j, k)``,
where ``k`` is the *bonded* parent of ``l`` (``bond (k, l)``), the bond angle is
at ``k`` (``bond_angles[(j, k, l)]``), and the dihedral is
``dihedra[(i, j, k, l)]``.  The quads are listed in topological order: every
``i/j/k`` is either a backbone atom (``N``/``CA``/``C``/``O``, already placed by
the main-chain pass) or an earlier side-chain atom.

The table is translated 1:1 from Rosetta's ``fa_standard`` residue params
(``ICOOR_INTERNAL`` atom-tree).  Source files, release 408:

    main/database/chemical/residue_type_sets/fa_standard/residue_types/l-caa/<AA>.params

Each ``ICOOR_INTERNAL <child> <dih> <ang> <len> <P1> <P2> <P3>`` row places the
child from the stub frame ``(P3, P2, P1)``: ``P1`` is the *bonded* parent
(bond ``(P1, child)``), the bond angle is at ``P1`` between ``P2`` and ``child``,
and the dihedral is ``dihedral(P3, P2, P1, child)`` (Rosetta measures it around
the ``P2-P1`` axis in the ``(P3, P2, P1)`` frame).  The IC quad therefore keeps
**the official order** ``(P3, P2, P1, child)``: the bonded parent ``P1`` stays
in the third slot ``k``, and the stored dihedral equals the official torsion
definition (e.g. ``chi1 = (N, CA, CB, CG)``).

Backbone ``O``/``OXT`` are placed by the main-chain pass (Biopython convention
``(N, CA, C, O)``); hydrogens and Pro's virtual ring-closure atom ``NV`` are not
represented (heavy-atom structure only).  ``CH2`` (Trp indole carbon) is a
heavy atom despite the ``H`` in its name.
"""

from __future__ import annotations

#: Backbone atoms placed by the uniform main-chain pass (not in the side-chain
#: table).  ``O``/``OXT`` are branches off ``C``; ``H_n`` backbone protons exist
#: only in explicit-H structures.
MAINCHAIN_ATOMS = frozenset(("N", "CA", "C", "O", "OXT"))

#: Side-chain grow-path per residue: ``{res_name: ((i, j, k, l), ...)}``, one
#: ``(i, j, k, l)`` atom-name quad per side-chain heavy atom ``l`` (atom ``k``
#: is its bonded parent).  Quads are in the **official dihedral order**
#: ``(P3, P2, P1, child)`` from Rosetta ICOOR, so the first side-chain quads
#: equal the official chi definitions (``chi1 = (N, CA, CB, CG)``,
#: ``chi2 = (CA, CB, CG, CD)``, ...).  ``GLY`` has no side chain so its tuple
#: is empty.
SIDE_CHAIN_IC_PATH = {
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
