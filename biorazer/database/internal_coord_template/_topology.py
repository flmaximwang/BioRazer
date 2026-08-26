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
(``ICOOR_INTERNAL`` atom-tree) -- the first parent there is the bonded parent
``k``; we swap it to the third slot.  Source files, release 408:

    main/database/chemical/residue_type_sets/fa_standard/residue_types/l-caa/<AA>.params

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
#: is its bonded parent).  ``GLY`` has no side chain so its tuple is empty.
SIDE_CHAIN_IC_PATH = {
    "ALA": (
        ("N", "C", "CA", "CB"),
    ),
    "ARG": (
        ("N", "C", "CA", "CB"),
        ("CA", "N", "CB", "CG"),
        ("CB", "CA", "CG", "CD"),
        ("CG", "CB", "CD", "NE"),
        ("CD", "CG", "NE", "CZ"),
        ("NE", "CD", "CZ", "NH1"),
        ("NE", "NH1", "CZ", "NH2"),
    ),
    "ASN": (
        ("N", "C", "CA", "CB"),
        ("CA", "N", "CB", "CG"),
        ("CB", "CA", "CG", "OD1"),
        ("CB", "OD1", "CG", "ND2"),
    ),
    "ASP": (
        ("N", "C", "CA", "CB"),
        ("CA", "N", "CB", "CG"),
        ("CB", "CA", "CG", "OD1"),
        ("CB", "OD1", "CG", "OD2"),
    ),
    "CYS": (
        ("N", "C", "CA", "CB"),
        ("CA", "N", "CB", "SG"),
    ),
    "GLN": (
        ("N", "C", "CA", "CB"),
        ("CA", "N", "CB", "CG"),
        ("CB", "CA", "CG", "CD"),
        ("CG", "CB", "CD", "OE1"),
        ("CG", "OE1", "CD", "NE2"),
    ),
    "GLU": (
        ("N", "C", "CA", "CB"),
        ("CA", "N", "CB", "CG"),
        ("CB", "CA", "CG", "CD"),
        ("CG", "CB", "CD", "OE1"),
        ("CG", "OE1", "CD", "OE2"),
    ),
    "GLY": (),
    "HIS": (
        ("N", "C", "CA", "CB"),
        ("CA", "N", "CB", "CG"),
        ("CB", "CA", "CG", "ND1"),
        ("CG", "CB", "ND1", "CE1"),
        ("ND1", "CG", "CE1", "NE2"),
        ("CE1", "ND1", "NE2", "CD2"),
    ),
    "ILE": (
        ("N", "C", "CA", "CB"),
        ("CA", "N", "CB", "CG1"),
        ("CB", "CA", "CG1", "CD1"),
        ("CA", "CG1", "CB", "CG2"),
    ),
    "LEU": (
        ("N", "C", "CA", "CB"),
        ("CA", "N", "CB", "CG"),
        ("CB", "CA", "CG", "CD1"),
        ("CB", "CD1", "CG", "CD2"),
    ),
    "LYS": (
        ("N", "C", "CA", "CB"),
        ("CA", "N", "CB", "CG"),
        ("CB", "CA", "CG", "CD"),
        ("CG", "CB", "CD", "CE"),
        ("CD", "CG", "CE", "NZ"),
    ),
    "MET": (
        ("N", "C", "CA", "CB"),
        ("CA", "N", "CB", "CG"),
        ("CB", "CA", "CG", "SD"),
        ("CG", "CB", "SD", "CE"),
    ),
    "PHE": (
        ("N", "C", "CA", "CB"),
        ("CA", "N", "CB", "CG"),
        ("CB", "CA", "CG", "CD1"),
        ("CG", "CB", "CD1", "CE1"),
        ("CD1", "CG", "CE1", "CZ"),
        ("CE1", "CD1", "CZ", "CE2"),
        ("CZ", "CE1", "CE2", "CD2"),
    ),
    "PRO": (
        ("N", "C", "CA", "CB"),
        ("CA", "N", "CB", "CG"),
        ("CB", "CA", "CG", "CD"),
    ),
    "SER": (
        ("N", "C", "CA", "CB"),
        ("CA", "N", "CB", "OG"),
    ),
    "THR": (
        ("N", "C", "CA", "CB"),
        ("CA", "N", "CB", "OG1"),
        ("CA", "OG1", "CB", "CG2"),
    ),
    "TRP": (
        ("N", "C", "CA", "CB"),
        ("CA", "N", "CB", "CG"),
        ("CB", "CA", "CG", "CD1"),
        ("CG", "CB", "CD1", "NE1"),
        ("CD1", "CG", "NE1", "CE2"),
        ("NE1", "CD1", "CE2", "CZ2"),
        ("CE2", "NE1", "CZ2", "CH2"),
        ("CZ2", "CE2", "CH2", "CZ3"),
        ("CH2", "CZ2", "CZ3", "CE3"),
        ("CZ3", "CH2", "CE3", "CD2"),
    ),
    "TYR": (
        ("N", "C", "CA", "CB"),
        ("CA", "N", "CB", "CG"),
        ("CB", "CA", "CG", "CD1"),
        ("CG", "CB", "CD1", "CE1"),
        ("CB", "CD1", "CG", "CD2"),
        ("CG", "CB", "CD2", "CE2"),
        ("CD2", "CG", "CE2", "CZ"),
        ("CE2", "CD2", "CZ", "OH"),
    ),
    "VAL": (
        ("N", "C", "CA", "CB"),
        ("CA", "N", "CB", "CG1"),
        ("CA", "CG1", "CB", "CG2"),
    ),
}