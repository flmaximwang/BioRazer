# -*- coding: utf-8 -*-
"""Protein bond-length reference data (merged into a single ``by_residue`` module).

``biorazer.database.molecule.bond.length.protein`` is a package with a single
module :mod:`.by_residue`, which holds *all* protein bond-length tables:

* ``AMINO_ACID_BOND_LENGTH`` -- generic main table (Engh-Huber common values,
  keyed by atom-name pair; every key referenced by ≥2 residues).
* ``AMINO_ACID_BOND_LENGTH_BY_RESIDUE`` -- per-residue refinements
  (Gly/Pro/Ala / VIT group) that override the generic values.
* ``AMINO_ACID_SIDECHAIN_BOND`` + ``BOND_SIDECHAIN_REFS`` -- per-residue
  side-chain ideal bond lengths (Rosetta 408 ICOOR ideals).

The former flat modules ``length/generic.py`` and ``length/protein.py`` were
consolidated here; backward-compatible re-exports let every consumer keep
working unchanged.  All values are ``{mean, std, lb, up, source}`` records
(``np.nan`` when no spread is known).
"""

from .by_residue import (  # noqa: F401
    BOND_REFS,
    AMINO_ACID_BOND_LENGTH,
    AMINO_ACID_BOND_LENGTH_BY_RESIDUE,
    AMINO_ACID_SIDECHAIN_BOND,
    BOND_SIDECHAIN_REFS,
)

__all__ = [
    "BOND_REFS",
    "AMINO_ACID_BOND_LENGTH",
    "AMINO_ACID_BOND_LENGTH_BY_RESIDUE",
    "AMINO_ACID_SIDECHAIN_BOND",
    "BOND_SIDECHAIN_REFS",
]
