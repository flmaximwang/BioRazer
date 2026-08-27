# -*- coding: utf-8 -*-
"""Ideal covalent bond lengths (Å).

* :mod:`.protein` -- a package whose single module :mod:`.protein.by_residue`
  holds all protein bond-length reference data: the generic main table
  (``AMINO_ACID_BOND_LENGTH``, Engh-Huber common values keyed by atom-name
  pair, each key referenced by ≥2 residues), the per-residue refinements
  (``AMINO_ACID_BOND_LENGTH_BY_RESIDUE``), and the per-residue side-chain
  bond lengths (``AMINO_ACID_SIDECHAIN_BOND``, Rosetta ICOOR ideals).

Every entry is ``{mean, std, lb, up, source}`` (``np.nan`` when unknown).
"""
