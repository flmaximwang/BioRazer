# -*- coding: utf-8 -*-
"""Ideal covalent bond angles (°).

* :mod:`.generic` -- molecule-agnostic **backbone** main table
  (Engh-Huber common angle values, keyed by the atom-name triple, vertex in
  the middle; flat, no residue key).
* :mod:`.protein` -- protein-specific, **residue-keyed** bond angles: the
  merged ``AMINO_ACID_BOND_ANGLE`` holds per-residue backbone refinements
  (Gly/Pro/Ala/VIT, Engh-Huber) plus the 20 AAs' side-chain bond angles
  (Rosetta ICOOR ideals).

Every entry is ``{mean, std, lb, up, source}`` (``np.nan`` when unknown).
"""
