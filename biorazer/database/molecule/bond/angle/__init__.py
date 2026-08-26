# -*- coding: utf-8 -*-
"""Ideal covalent bond angles (°).

* :mod:`.generic` -- molecule-agnostic main table (Engh-Huber common
  backbone angle values, keyed by the atom-name triple, vertex in the
  middle).
* :mod:`.protein` -- protein-specific refinements: per-residue backbone
  values (``*_BY_RESIDUE``) and the per-residue side-chain bond angles
  (Rosetta ICOOR ideals).

Every entry is ``{mean, std, lb, up, source}`` (``np.nan`` when unknown).
"""
