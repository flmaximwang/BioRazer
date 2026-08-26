# -*- coding: utf-8 -*-
"""Ideal covalent bond lengths (Å).

* :mod:`.generic` -- molecule-agnostic main table (Engh-Huber common
  backbone bond values, keyed by the atom-name pair).
* :mod:`.protein` -- protein-specific refinements: per-residue backbone
  values (``*_BY_RESIDUE``) and the per-residue side-chain bond lengths
  (Rosetta ICOOR ideals).

Every entry is ``{mean, std, lb, up, source}`` (``np.nan`` when unknown).
"""
