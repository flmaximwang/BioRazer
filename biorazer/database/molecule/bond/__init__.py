# -*- coding: utf-8 -*-
"""Covalent bond geometry (length / angle / dihedral), generic vs protein.

``biorazer.database.molecule.bond`` is a package split by geometry type and
by scope:

* :mod:`.length`  -- ideal covalent bond lengths (Å).
* :mod:`.angle`   -- ideal covalent bond angles (°).
* :mod:`.dihedral`-- ideal torsion angles (°) and torsion definitions.

Each geometry type further splits into ``generic`` (molecule-agnostic /
non-protein main tables, e.g. the Engh-Huber common bond values) and
``protein`` (protein-residue-specific tables: per-residue refinements,
side-chain geometry and the backbone/side-chain torsion data).

Every numeric entry carries the uniform record
``{mean, std, lb, up, source}`` (std/lb/up = ``np.nan`` when the source
provides no spread).
"""
