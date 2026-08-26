# -*- coding: utf-8 -*-
"""Per-molecule physico-chemical reference data, split by property.

``biorazer.database.molecule`` groups all molecule-level geometry / property
tables under one package so future molecule classes (nucleic acids, ligands,
...) have room to grow alongside the protein data:

* :mod:`.atom`    -- per-atom properties (van der Waals radius, charge).
* :mod:`.bond`    -- covalent bond geometry: length / angle / dihedral.
* :mod:`.icoor`   -- internal-coordinate reference data (grow-path topology
  and ideal per-residue templates for proteins).

Within :mod:`.bond`, each geometry type splits into ``generic``
(molecule-agnostic main tables) and ``protein`` (protein-residue-specific
refinements).  Every numeric entry carries the uniform record
``{mean, std, lb, up, source}`` (std/lb/up = ``np.nan`` when the source
provides no spread).
"""
