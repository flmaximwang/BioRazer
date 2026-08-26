# -*- coding: utf-8 -*-
"""Backbone / side-chain torsion angles for amino acids.

``biorazer.database.torsion_angle`` is now a package with two modules:

* :mod:`.backbone`  -- secondary-structure-classified phi/psi/omega, beta-turn
  dihedrals (PROMOTIF), SS codes and the literature/ref registry.
* :mod:`.sidechain` -- the Dunbrack rotamer library: per-residue chi_torsions
  (definitions + common rotamer mean angles).

Backward-compatible re-exports of the old ``biorazer.database.torsion_angle``
module's names are provided here (from :mod:`.backbone`).
"""

from .backbone import (  # noqa: F401
    SS_BB_TORSION_ANGLE,
    DSSP_SS_CODE,
    BB_TORSION_TURNS,
    BB_TORSION_REFS,
    OMEGA_TRANS,
    OMEGA_CIS,
)