# -*- coding: utf-8 -*-
"""Protein torsion (dihedral) data, split by classification axis.

``biorazer.database.molecule.bond.dihedral.protein`` is a package with two
modules:

* :mod:`.by_ss`      -- backbone torsion data classified by secondary
  structure: phi/psi/omega per DSSP 9-state class
  (:data:`.by_ss.SS_BB_TORSION_ANGLE`), the DSSP code map, beta-turn
  dihedrals (PROMOTIF), the cis/trans omega constants and the official
  main-chain torsion definitions.
* :mod:`.by_residue` -- per-residue (side-chain) torsion data: the official
  chi definitions (:data:`.by_residue.SIDECHAIN_CHI`), the canonical
  IC-frame side-chain dihedrals (Rosetta ICOOR ideals), and the Dunbrack
  rotamer bin framework.

Backward-compatible re-exports of the former single ``protein`` module's
names are provided here; every public name also resolves from its per-module
home.

所有角度单位 **度 (degree)**。数值记录统一为 ``{mean, std, lb, up,
source}``; 查不到 spread 的字段为 ``np.nan``。
"""

from .by_ss import (  # noqa: F401
    SS_BB_TORSION_ANGLE,
    DSSP_SS_CODE,
    BB_TORSION_TURNS,
    BB_TORSION_REFS,
    OMEGA_TRANS,
    OMEGA_CIS,
    MAINCHAIN_TORSION_DEFINITIONS,
)
from .by_residue import (  # noqa: F401
    AAS,
    SIDECHAIN_CHI,
    SIDECHAIN_IC_DIHEDRAL,
    ROTAMER_BIN,
    NON_ROTAMERIC_BIN_WIDTH,
    DUNBRACK_ROTAMERS,
    SIDECHAIN_ROTAMER_LIB,
    SIDECHAIN_NON_ROTAMERIC_BINS,
    SIDECHAIN_DIHE_REFS,
)

__all__ = [
    # backbone (by secondary structure)
    "SS_BB_TORSION_ANGLE", "DSSP_SS_CODE", "BB_TORSION_TURNS",
    "BB_TORSION_REFS", "OMEGA_TRANS", "OMEGA_CIS",
    "MAINCHAIN_TORSION_DEFINITIONS",
    # side-chain (by residue)
    "AAS", "SIDECHAIN_CHI", "SIDECHAIN_IC_DIHEDRAL", "ROTAMER_BIN",
    "NON_ROTAMERIC_BIN_WIDTH", "DUNBRACK_ROTAMERS",
    "SIDECHAIN_ROTAMER_LIB", "SIDECHAIN_NON_ROTAMERIC_BINS",
    "SIDECHAIN_DIHE_REFS",
]
