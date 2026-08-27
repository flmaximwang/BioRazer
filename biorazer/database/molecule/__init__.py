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

Convenience re-exports (the full public names also resolve from their
per-module homes)::

    from biorazer.database.molecule import (
        AMINO_ACID_BOND_LENGTH,            # bond.length.protein (by_residue)
        AMINO_ACID_BOND_ANGLE,             # bond.angle.generic
        SS_BB_TORSION_ANGLE,               # bond.dihedral.protein.by_ss
        IC_PATH,                           # icoor.protein.topology
        BACKBONE_IC_PATH,                  # icoor.protein.topology
        ATOM_RADIUS,                       # atom.radius
    )
"""

# bond geometry
from .bond.length.protein import (  # noqa: F401
    BOND_REFS,
    AMINO_ACID_BOND_LENGTH,
    AMINO_ACID_BOND_LENGTH_BY_RESIDUE,
    AMINO_ACID_SIDECHAIN_BOND,
    BOND_SIDECHAIN_REFS,
)
from .bond.angle.generic import AMINO_ACID_BOND_ANGLE  # noqa: F401
from .bond.angle.protein import (  # noqa: F401
    AMINO_ACID_BOND_ANGLE_BY_RESIDUE,
    AMINO_ACID_SIDECHAIN_BOND_ANGLE,
)
from .bond.dihedral.protein import (  # noqa: F401
    SS_BB_TORSION_ANGLE,
    DSSP_SS_CODE,
    BB_TORSION_TURNS,
    BB_TORSION_REFS,
    OMEGA_TRANS,
    OMEGA_CIS,
    MAINCHAIN_TORSION_DEFINITIONS,
    SIDECHAIN_CHI,
    ROTAMER_BIN,
    NON_ROTAMERIC_BIN_WIDTH,
    SIDECHAIN_ROTAMER_LIB,
    SIDECHAIN_NON_ROTAMERIC_BINS,
    SIDECHAIN_DIHE_REFS,
    AAS,
)

# atom properties
from .atom.radius import ATOM_RADIUS, vdw_dict, vdw_radii  # noqa: F401

# internal-coordinate reference data
from .icoor.protein.topology import (  # noqa: F401
    IC_PATH,
    MAINCHAIN_ATOMS,
    BACKBONE_IC_PATH,
)

__all__ = [
    # bond length
    "AMINO_ACID_BOND_LENGTH", "AMINO_ACID_BOND_LENGTH_BY_RESIDUE",
    "AMINO_ACID_SIDECHAIN_BOND", "BOND_REFS", "BOND_SIDECHAIN_REFS",
    # bond angle
    "AMINO_ACID_BOND_ANGLE", "AMINO_ACID_BOND_ANGLE_BY_RESIDUE",
    "AMINO_ACID_SIDECHAIN_BOND_ANGLE",
    # bond dihedral
    "SS_BB_TORSION_ANGLE", "DSSP_SS_CODE", "BB_TORSION_TURNS",
    "BB_TORSION_REFS", "OMEGA_TRANS", "OMEGA_CIS",
    "MAINCHAIN_TORSION_DEFINITIONS",
    "SIDECHAIN_CHI", "ROTAMER_BIN",
    "NON_ROTAMERIC_BIN_WIDTH",
    "SIDECHAIN_ROTAMER_LIB", "SIDECHAIN_NON_ROTAMERIC_BINS",
    "SIDECHAIN_DIHE_REFS", "AAS",
    # atom properties
    "ATOM_RADIUS", "vdw_dict", "vdw_radii",
    # icoor topology
    "IC_PATH", "MAINCHAIN_ATOMS", "BACKBONE_IC_PATH",
]
