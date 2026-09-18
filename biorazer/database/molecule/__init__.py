# -*- coding: utf-8 -*-
"""Per-molecule physico-chemical reference data, split by property.

``biorazer.database.molecule`` groups all molecule-level geometry / property
tables under one package so future molecule classes (nucleic acids, ligands,
...) have room to grow alongside the protein data:

* :mod:`.atom`    -- per-atom properties (van der Waals radius, charge).
* :mod:`.bond`    -- covalent bond geometry: length / angle / dihedral.
* :mod:`.icoor`   -- internal-coordinate reference data (grow-path topology
  and ideal per-residue templates for proteins).
* :mod:`.rotamer` -- **readers for external** rotamer libraries (PyMOL's
  bundled Dunbrack pickles, Rosetta's Shapovalov text libraries).  Those
  libraries are large and variously licensed, so they are not vendored
  here; the readers load them from the local installation.

Within :mod:`.bond`, each geometry type splits into ``generic``
(molecule-agnostic main tables) and ``protein`` (protein-residue-specific
refinements).  Every numeric entry carries the uniform record
``{mean, std, lb, up, source}`` (std/lb/up = ``np.nan`` when the source
provides no spread).

Convenience re-exports (the full public names also resolve from their
per-module homes)::

    from biorazer.database.molecule import (
        AMINO_ACID_BOND_LENGTH,            # bond.length.protein (by_residue)
        AMINO_ACID_BOND_ANGLE,             # bond.angle.protein (residue-keyed)
        AMINO_ACID_BACKBONE_BOND_ANGLE,    # bond.angle.generic (flat backbone)
        SS_BB_TORSION_ANGLE,               # bond.dihedral.protein.by_ss
        IC_PATH,                           # icoor.protein.topology
        BACKBONE_IC_PATH,                  # icoor.protein.topology
        ATOM_RADIUS,                       # atom.radius
        read_pymol_dep,                    # rotamer.pymol (external readers)
        read_shapovalov,                   # rotamer.rosetta
        RotamerRecord,                     # rotamer (shared record types)
        RotamerLibrary,
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
from .bond.angle.generic import AMINO_ACID_BACKBONE_BOND_ANGLE  # noqa: F401
from .bond.angle.protein import AMINO_ACID_BOND_ANGLE  # noqa: F401
from .bond.dihedral.protein import (  # noqa: F401
    SS_BB_TORSION_ANGLE,
    DSSP_SS_CODE,
    BB_TORSION_TURNS,
    BB_TORSION_REFS,
    OMEGA_TRANS,
    OMEGA_CIS,
    ALIAS_QUAD,
    QUAD_ALIAS,
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

# external rotamer library readers (split by source)
from .rotamer import (  # noqa: F401
    RotamerRecord,
    RotamerLibrary,
    N_MAINCHAIN as ROTAMER_N_MAINCHAIN,
    DEFAULT_BIN_GRID,
)
from .rotamer.rosetta import (  # noqa: F401
    ROSETTA_ROTAMER_DIR,
    DEFAULT_STEPDOWN,
    read_rosetta_text,
    read_shapovalov,
    read_bbdep02,
    default_bbdep02_path,
    default_shapovalov_path,
)
from .rotamer.pymol import (  # noqa: F401
    PYMOL_ROTAMER_DIR,
    PYMOL_SC_BB_IND,
    PYMOL_SC_BB_DEP,
    PYMOL_SC_LIBRARY,
    read_pymol_ind,
    read_pymol_dep,
    read_pymol_library,
)

__all__ = [
    # bond length
    "AMINO_ACID_BOND_LENGTH", "AMINO_ACID_BOND_LENGTH_BY_RESIDUE",
    "AMINO_ACID_SIDECHAIN_BOND", "BOND_REFS", "BOND_SIDECHAIN_REFS",
    # bond angle
    "AMINO_ACID_BOND_ANGLE", "AMINO_ACID_BACKBONE_BOND_ANGLE",
    # bond dihedral
    "SS_BB_TORSION_ANGLE", "DSSP_SS_CODE", "BB_TORSION_TURNS",
    "BB_TORSION_REFS", "OMEGA_TRANS", "OMEGA_CIS",
    "ALIAS_QUAD", "QUAD_ALIAS",
    "SIDECHAIN_CHI", "ROTAMER_BIN",
    "NON_ROTAMERIC_BIN_WIDTH",
    "SIDECHAIN_ROTAMER_LIB", "SIDECHAIN_NON_ROTAMERIC_BINS",
    "SIDECHAIN_DIHE_REFS", "AAS",
    # atom properties
    "ATOM_RADIUS", "vdw_dict", "vdw_radii",
    # icoor topology
    "IC_PATH", "MAINCHAIN_ATOMS", "BACKBONE_IC_PATH",
    # external rotamer readers -- shared record types
    "RotamerRecord", "RotamerLibrary", "ROTAMER_N_MAINCHAIN", "DEFAULT_BIN_GRID",
    # external rotamer readers -- Rosetta / Dunbrack
    "ROSETTA_ROTAMER_DIR", "DEFAULT_STEPDOWN",
    "read_rosetta_text", "read_shapovalov", "read_bbdep02",
    "default_bbdep02_path", "default_shapovalov_path",
    # external rotamer readers -- PyMOL
    "PYMOL_ROTAMER_DIR", "PYMOL_SC_BB_IND", "PYMOL_SC_BB_DEP", "PYMOL_SC_LIBRARY",
    "read_pymol_ind", "read_pymol_dep", "read_pymol_library",
]
