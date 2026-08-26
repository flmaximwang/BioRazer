"""File-format glue for the biopython internal-coordinate converters.

The SMCRA <-> internal-coordinate bridge lives in
``biorazer.structure.bridge.internal_coords`` (:class:`SMCRA_InternalCoord` /
:class:`InternalCoord_SMCRA`). This module only re-exports those bridge
converters plus the Bio.PDB file parsers/writers that the protein package's
internal-coordinate converters (Pdb_InternalCoord / Cif_InternalCoord /
InternalCoord_Pdb / InternalCoord_Cif) need to move structures to/from files.
"""

from Bio.PDB import MMCIFIO, MMCIFParser, PDBIO, PDBParser

from biorazer.structure.bridge.internal_coords import (
    InternalCoord_SMCRA,
    SMCRA_InternalCoord,
)

__all__ = [
    "MMCIFIO",
    "MMCIFParser",
    "PDBIO",
    "PDBParser",
    "SMCRA_InternalCoord",
    "InternalCoord_SMCRA",
]
