"""File-format glue for the biopython IC_Chain converters.

The SMCRA <-> IC_Chain bridge lives in
``biorazer.structure.bridge.icchain`` (:class:`SMCRA_ICChain` /
:class:`ICChain_SMCRA`).  This module only re-exports those bridge converters
plus the Bio.PDB file parsers/writers that the protein package's IC_Chain
converters (Pdb_ICChain / Cif_ICChain / ICChain_Pdb / ICChain_Cif) need to move
structures to/from files.

The parsers/writers (``PDBParser``, ``PDBIO``, ``MMCIFParser``, ``MMCIFIO``)
are I/O tools, not structure objects, so they stay imported from biopython
here; the structure objects they produce (``Structure``, ``Model``) come from
:mod:`biorazer.structure.objects`.
"""

from Bio.PDB import MMCIFIO, MMCIFParser, PDBIO, PDBParser

from biorazer.structure.bridge.icchain import ICChain_SMCRA, SMCRA_ICChain

__all__ = [
    "MMCIFIO",
    "MMCIFParser",
    "PDBIO",
    "PDBParser",
    "SMCRA_ICChain",
    "ICChain_SMCRA",
]
