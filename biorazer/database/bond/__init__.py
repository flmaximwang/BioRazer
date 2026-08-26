# -*- coding: utf-8 -*-
"""Covalent bond geometry for amino acids, split by backbone / side-chain.

Imports resolve as ``biorazer.database.bond` (package).  Backbone bond
lengths and angles live in :mod:`.backbone`; side-chain bonds/angles (atoms
beyond CA-CB) in :mod:`.sidechain`.  Backward-compatible re-exports of the old
``biorazer.database.bond`` module's names are provided here.
"""

from .backbone import (  # noqa: F401
    AMINO_ACID_BOND_LENGTH,
    AMINO_ACID_BOND_LENGTH_BY_RESIDUE,
    AMINO_ACID_BOND_ANGLE,
    AMINO_ACID_BOND_ANGLE_BY_RESIDUE,
    BOND_REFS,
)
from .sidechain import (  # noqa: F401
    AMINO_ACID_SIDECHAIN_BOND,
    BOND_SIDECHAIN_REFS,
)