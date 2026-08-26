"""Amino-acid alphabets and code maps.

Re-export shim over the submodules so ``from biorazer.database.alphabet import ...``
works for every public name (moved from ``biorazer.database.amino_acid``).
"""

from .aa_list import (
    AMINO_ACIDS_1LETTER,
    AMINO_ACIDS_1TO3_INITIAL_CAPITAL,
    AMINO_ACIDS_1TO3_LOWER,
    AMINO_ACIDS_1TO3_UPPER,
    AMINO_ACIDS_3LETTER,
    AMINO_ACIDS_3TO1_LOWER,
    AMINO_ACIDS_3TO1_UPPER,
)
from .aa_types import TYPES2AA
from .protein import (
    APOLAR_AA_ALPHABET,
    APOLAR_AA_NAME1,
    COMMON_AA_ALPHABET,
    COMMON_AA_NAME1,
    HYDROPHOBIC_AA_ALPHABET,
    HYDROPHOBIC_AA_NAME1,
    NEGATIVE_AA_ALPHABET,
    NEGATIVE_AA_NAME1,
    POLAR_AA_ALPHABET,
    POLAR_AA_NAME1,
    POSITIVE_AA_ALPHABET,
    POSITIVE_AA_NAME1,
    sequences_to_symbols,
)

__all__ = [
    "AMINO_ACIDS_1LETTER",
    "AMINO_ACIDS_3LETTER",
    "AMINO_ACIDS_1TO3_UPPER",
    "AMINO_ACIDS_1TO3_LOWER",
    "AMINO_ACIDS_1TO3_INITIAL_CAPITAL",
    "AMINO_ACIDS_3TO1_UPPER",
    "AMINO_ACIDS_3TO1_LOWER",
    "TYPES2AA",
    "COMMON_AA_NAME1",
    "COMMON_AA_ALPHABET",
    "HYDROPHOBIC_AA_NAME1",
    "HYDROPHOBIC_AA_ALPHABET",
    "POLAR_AA_NAME1",
    "POLAR_AA_ALPHABET",
    "POSITIVE_AA_NAME1",
    "POSITIVE_AA_ALPHABET",
    "NEGATIVE_AA_NAME1",
    "NEGATIVE_AA_ALPHABET",
    "APOLAR_AA_NAME1",
    "APOLAR_AA_ALPHABET",
    "sequences_to_symbols",
]
