"""Alphabet package refactor identity tests.

``biorazer/database/amino_acid.py`` was split into the
``biorazer/database/alphabet`` package (``aa_list`` / ``aa_types`` /
``protein`` submodules) with a re-export shim ``__init__``.  These tests pin
the refactor: every public name is still importable from the package, the
re-exports are the SAME objects as the defining submodules, and the old module
path is gone.
"""

import importlib

import pytest


def test_public_names_reexported():
    """Every name the package exposes is importable from the package path."""
    from biorazer.database import alphabet as A

    for name in A.__all__:
        assert hasattr(A, name), f"missing public name: {name}"


def test_reexports_are_the_defining_objects():
    """The package-level names ARE the submodule objects (identity, not copies)."""
    from biorazer.database import alphabet as A
    from biorazer.database.alphabet import aa_list, aa_types, protein

    assert A.AMINO_ACIDS_1LETTER is aa_list.AMINO_ACIDS_1LETTER
    assert A.AMINO_ACIDS_3LETTER is aa_list.AMINO_ACIDS_3LETTER
    assert A.AMINO_ACIDS_1TO3_UPPER is aa_list.AMINO_ACIDS_1TO3_UPPER
    assert A.AMINO_ACIDS_1TO3_LOWER is aa_list.AMINO_ACIDS_1TO3_LOWER
    assert A.AMINO_ACIDS_1TO3_INITIAL_CAPITAL is aa_list.AMINO_ACIDS_1TO3_INITIAL_CAPITAL
    assert A.AMINO_ACIDS_3TO1_UPPER is aa_list.AMINO_ACIDS_3TO1_UPPER
    assert A.AMINO_ACIDS_3TO1_LOWER is aa_list.AMINO_ACIDS_3TO1_LOWER
    assert A.TYPES2AA is aa_types.TYPES2AA
    assert A.COMMON_AA_ALPHABET is protein.COMMON_AA_ALPHABET
    assert A.sequences_to_symbols is protein.sequences_to_symbols


def test_old_module_path_gone():
    """The pre-refactor module path no longer imports."""
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("biorazer.database.amino_acid")


def test_maps_values():
    """Spot-check the data values survived the split unchanged."""
    from biorazer.database.alphabet import (
        AMINO_ACIDS_1LETTER,
        AMINO_ACIDS_1TO3_UPPER,
        AMINO_ACIDS_3TO1_UPPER,
        TYPES2AA,
    )

    assert AMINO_ACIDS_1LETTER == tuple("ACDEFGHIKLMNPQRSTVWY")
    assert AMINO_ACIDS_1TO3_UPPER["M"] == "MET"
    assert AMINO_ACIDS_3TO1_UPPER["Ala"] == "A"
    assert TYPES2AA["hydrophobic"] == "AFILMPVWY"
