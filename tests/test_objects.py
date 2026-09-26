"""Catalogue tests for ``biorazer.structure.objects`` / ``biorazer.sequence.objects``.

The convention (user-established) is that **every object biorazer works with is
enumerated once** in the catalogue of its layer -- ``biorazer.structure.objects``
for structure objects (its own and the ones borrowed from biotite / rdkit /
biopython / pyrosetta), ``biorazer.sequence.objects`` for the sequence objects
borrowed from biotite -- and that every other module imports them from there
instead of importing the provider package directly.  These tests pin both
halves of that:

* the catalogue re-exports the *same* classes as their provider packages
  (identity, not copies);
* importing the structure catalogue does not drag in an optional dependency;
* no module outside the catalogues references a borrowed object through the
  provider (the "one home per object name" rule), checked on the source tree.
"""

import ast
import os
import pathlib
import subprocess
import sys

import pytest

PACKAGE_ROOT = pathlib.Path(__file__).resolve().parents[1] / "biorazer"
CATALOGUES = (
    PACKAGE_ROOT / "structure" / "objects",
    PACKAGE_ROOT / "sequence" / "objects",
)
REPO_ROOT = PACKAGE_ROOT.parent

#: Which names of a provider module are *objects* (and must therefore be
#: imported from the catalogue).  The rest of each provider is deliberately not
#: listed: helper functions (``biotite.structure.sasa``,
#: ``biotite.sequence.align.align_multiple``), file objects (``PDBFile``,
#: ``FastaFile``) and parsers/writers (``PDBParser``, ``PDBIO``) are not objects
#: and stay imported from their provider where they are used.
BANNED_FROM_IMPORTS = {
    "biotite.structure": {"AtomArray", "AtomArrayStack", "BondList", "BondType"},
    "biotite.sequence": {
        "ProteinSequence", "NucleotideSequence", "SequenceProfile",
        "LetterAlphabet", "AlphabetError",
    },
    "biotite.sequence.profile": {"SequenceProfile"},
    "biotite.sequence.align": {"Alignment", "SubstitutionMatrix"},
    "rdkit.Chem": {"Mol"},
    "Bio.PDB": {"Structure", "Model"},
    "Bio.PDB.internal_coords": {"IC_Chain"},
}

#: Provider modules whose attribute access (``bio_struct.AtomArray``) is checked
#: too -- the name can be *reached* through a module alias instead of imported.
BANNED_MODULES = frozenset(BANNED_FROM_IMPORTS)
BANNED_ATTRIBUTES = frozenset(
    name for names in BANNED_FROM_IMPORTS.values() for name in names)


def _module_aliases(tree):
    """``{local name: dotted module}`` for the ``import`` statements in a file.

    Both forms are collected: ``import biotite.structure as bio_struct``
    (``ast.Import``) and ``from biotite import structure as bio_struct``
    (``ast.ImportFrom``, where the imported *name* is a module).
    """
    aliases = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                aliases[alias.asname or alias.name.split(".")[0]] = alias.name
        elif isinstance(node, ast.ImportFrom) and node.module:
            for alias in node.names:
                aliases[alias.asname or alias.name] = f"{node.module}.{alias.name}"
    return aliases


def _dotted(node):
    """``"a.b.c"`` for a (possibly nested) attribute expression, else ``None``."""
    parts = []
    while isinstance(node, ast.Attribute):
        parts.append(node.attr)
        node = node.value
    if isinstance(node, ast.Name):
        parts.append(node.id)
        return ".".join(reversed(parts))
    return None


def _borrowed_object_uses(path):
    """Lines of `path` that reach a borrowed object anywhere but the catalogue.

    Parsed with :mod:`ast`, so prose (docstrings, comments) that merely *names*
    a provider class -- e.g. "``Bio.PDB.Structure``" in a module docstring -- is
    not an offence; only imports and attribute access are.
    """
    tree = ast.parse(path.read_text(errors="replace"))
    aliases = _module_aliases(tree)
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            names = BANNED_FROM_IMPORTS.get(node.module or "", ())
            for alias in node.names:
                if alias.name in names:
                    yield (node.lineno,
                           f"from {node.module} import {alias.name}")
        elif isinstance(node, ast.Attribute) and node.attr in BANNED_ATTRIBUTES:
            dotted = _dotted(node)
            if dotted is None:
                continue
            root = dotted.split(".")[0]
            module = aliases.get(root)
            if module in BANNED_MODULES or ".".join(
                    dotted.split(".")[:-1]) in BANNED_MODULES:
                yield node.lineno, dotted


def test_public_names_reexported():
    """Every name each catalogue advertises is importable from its package."""
    from biorazer.sequence import objects as SO
    from biorazer.structure import objects as O

    for catalogue in (O, SO):
        for name in catalogue.__all__:
            assert hasattr(catalogue, name), f"missing public name: {name}"


def test_reexports_are_the_defining_objects():
    """The catalogue names ARE the provider classes (identity, not copies)."""
    import biotite.structure as bt
    from rdkit.Chem import Mol as RdMol

    from biorazer.structure import objects as O
    from biorazer.structure.objects import bp_icchain, bt_atom_array, bt_bond, rd_mol

    assert O.AtomArray is bt_atom_array.AtomArray is bt.AtomArray
    assert O.AtomArrayStack is bt_atom_array.AtomArrayStack is bt.AtomArrayStack
    assert O.BondList is bt_bond.BondList is bt.BondList
    assert O.BondType is bt_bond.BondType is bt.BondType
    assert O.Mol is rd_mol.Mol is RdMol
    assert O.IC_Chain is bp_icchain.IC_Chain
    assert O.InternalCoord is not bt.AtomArray      # biorazer's own object


def test_sequence_reexports_are_the_defining_objects():
    """Same identity rule for the sequence catalogue."""
    import biotite.sequence as bt_seq
    import biotite.sequence.align as bt_align

    from biorazer.sequence import objects as SO
    from biorazer.sequence.objects import (
        bt_align as seq_bt_align,
        bt_alphabet,
        bt_profile,
        bt_sequence,
    )

    assert SO.ProteinSequence is bt_sequence.ProteinSequence is bt_seq.ProteinSequence
    assert (SO.NucleotideSequence is bt_sequence.NucleotideSequence
            is bt_seq.NucleotideSequence)
    assert SO.SequenceProfile is bt_profile.SequenceProfile is bt_seq.SequenceProfile
    assert SO.LetterAlphabet is bt_alphabet.LetterAlphabet is bt_seq.LetterAlphabet
    assert SO.AlphabetError is bt_alphabet.AlphabetError is bt_seq.AlphabetError
    assert SO.Alignment is seq_bt_align.Alignment is bt_align.Alignment
    assert (SO.SubstitutionMatrix is seq_bt_align.SubstitutionMatrix
            is bt_align.SubstitutionMatrix)


def test_optional_dependency_is_not_imported():
    """The lazy PyRosetta accessor keeps pyrosetta out of ``import objects``."""
    code = (
        "import sys, biorazer.structure.objects as O;"
        "print('pyrosetta' in sys.modules);"
        "print(callable(O.pose_class))"
    )
    env = dict(os.environ)
    env["PYTHONPATH"] = str(REPO_ROOT) + os.pathsep + env.get("PYTHONPATH", "")
    proc = subprocess.run(
        [sys.executable, "-c", code],
        cwd=REPO_ROOT, env=env, capture_output=True, text=True, check=True,
    )
    assert proc.stdout.split() == ["False", "True"], proc.stdout


def test_bridge_converter_classes_are_reexported():
    """The bridge package exposes its ``A_B`` converters under one name each."""
    from biorazer.structure import bridge as B

    assert set(B.__all__) == {
        "AtomArray_InternalCoord",
        "InternalCoord_AtomArray",
        "SMCRA_ICChain",
        "ICChain_SMCRA",
        "AtomArray_ProteinSequence",
        "AtomArraySelection_AtomArrayMask",
        "AtomArraySelection_AtomArrayIndices",
        "RuleCsv_AtomArraySelection",
        "AtomArraySelection_RuleCsv",
        "SelectionCsv_AtomArraySelection",
        "AtomArraySelection_SelectionCsv",
        "AtomArraySelection_PyMOLSelection",
        "PyMOLSelection_AtomArraySelection",
        "PyMOLSelection_SelectionCsv",
        "SelectionCsv_PyMOLSelection",
        "Mask_PyMOLSelection",
        "Indices_PyMOLSelection",
    }
    for name in B.__all__:
        assert isinstance(getattr(B, name), type), name


#: Synthetic sources the guard must flag / must not flag -- the detector is
#: itself tested, so a guard that silently stopped matching is caught here.
OFFENDING_SOURCES = [
    "from biotite.structure import AtomArray\n",
    "from biotite.structure import AtomArrayStack as Stack\n",
    "import biotite.structure as bio_struct\nx = bio_struct.BondList(3)\n",
    "from biotite import structure as bio_struct\nx: bio_struct.AtomArray = None\n",
    "from biotite.sequence import ProteinSequence\n",
    "from biotite.sequence import ProteinSequence, SequenceProfile\n",
    "from biotite.sequence.profile import SequenceProfile\n",
    "from biotite.sequence.align import Alignment, SubstitutionMatrix\n",
    "import biotite.sequence as bio_seq\nx = bio_seq.LetterAlphabet('ACDE')\n",
    "from rdkit.Chem import Mol, SDWriter\n",
    "from Bio.PDB import Model as BioModel, Structure as BioStructure\n",
    "from Bio.PDB.internal_coords import IC_Chain\n",
    "import Bio.PDB\nx = Bio.PDB.Structure.Structure('x')\n",
]
ALLOWED_SOURCES = [
    '"""A docstring naming ``biotite.structure.AtomArray`` and Bio.PDB.Structure."""\n',
    "# a comment: from biotite.structure import AtomArray\n",
    "import biotite.structure as bio_struct\nx = bio_struct.sasa(structure)\n",
    "from biotite.structure.io import pdb, pdbx\n",
    "from biotite.structure.io.pdb.hybrid36 import encode_hybrid36\n",
    "from biotite.sequence.io import fasta\n",
    "from biotite.sequence.io.fasta import FastaFile\n",
    "from biotite.sequence.align import align_multiple\n",
    "from biotite.sequence.graphics import plot_sequence_logo\n",
    "from Bio.PDB import MMCIFIO, MMCIFParser, PDBIO, PDBParser\n",
    "from biorazer.structure.objects import AtomArray, BondList, Mol, IC_Chain\n",
    "from biorazer.sequence.objects import Alignment, ProteinSequence\n",
]


@pytest.mark.parametrize("source", OFFENDING_SOURCES)
def test_guard_flags_a_provider_reference(tmp_path, source):
    """The detector fires on every banned import / attribute form."""
    path = tmp_path / "sample.py"
    path.write_text(source)
    assert list(_borrowed_object_uses(path)), source


@pytest.mark.parametrize("source", ALLOWED_SOURCES)
def test_guard_allows_prose_helpers_and_the_catalogue(tmp_path, source):
    """Prose, provider helpers/parsers and catalogue imports are not offences."""
    path = tmp_path / "sample.py"
    path.write_text(source)
    assert list(_borrowed_object_uses(path)) == [], source


def test_borrowed_objects_are_referenced_only_from_the_catalogue():
    """No module outside ``objects/`` reaches a borrowed object via its provider."""

    def excluded(path):
        # The catalogues themselves, and ``sequence/archive`` -- frozen legacy
        # code kept for reference, never migrated (see the archive modules).
        return (any(catalogue in path.parents for catalogue in CATALOGUES)
                or "archive" in path.parts)

    offences = []
    for path in sorted(PACKAGE_ROOT.rglob("*.py")):
        if excluded(path):
            continue
        for lineno, what in _borrowed_object_uses(path):
            offences.append(f"{path.relative_to(PACKAGE_ROOT)}:{lineno}: {what}")
    assert not offences, (
        "borrowed objects reached outside biorazer/structure/objects / "
        "biorazer/sequence/objects (import them from there instead):\n"
        + "\n".join(offences))
