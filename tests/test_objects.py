"""``biorazer.structure.objects`` catalogue tests.

The convention (user-established) is that **every structure object biorazer
works with is enumerated once** in ``biorazer.structure.objects`` -- its own
objects and the ones borrowed from biotite / rdkit / biopython / pyrosetta --
and that every other module imports them from there instead of importing the
provider package directly.  These tests pin both halves of that:

* the catalogue re-exports the *same* classes as their provider packages
  (identity, not copies);
* importing the catalogue does not drag in an optional dependency;
* no module outside the catalogue references a borrowed object through the
  provider (the "one home per object name" rule), checked on the source tree.
"""

import os
import pathlib
import re
import subprocess
import sys

PACKAGE_ROOT = pathlib.Path(__file__).resolve().parents[1] / "biorazer"
CATALOGUE = PACKAGE_ROOT / "structure" / "objects"
REPO_ROOT = PACKAGE_ROOT.parent

#: ``(regex, why)`` -- a borrowed **object** referenced outside the catalogue.
#: Biorazer's *functions* of those packages are not listed: helpers such as
#: ``biostruc.sasa`` or ``pdb.PDBFile`` are not structure objects and stay
#: imported from their provider where they are used.
BANNED_OBJECT_REFERENCES = [
    (
        r"from biotite\.structure import (?:AtomArray|AtomArrayStack|BondList|BondType)\b",
        "import it from biorazer.structure.objects instead",
    ),
    (
        r"\b(?:bio_struct|bio_struc|bt_struct|bt_struc)"
        r"\.(?:AtomArray|AtomArrayStack|BondList|BondType)\b",
        "use the name imported from biorazer.structure.objects, "
        "not the biotite module alias",
    ),
    (
        r"from rdkit\.Chem import [^#\n]*\bMol\b",
        "import Mol from biorazer.structure.objects instead",
    ),
]


def test_public_names_reexported():
    """Every name the catalogue advertises is importable from the package."""
    from biorazer.structure import objects as O

    for name in O.__all__:
        assert hasattr(O, name), f"missing public name: {name}"


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


def test_borrowed_objects_are_referenced_only_from_the_catalogue():
    """No module outside ``objects/`` names a borrowed object via its provider."""
    offences = []
    for path in sorted(PACKAGE_ROOT.rglob("*.py")):
        if CATALOGUE in path.parents:
            continue
        for lineno, line in enumerate(
                path.read_text(errors="replace").splitlines(), start=1):
            for pattern, why in BANNED_OBJECT_REFERENCES:
                if re.search(pattern, line):
                    offences.append(
                        f"{path.relative_to(PACKAGE_ROOT)}:{lineno}: "
                        f"{line.strip()}  -- {why}")
    assert not offences, (
        "borrowed structure objects referenced outside "
        "biorazer/structure/objects:\n" + "\n".join(offences))
