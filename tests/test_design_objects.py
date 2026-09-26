"""``biorazer.design.objects`` catalogue tests.

The convention (user-established, same as ``biorazer.structure.objects``) is
that **every design object biorazer defines is enumerated once** in
``biorazer.design.objects``, and that every other module imports them from
there (``from biorazer.design.objects import Entry``) instead of from the
module that happens to define them.  These tests pin both halves of that:

* the catalogue advertises every class its modules define -- a class added to
  ``objects/`` but not to ``__all__`` fails here;
* the advertised names ARE the defining classes (identity, not copies), and
  the inheritance chain the catalogue documents is the real one.
"""

import ast
import pathlib

import biorazer.design
from biorazer.design import objects as O

CATALOGUE = pathlib.Path(biorazer.design.__file__).resolve().parent / "objects"


def _top_level_classes(path):
    """Names of the classes defined at the top level of `path`."""
    tree = ast.parse(path.read_text(errors="replace"))
    return [node.name for node in tree.body if isinstance(node, ast.ClassDef)]


def test_public_names_reexported():
    """Every name the catalogue advertises is importable from the package."""
    for name in O.__all__:
        assert hasattr(O, name), f"missing public name: {name}"


def test_every_class_defined_here_is_advertised():
    """A class defined in ``objects/`` but left out of ``__all__`` is a bug."""
    defined = set()
    for path in sorted(CATALOGUE.glob("*.py")):
        if path.name.startswith("_"):  # _shared.py: private helpers, not objects
            continue
        defined.update(_top_level_classes(path))

    assert defined, "no classes found -- the detector is looking in the wrong place"
    assert defined <= set(O.__all__), sorted(defined - set(O.__all__))


def test_reexports_are_the_defining_objects():
    """The catalogue names ARE the classes of the ``objects`` submodules."""
    from biorazer.design.objects import entry, library, sequence, single_test

    assert O.SingleTest is single_test.SingleTest
    assert O.EntryProperty is entry.EntryProperty
    assert O.EntryFormatter is entry.EntryFormatter
    assert O.EntryPropertyFormatted is entry.EntryPropertyFormatted
    assert O.EntryIO is entry.EntryIO
    assert O.Entry is entry.Entry
    assert O.Library is library.Library
    assert O.SequenceEntry is sequence.SequenceEntry
    assert O.SequenceLibrary is sequence.SequenceLibrary


def test_inheritance_chain_matches_the_catalogue():
    """Property -> Formatter -> PropertyFormatted -> IO -> Entry is real."""
    assert issubclass(O.EntryFormatter, O.EntryProperty)
    assert issubclass(O.EntryPropertyFormatted, O.EntryFormatter)
    assert issubclass(O.EntryIO, O.EntryPropertyFormatted)
    assert issubclass(O.Entry, O.EntryIO)
    assert issubclass(O.SequenceEntry, O.Entry)
    assert issubclass(O.SequenceLibrary, O.Library)
    assert O.SequenceLibrary.entry_type is O.SequenceEntry
