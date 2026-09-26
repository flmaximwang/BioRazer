"""Stationary (object) representation of protein design results.

This package is **the catalogue of the design objects biorazer works with**,
and it is the only place they are enumerated.  Every other module imports them
from here (``from biorazer.design.objects import Entry``) rather than from the
module that happens to define them, so that *which* object is in use is
readable in one file.  It mirrors :mod:`biorazer.structure.objects`.

Biorazer's own objects
----------------------
- :mod:`.single_test` -- :class:`SingleTest`; the results of one test run
  (one design, one seed) plus the metadata files that go with it.
- :mod:`.entry` -- the :class:`Entry` hierarchy, one subclass per step of the
  pipeline: :class:`EntryProperty` (what every entry carries) ->
  :class:`EntryFormatter` (how a program's output is reformatted into an
  entry directory) -> :class:`EntryPropertyFormatted` ->
  :class:`EntryIO` (the directory / DataFrame entry points) ->
  :class:`Entry` (one design).
- :mod:`.library` -- :class:`Library`; a collection of entries with the
  DataFrame merge and the plotting helpers.
- :mod:`.sequence` -- :class:`SequenceEntry` / :class:`SequenceLibrary`; the
  entry and library for designs that are a bare sequence.

What is *not* here
------------------
- **Directory and DataFrame conversion** -- ``Entry.from_dir``,
  ``Entry.from_dataframe``, ``Entry.to_dir``, ``SingleTest.from_dir`` -- is
  still a classmethod on the objects themselves.
- ``._shared`` is a private helper module (``_normalize_dir``) and is
  deliberately not re-exported.
"""

from biorazer.design.objects.entry import (
    Entry,
    EntryFormatter,
    EntryIO,
    EntryProperty,
    EntryPropertyFormatted,
)
from biorazer.design.objects.library import Library
from biorazer.design.objects.sequence import SequenceEntry, SequenceLibrary
from biorazer.design.objects.single_test import SingleTest

__all__ = [
    # one test run (one design, one seed)
    "SingleTest",
    # the Entry hierarchy, base -> concrete
    "EntryProperty",
    "EntryFormatter",
    "EntryPropertyFormatted",
    "EntryIO",
    "Entry",
    # a collection of entries
    "Library",
    # sequence-specialised entry / library
    "SequenceEntry",
    "SequenceLibrary",
]
