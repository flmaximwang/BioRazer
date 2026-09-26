"""Stationary (object) representation of biological sequences.

This package is **the catalogue of the sequence objects biorazer works
with**: the ones biorazer defines itself, and the ones borrowed from other
packages.  Every other module imports sequence objects from here
(``from biorazer.sequence.objects import ProteinSequence``) rather than
importing biotite directly, so that *which* object is in use -- and *where it
comes from* -- is readable in one file.  It mirrors
:mod:`biorazer.structure.objects`, one layer down.

Biorazer's own objects
----------------------
None yet.  A sequence object of biorazer's own (an annotated sequence, an
alignment wrapper) belongs here when one exists; today the sequence layer
works directly on biotite's objects.

Borrowed objects
----------------
- :mod:`.bt_sequence` -- biotite ``ProteinSequence`` / ``NucleotideSequence``;
  the core sequence objects.
- :mod:`.bt_alphabet` -- biotite ``LetterAlphabet`` (the symbol set a sequence
  is defined over) and ``AlphabetError`` (raised for a letter outside it).
- :mod:`.bt_profile` -- biotite ``SequenceProfile``; per-column symbol counts.
- :mod:`.bt_align` -- biotite ``Alignment`` and ``SubstitutionMatrix``; the
  alignment object and the scoring matrix it is computed with.

What is *not* here
------------------
- **Helper functions of the provider** -- ``biotite.sequence.align.align_multiple``,
  ``biotite.sequence.graphics.plot_sequence_logo`` -- are not objects and stay
  imported from biotite where they are used.
- **File objects and their parsers/writers** live in
  :mod:`biorazer.sequence.io` (e.g. ``biotite.sequence.io.fasta.FastaFile``).
- **Conversions between these objects** live in
  :mod:`biorazer.sequence.bridge` (``central_dogma``, ``translation``) and in
  :mod:`biorazer.structure.bridge.sequence` (``AtomArray`` ->
  ``ProteinSequence``), as ``A_B`` converter classes.  The objects themselves
  carry no ``from_a`` / ``to_b`` methods.
"""

from biorazer.sequence.objects.bt_align import Alignment, SubstitutionMatrix
from biorazer.sequence.objects.bt_alphabet import AlphabetError, LetterAlphabet
from biorazer.sequence.objects.bt_profile import SequenceProfile
from biorazer.sequence.objects.bt_sequence import (
    NucleotideSequence,
    ProteinSequence,
)

__all__ = [
    # biotite sequence
    "ProteinSequence",
    "NucleotideSequence",
    # biotite alphabet
    "LetterAlphabet",
    "AlphabetError",
    # biotite profile
    "SequenceProfile",
    # biotite alignment
    "Alignment",
    "SubstitutionMatrix",
]
