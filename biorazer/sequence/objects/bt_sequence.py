"""biotite's sequence objects, re-exported for the object layer.

``ProteinSequence`` and ``NucleotideSequence`` are biorazer's core sequence
objects (see the :mod:`biorazer.sequence` docstring): a sequence of letters
over a fixed alphabet, holding no annotation and no coordinates.  They are the
counterpart of :class:`~biorazer.structure.objects.AtomArray` one layer up --
the structure layer resolves a residue to 3-D coordinates, this layer to a
single letter.

Both are imported once, here, and composed by the rest of the package: every
other module imports them from ``biorazer.sequence.objects`` instead of from
biotite directly, so which object is in use is readable in one file.
"""

from biotite.sequence import NucleotideSequence, ProteinSequence

__all__ = [
    "ProteinSequence",
    "NucleotideSequence",
]
