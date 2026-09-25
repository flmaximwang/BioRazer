"""Bridges: in-memory converters between two object representations.

Each module holds the ``A_B`` :class:`~biorazer.io.Converter` classes that move
an in-memory source object (``A``) to an in-memory target object (``B``) -- no
files involved, that is :mod:`biorazer.structure.io`.  The object classes
themselves carry no ``from_a`` / ``to_b`` methods (see
:mod:`biorazer.structure.objects`).

- :mod:`~biorazer.structure.bridge.atom_array` -- biotite ``AtomArray`` <->
  biorazer ``InternalCoord``.
- :mod:`~biorazer.structure.bridge.internal_coords` -- biopython ``Structure``
  (SMCRA) <-> biopython ``IC_Chain``.
"""

from biorazer.structure.bridge.atom_array import (
    AtomArray_InternalCoord,
    InternalCoord_AtomArray,
)

__all__ = [
    "AtomArray_InternalCoord",
    "InternalCoord_AtomArray",
]
