"""Bridges: in-memory converters between two object representations.

Each module holds the ``A_B`` :class:`~biorazer.io.Converter` classes that move
an in-memory source object (``A``) to an in-memory target object (``B``) -- no
files involved, that is :mod:`biorazer.structure.io`.  The object classes
themselves carry no ``from_a`` / ``to_b`` methods (see
:mod:`biorazer.structure.objects`).

- :mod:`~biorazer.structure.bridge.atom_array` -- biotite ``AtomArray`` <->
  biorazer ``InternalCoord``.
- :mod:`~biorazer.structure.bridge.icchain` -- biopython ``Structure``
  (SMCRA) <-> biopython ``IC_Chain``.
- :mod:`~biorazer.structure.bridge.sequence` -- biotite ``AtomArray`` ->
  biotite ``ProteinSequence``.
- :mod:`~biorazer.structure.bridge.selector` -- biorazer ``AtomArraySelection``:
  its rule table / selection table csv (file-backed) and the boolean mask / atom
  indices it makes over an ``AtomArray`` (a *selection*, so the target array is
  an argument of ``convert()`` / ``write()``).
"""

from biorazer.structure.bridge.atom_array import (
    AtomArray_InternalCoord,
    InternalCoord_AtomArray,
)
from biorazer.structure.bridge.icchain import ICChain_SMCRA, SMCRA_ICChain
from biorazer.structure.bridge.selector import (
    AtomArraySelection_AtomArrayIndices,
    AtomArraySelection_AtomArrayMask,
    AtomArraySelection_RuleCsv,
    AtomArraySelection_SelectionCsv,
    RuleCsv_AtomArraySelection,
    SelectionCsv_AtomArraySelection,
)
from biorazer.structure.bridge.sequence import AtomArray_ProteinSequence

__all__ = [
    "AtomArray_InternalCoord",
    "InternalCoord_AtomArray",
    "SMCRA_ICChain",
    "ICChain_SMCRA",
    "AtomArray_ProteinSequence",
    "RuleCsv_AtomArraySelection",
    "AtomArraySelection_RuleCsv",
    "SelectionCsv_AtomArraySelection",
    "AtomArraySelection_SelectionCsv",
    "AtomArraySelection_AtomArrayMask",
    "AtomArraySelection_AtomArrayIndices",
]
