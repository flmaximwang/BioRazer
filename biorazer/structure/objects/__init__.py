"""Stationary (object) representation of biomolecular structures.

This package holds the user-defined "object" layer of
:mod:`biorazer.structure` -- representations of a structure as Python objects
as opposed to the converters in ``biorazer.structure.io`` and the
analytics / selection / manipulation modules that operate on them.

Third-party structure objects (e.g. biotite's ``AtomArray``) are re-exported
once, in :mod:`biorazer.structure.objects.external`, and the modules in this
package build on those instead of importing biotite directly.
"""

from biorazer.structure.objects.external import AtomArray, AtomArrayStack
from biorazer.structure.objects.internal_coords import (
    InternalCoordAtom,
    InternalCoord,
)

__all__ = [
    "AtomArray",
    "AtomArrayStack",
    "InternalCoordAtom",
    "InternalCoord",
]