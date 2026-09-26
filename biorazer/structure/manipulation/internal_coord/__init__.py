"""Manipulation of biorazer's own internal-coordinate container.

The operand here is ``InternalCoord``
(:mod:`biorazer.structure.objects.internal_coords`), not biotite's
``AtomArray`` -- the AtomArray operations live in the sibling package
:mod:`biorazer.structure.manipulation.atom_array`.

- :mod:`~biorazer.structure.manipulation.internal_coord.modification` --
  connect two ``InternalCoord`` fragments with a new peptide bond.
- :mod:`~biorazer.structure.manipulation.internal_coord.mutation` --
  swap a residue's side chain for another residue's template.
"""

from biorazer.structure.manipulation.internal_coord.modification import (
    connect_internal_coords,
)
from biorazer.structure.manipulation.internal_coord.mutation import (
    mutate,
)

__all__ = [
    "connect_internal_coords",
    "mutate",
]
