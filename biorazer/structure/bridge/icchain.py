"""Converters between biopython's SMCRA hierarchy and its internal coordinates.

SMCRA (Structure-Model-Chain-Residue-Atom) is the hierarchy of
``Bio.PDB.Structure``.  This module bridges that in-memory hierarchy to
biopython's internal-coordinate (torsion) representation:

- :class:`SMCRA_ICChain` -- ``Bio.PDB.Structure`` -> list of
  :class:`Bio.PDB.internal_coords.IC_Chain` (one per chain of the first
  model), with ``atom_to_internal_coordinates`` already applied.
- :class:`ICChain_SMCRA` -- list of ``IC_Chain`` -> ``Bio.PDB.Structure``
  (Cartesian coords regenerated from the internal coords).

Both are pure in-memory bridges, so they implement the parent's
:meth:`~biorazer.io.Converter.convert` (a single transform, no ``read()`` /
``write()``): the source object is held in ``input_io`` and ``convert()``
returns the result.

Naming
------
``IC_Chain`` is **biopython's** internal-coordinate object.  biorazer's own
internal-coordinate object is ``InternalCoord`` (see
:mod:`biorazer.structure.objects.internal_coords`), so every converter in this
package named ``InternalCoord_*`` / ``*_InternalCoord`` moves *biorazer's*
object -- this module's pair, which moves biopython's, is named after
``IC_Chain``.

The object classes come from :mod:`biorazer.structure.objects` rather than
from biopython directly; biopython is a core dependency, so they import at
module level.
"""

from typing import Any, cast

from biorazer.io import Converter
from biorazer.structure.objects import IC_Chain, Model, Structure


class SMCRA_ICChain(Converter):
    """
    Converts a ``Bio.PDB.Structure`` (SMCRA hierarchy) to biopython's
    internal-coordinate representation.

    ``input_io`` holds a ``Bio.PDB.Structure``; :meth:`convert` returns a list
    of :class:`Bio.PDB.internal_coords.IC_Chain`, one per chain of the first
    model, each with its Cartesian coordinates already converted to
    internal/torsion coordinates via ``atom_to_internal_coordinates``.
    """

    def convert(self, **kwargs):
        structure = cast(Structure, self.input_io)
        chains = []
        model = next(iter(structure))
        for chain in model:
            ic = IC_Chain(chain)
            ic.atom_to_internal_coordinates()
            chains.append(ic)
        return chains


class ICChain_SMCRA(Converter):
    """
    Regenerates a ``Bio.PDB.Structure`` (SMCRA hierarchy) from biopython's
    internal-coordinate representation (the inverse of
    :class:`SMCRA_ICChain`).

    ``input_io`` holds an :class:`Bio.PDB.internal_coords.IC_Chain` or a list
    of them; :meth:`convert` rebuilds each chain back to Cartesian coordinates
    from its internal/torsion coordinates via
    ``internal_to_atom_coordinates`` and returns the rebuilt structure.
    """

    def convert(self, **kwargs):
        chains: Any = self.input_io
        if not isinstance(chains, (list, tuple)):
            chains = [chains]
        for ic in chains:
            ic.internal_to_atom_coordinates()
        structure = Structure("biorazer")
        model = Model(0)
        for ic in chains:
            chain = ic.chain
            if chain.parent is not None:
                chain.parent.detach_child(chain.id)
            model.add(chain)
        structure.add(model)
        return structure
