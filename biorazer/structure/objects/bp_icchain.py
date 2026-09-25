"""biopython's structure objects, re-exported for the object layer.

Two representations of the same molecule are enumerated here:

- the **SMCRA hierarchy** roots -- ``Bio.PDB.Structure`` and ``Bio.PDB.Model``
  (Structure-Model-Chain-Residue-Atom): biopython's tree of nested Python
  objects, one ``Chain`` object per chain, each holding ``Residue`` and
  ``Atom`` objects;
- ``Bio.PDB.internal_coords.IC_Chain`` -- biopython's *internal coordinate*
  (Z-matrix / torsion) representation of one chain, as produced by
  ``IC_Chain.atom_to_internal_coordinates()``.

Both are converted to/from each other by
:mod:`biorazer.structure.bridge.icchain`, and to/from files by the converter
classes in :mod:`biorazer.structure.io.protein`.  biopython is a core
dependency, so both import at module level.

Naming note
-----------
The class named ``InternalCoord`` in
:mod:`biorazer.structure.objects.internal_coords` is **biorazer's own** object.
biopython's internal-coordinate object is ``IC_Chain`` -- here and everywhere
else in the package -- so that ``InternalCoord*`` names never refer to
biopython's class.

Import paths
------------
``Structure`` and ``Model`` are the **classes**, imported from their
submodules (``Bio.PDB.Structure.Structure`` / ``Bio.PDB.Model.Model``).  Note
that ``from Bio.PDB import Structure`` would bind the *submodule* of that name,
not the class, which is why the submodule paths are used here.
"""

from Bio.PDB.internal_coords import IC_Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure

__all__ = [
    "Structure",
    "Model",
    "IC_Chain",
]
