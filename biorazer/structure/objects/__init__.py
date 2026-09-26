"""Stationary (object) representation of biomolecular structures.

This package is **the catalogue of the structure objects biorazer works
with**: the ones biorazer defines itself, and the ones borrowed from other
packages.  Every other module imports structure objects from here
(``from biorazer.structure.objects import AtomArray``) rather than importing
biotite / rdkit / biopython / pyrosetta directly, so that *which* object is in
use -- and *where it comes from* -- is readable in one file.

Biorazer's own objects
----------------------
- :mod:`.internal_coords` -- :class:`InternalCoord` (a generative
  internal-coordinate description of a structure) and
  :class:`InternalCoordAtom` (its per-atom annotation record).
- :mod:`.selector` -- :class:`AtomArraySelection`: a set of per-atom selection rules
  over ``ins_code`` / ``chain`` / ``resi`` / ``name`` / ``altloc``, with its rule
  table / selection table csv (read and written by
  :mod:`biorazer.structure.bridge`) and the tkinter editor that edits them.

Borrowed objects
----------------
- :mod:`.bt_atom_array` -- biotite ``AtomArray`` / ``AtomArrayStack``; the
  core structure object.
- :mod:`.bt_bond` -- biotite ``BondList`` / ``BondType``; the bond table an
  ``AtomArray`` carries.
- :mod:`.rd_mol` -- RDKit ``Mol``; small molecules.
- :mod:`.bp_icchain` -- biopython ``Structure`` / ``Model`` (the SMCRA
  hierarchy) and ``IC_Chain`` (biopython's internal coordinates -- **not**
  biorazer's ``InternalCoord``).
- :mod:`.pr_pose` -- PyRosetta ``Pose``; optional dependency, exposed through
  the lazy accessor :func:`pose_class`.

What is *not* here
------------------
- **Conversions between these objects** live in
  :mod:`biorazer.structure.bridge`, as ``A_B`` converter classes (source ->
  target, e.g. ``AtomArray_InternalCoord``).  The objects themselves carry no
  ``from_a`` / ``to_b`` methods.
- **File I/O** lives in :mod:`biorazer.structure.io` (parsers, writers and the
  file-object types they need, e.g. ``biotite.structure.io.pdb.PDBFile``).
"""

from biorazer.structure.objects.bp_icchain import IC_Chain, Model, Structure
from biorazer.structure.objects.bt_atom_array import AtomArray, AtomArrayStack
from biorazer.structure.objects.bt_bond import BondList, BondType
from biorazer.structure.objects.internal_coords import (
    InternalCoordAtom,
    InternalCoord,
)
from biorazer.structure.objects.pr_pose import pose_class
from biorazer.structure.objects.rd_mol import Mol
from biorazer.structure.objects.selector import AtomArraySelection

__all__ = [
    # biorazer's own
    "InternalCoord",
    "InternalCoordAtom",
    "AtomArraySelection",
    # biotite
    "AtomArray",
    "AtomArrayStack",
    "BondList",
    "BondType",
    # rdkit
    "Mol",
    # biopython
    "Structure",
    "Model",
    "IC_Chain",
    # pyrosetta (lazy: optional dependency)
    "pose_class",
]
