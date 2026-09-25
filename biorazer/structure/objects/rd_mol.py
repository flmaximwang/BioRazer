"""RDKit's ``Mol``, re-exported for the object layer.

A ``Mol`` is the small-molecule counterpart of biotite's
:class:`~biorazer.structure.objects.AtomArray`: RDKit's own in-memory molecule,
carrying explicit atoms, bonds, bond orders and formal charges (and, for
docking-style work, conformers).  The organic-molecule converters in
:mod:`biorazer.structure.io.organic` read into and write from it, so it is
enumerated here with the other structure objects.

RDKit's *writers* (``SDWriter``) and the ``rdkit.Chem`` module functions are
I/O helpers, not objects; they stay imported from rdkit where they are used.
"""

from rdkit.Chem import Mol

__all__ = [
    "Mol",
]
