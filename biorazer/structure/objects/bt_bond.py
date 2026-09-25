"""biotite's bond objects, re-exported for the object layer.

``BondList`` is the per-structure bond table an ``AtomArray`` carries as its
``bonds`` annotation (``array.bonds``): pairs of atom indices plus a bond type
per row, queried with ``as_array()``.  ``BondType`` is the enum of bond orders
in that table (``BondType.SINGLE`` / ``DOUBLE`` / ``AROMATIC`` ...).

Both are structure objects, so -- like ``AtomArray`` in
:mod:`biorazer.structure.objects.bt_atom_array` -- they are enumerated here and
imported from here by the rest of the package.  biotite's bond *functions*
(``connect_via_distances``, ``connect_via_residue_names``) are helpers, not
objects; they stay imported from ``biotite.structure`` where they are used.
"""

from biotite.structure import BondList, BondType

__all__ = [
    "BondList",
    "BondType",
]
