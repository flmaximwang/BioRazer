"""biotite's ``AtomArray`` / ``AtomArrayStack``, re-exported for the object layer.

``AtomArray`` is biorazer's core structure object: one flat array of atoms (a
record per row, annotations as parallel arrays) plus optional bonds.  The
convention (user-established) is that the third-party structure objects are
imported once, here inside :mod:`biorazer.structure.objects`, and composed by
the rest of the package -- every other module imports them from
``biorazer.structure.objects`` instead of importing biotite directly.  That
keeps the surface of the third-party dependency in one place and makes it
visible in one file which object is in use.

``AtomArrayStack`` (the multi-model variant) is re-exported alongside it, but
nothing in biorazer uses it yet.
"""

from biotite.structure import AtomArray, AtomArrayStack

__all__ = [
    "AtomArray",
    "AtomArrayStack",
]
