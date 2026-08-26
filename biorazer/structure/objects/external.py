"""External object re-exports for :mod:`biorazer.structure.objects`.

This package composes the "object" layer of ``biorazer.structure``.  The
convention (user-established) is that third-party structure objects are
imported here, in :mod:`.external`, and composed/used by the modules in this
package -- other modules inside :mod:`biorazer.structure` should import those
objects from here (or from the individual ``objects`` submodules) instead of
importing the third-party library (biotite) directly.

This keeps the surface of third-party dependencies on one place and makes the
packages below reusable without re-binding biotite symbols.
"""

from biotite.structure import AtomArray, AtomArrayStack

__all__ = [
    "AtomArray",
    "AtomArrayStack",
]