"""Converters from an ``AtomArraySelector`` to a per-atom selection on an ``AtomArray``.

``AtomArraySelector_AtomArrayMask`` and ``AtomArraySelector_AtomArrayIndices`` are
biorazer's ``A_B`` in-memory bridge converters -- ``A`` is the source object, ``B``
the target -- so they implement the parent's
:meth:`~biorazer.io.Converter.convert` (one transform, no ``read()`` / ``write()``):
the selector is held in ``input_io`` and ``convert()`` returns the selection.

Both need the array the selection is *made against*: a selector is a set of patterns
over ``ins_code`` / ``chain`` / ``resi`` / ``name``, and only an ``AtomArray`` says
which atoms exist.  The array is therefore a parameter of ``convert()`` -- the shape
the ``quads`` / ``anchor`` parameters of
:class:`~biorazer.structure.bridge.atom_array.AtomArray_InternalCoord` already use:
the source object comes from ``input_io``, the per-call data is an argument.

The two targets differ in ordering, which is why both exist:

* a **mask** is aligned with the array, so the rule order is invisible in it, and
* **indices** are flat (rule order first, then atom order within a rule), an order
  a mask cannot express.  ``dedupe=False`` keeps atoms matched by more than one
  rule once per match.
"""

from __future__ import annotations

import numpy as np

from biorazer.io import Converter
from biorazer.structure.objects import AtomArray, AtomArraySelector


class AtomArraySelector_AtomArrayMask(Converter):
    """Builds a boolean mask from an :class:`AtomArraySelector`.

    The mask is aligned with the input array (``True`` at the atoms any rule
    matches), so it can be used directly as ``atom_array[mask]`` and composed with
    the helpers in :mod:`biorazer.structure.selection.mask`.
    """

    def convert(self, atom_array: AtomArray) -> np.ndarray:
        """Select on ``atom_array`` with ``self.input_io``.

        Parameters
        ----------
        atom_array : AtomArray
            The array the patterns are matched against, and the array the
            returned mask is aligned with.

        Returns
        -------
        numpy.ndarray
            1D boolean mask of shape ``atom_array.shape``, ``True`` at every atom
            matched by at least one rule.
        """
        return self.input_io.mask(atom_array)


class AtomArraySelector_AtomArrayIndices(Converter):
    """Builds flat atom indices from an :class:`AtomArraySelector`.

    The order is the one a mask cannot express: rule by rule, and within a rule
    the atom order of the array.  It is the order
    :meth:`~biorazer.structure.objects.AtomArraySelector.to_csv` writes into a
    selection table, so indices and the exported table agree row by row.
    """

    def convert(self, atom_array: AtomArray, dedupe: bool = True) -> np.ndarray:
        """Select on ``atom_array`` with ``self.input_io``.

        Parameters
        ----------
        atom_array : AtomArray
            The array the patterns are matched against.
        dedupe : bool
            Keep an atom matched by more than one rule only at its first match
            (the default), or once per matching rule.

        Returns
        -------
        numpy.ndarray
            1D ``int`` array of atom indices.
        """
        return self.input_io.indices(atom_array, dedupe=dedupe)
