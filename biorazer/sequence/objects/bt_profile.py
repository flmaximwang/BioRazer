"""biotite's ``SequenceProfile``, re-exported for the object layer.

A ``SequenceProfile`` is the per-column counts of a set of equal-length
sequences: a ``symbols`` matrix (columns x alphabet) of how often each letter
occurs at each column, a ``gaps`` vector, and the ``alphabet`` the columns are
indexed by.  It is the object every conservation / entropy / sequence-logo
routine in :mod:`biorazer.sequence.analysis.alignment` takes, either directly
or via ``SequenceProfile.from_alignment``.
"""

from biotite.sequence.profile import SequenceProfile

__all__ = [
    "SequenceProfile",
]
