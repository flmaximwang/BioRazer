"""biotite's alignment objects, re-exported for the object layer.

``Alignment`` is the multiple-sequence-alignment object: the ungapped
sequences plus a ``trace`` matrix that says, per sequence, which columns are
letters and which are gaps.  ``SubstitutionMatrix`` is the scoring matrix an
alignment is computed with (``SubstitutionMatrix.std_protein_matrix()`` for the
BLOSUM62 default) -- an object rather than a function, so it is enumerated
here.

biotite's alignment *functions* (``align_multiple``, ``align_optimal``,
``KmerTable``-based seeds) are helpers, not objects; they stay imported from
``biotite.sequence.align`` where they are used.
"""

from biotite.sequence.align import Alignment, SubstitutionMatrix

__all__ = [
    "Alignment",
    "SubstitutionMatrix",
]
