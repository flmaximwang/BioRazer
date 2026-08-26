# -*- coding: utf-8 -*-
"""Per-atom physical property tables.

* :mod:`.radius` -- van der Waals radii (Å) by element.
* :mod:`.charge` -- formal / partial atomic charges (placeholder; no data
  yet).

Every numeric entry is ``{mean, std, lb, up, source}`` (``np.nan`` when
unknown).
"""
