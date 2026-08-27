# -*- coding: utf-8 -*-
"""Ideal per-residue ``InternalCoord`` templates.

The whole public surface lives in :mod:`.template` (``build_template`` and
``get_available_specs``).  Prefer importing from the package path::

    from biorazer.database.molecule.icoor.protein import template

    template.get_available_specs("ALA")
    ic = template.build_template("ALA", "alpha-helix", "canonical")
"""

from biorazer.database.molecule.icoor.protein.template.template import (
    build_template,
    get_available_specs,
)

__all__ = ["build_template", "get_available_specs"]
