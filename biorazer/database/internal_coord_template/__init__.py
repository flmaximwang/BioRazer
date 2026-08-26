# -*- coding: utf-8 -*-
"""Per-amino-acid ideal ``InternalCoord`` templates, one module per residue.

This package is the **generative / template** store of the internal-coordinate
layer: for every standard amino acid it provides ideal residue geometries as
:class:`~biorazer.structure.objects.InternalCoord` objects, keyed by **two
axes**:

* **secondary structure** (the ``phi``/``psi``/``omega`` means from
  ``torsion_angle.backbone.SS_BB_TORSION_ANGLE``), carried on each template;
* **side-chain rotamer** (chi is a first-class axis; the chi1/chi2 bin-center
  ideals ``g-/t/g+`` at ``-60/180/+60`` from ``torsion_angle.sidechain``),
  with the side chain rebuilt to that rotamer and its IC-frame
  bond/angle/dihedral measured at that conformer.

Each residue lives in its own module (``ala.py``, ``arg.py``, ...) exposing
``RESN`` and ``<RESN>_TEMPLATES = {ss: {rotamer: InternalCoord}}``.  Modules
are imported lazily on access so ``import biorazer.database.internal_coord_template``
stays cheap.

Usage::

    from biorazer.database.internal_coord_template import get_template, templates

    ic = get_template("SER", "alpha-helix", "g-")   # one template
    ser = templates("SER")                           # {ss: {rotamer: ic}}
    ic.phi, ic.psi, ic.omega                          # -60.0, -45.0, 180.0
    ic.to_atomarray()                                 # build the residue

A residue's own atoms cannot determine ``phi/psi`` (those need its neighbours),
so the per-SS differentiation is carried as the ``ic.phi/psi/omega`` attrs with
the geometry anchored at ``{N, CA, C}`` --- the Rosetta residue-params model.
See :mod:`._builder` for the exact build rule and the per-conformer caveat.

Current build: 20 residues, ``templates(resn)`` = 12 SS classes x (1/3/9)
common rotamers.
"""

from __future__ import annotations

import importlib

from biorazer.database.torsion_angle.backbone import SS_BB_TORSION_ANGLE

from ._builder import (  # noqa: F401  (re-exported helpers)
    build_template,
    make_residue_templates,
    rotamer_names,
    ss_torsions,
)
from ._topology import MAINCHAIN_ATOMS, SIDE_CHAIN_IC_PATH  # noqa: F401

#: Standard amino-acid residues (three-letter, upper-case), module-per-residue.
RESIDUES = (
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
)

#: Secondary-structure classes (keys of ``SS_BB_TORSION_ANGLE``).
SS_CLASSES = tuple(sorted(SS_BB_TORSION_ANGLE))

_MODULES = {}


def _module(resn):
    """Import (once) the residue module ``<resn>.py`` and return it."""
    resn = resn.upper()
    if resn not in RESIDUES:
        raise ValueError(f"unknown residue {resn!r}; known: {list(RESIDUES)}")
    if resn not in _MODULES:
        _MODULES[resn] = importlib.import_module(
            f"biorazer.database.internal_coord_template.{resn.lower()}")
    return _MODULES[resn]


def templates(resn):
    """The ``{ss: {rotamer: InternalCoord}}`` map for a residue (lazy)."""
    return getattr(_module(resn), f"{resn}_TEMPLATES")


def get_template(resn, ss, rotamer="canonical"):
    """A single ``InternalCoord`` template for ``resn`` at ``ss`` x ``rotamer``."""
    try:
        return templates(resn)[ss][rotamer]
    except KeyError:
        raise KeyError(f"no {resn} template for ss={ss!r} rotamer={rotamer!r} "
                       f"(ss in {SS_CLASSES}; rotamers of {resn}: {rotamer_names(resn)})") from None


__all__ = [
    "RESIDUES", "SS_CLASSES",
    "MAINCHAIN_ATOMS", "SIDE_CHAIN_IC_PATH",
    "build_template", "make_residue_templates", "rotamer_names", "ss_torsions",
    "templates", "get_template", "_module",
]