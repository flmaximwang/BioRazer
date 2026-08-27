# -*- coding: utf-8 -*-
"""Per-amino-acid ideal ``InternalCoord`` templates, one module per residue.

This package is the **generative / template** store of the internal-coordinate
layer: for every standard amino acid it provides ideal residue geometries as
:class:`~biorazer.structure.objects.InternalCoord` objects, keyed by **two
axes**:

* **secondary structure** (the ``phi``/``psi``/``omega`` means from
  :data:`biorazer.database.molecule.bond.dihedral.protein.SS_BB_TORSION_ANGLE`),
  carried on each template;
* **side-chain rotamer** (chi is a first-class axis; the chi1/chi2 bin-center
  ideals ``g-/t/g+`` at ``-60/180/+60`` from
  :data:`biorazer.database.molecule.bond.dihedral.protein.ROTAMER_BIN`),
  with the side chain rebuilt to that rotamer and its IC-frame
  bond/angle/dihedral taken directly from the database tables (no
  measure-from-coordinates round-trip; see :mod:`._builder`).

Each residue lives in its own module (``gly.py``, ``ser.py``, ...) exposing
``RESN`` and its **own named template constants** ``IC_<Res>_<SS>[_<rotamer>]``,
each built directly from :func:`build_template` (no intermediate
``{ss: {rotamer}}`` dict), so the named templates live right next to their
residue.

Usage::

    from biorazer.database.molecule.icoor.protein.template import get_template

    ic = get_template("SER", "alpha-helix", "g-")   # one template
    ic.phi, ic.psi, ic.omega                          # -60.0, -45.0, 180.0
    ic.to_atomarray()                                 # build the residue

Named single-template constants (real, visible, per residue)::

    from biorazer.database.molecule.icoor.protein.template import IC_Gly_HelixAlpha, IC_Ser_HelixAlpha_gminus

    IC_Gly_HelixAlpha            # == get_template("GLY", "alpha-helix")
    IC_Ser_HelixAlpha_gminus     # == get_template("SER", "alpha-helix", "g-")

Each is a literal constant ``IC_<Res>_<SS>[_<rotamer>]`` defined in its residue
module (``IC_Gly_*`` in ``gly.py``, ``IC_Ser_*`` in ``ser.py``, ...), built
directly from :func:`build_template` -- ``<Res>`` the residue title token
(``Gly``), ``<SS>`` a secondary-structure token (``HelixAlpha`` for
``alpha-helix``, see :data:`SS_TOKEN`), plus an optional sanitized rotamer
suffix (``gminus``/``trans``/``gplus``, 2-axis rotamers joined by ``_``).  The
canonical representative is ``IC_<Res>_<SS>``.  They are re-exported here so the
package path works.

A residue's own atoms cannot determine ``phi/psi`` (those need its neighbours),
so the per-SS differentiation is carried as the ``ic.phi/psi/omega`` attrs with
the geometry anchored at ``{N, CA, C}`` --- the Rosetta residue-params model.
See :mod:`._builder` for the exact build rule and the per-conformer caveat.

Current build: 20 residues, ``get_template`` over 12 SS classes x
(1 canonical + 0/3/9 common rotamers).  The named constants in every residue
module call :func:`build_template` eagerly, so importing this package builds all
templates up front.
"""

from __future__ import annotations

import importlib

from biorazer.database.molecule.bond.dihedral.protein import SS_BB_TORSION_ANGLE

from ._builder import (  # noqa: F401  (re-exported helpers)
    build_template,
    rotamer_names,
    ss_torsions,
)
from biorazer.database.molecule.icoor.protein.topology import (  # noqa: F401
    IC_PATH,
    MAINCHAIN_ATOMS,
)
from ._naming import SS_TOKEN, ROTAMER_TOKEN, rotamer_token, template_token  # noqa: F401

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
            f"biorazer.database.molecule.icoor.protein.template.{resn.lower()}")
    return _MODULES[resn]


def get_template(resn, ss, rotamer="canonical"):
    """A single ``InternalCoord`` template for ``resn`` at ``ss`` x ``rotamer``.

    Resolves to the residue module's literal named constant
    ``IC_<Res>_<SS>[_<rotamer>]``, which is built directly from
    :func:`build_template`.
    """
    name = template_token(resn, ss, rotamer)
    try:
        return getattr(_module(resn), name)
    except AttributeError:
        raise KeyError(
            f"no {resn} template for ss={ss!r} rotamer={rotamer!r} "
            f"(ss in {SS_CLASSES}; rotamers of {resn}: {rotamer_names(resn)})"
        ) from None


def templates(resn):
    """The ``{ss: {rotamer: InternalCoord}}`` map for a residue (built on the
    fly from its named :func:`build_template` constants)."""
    return {
        ss: {rot: get_template(resn, ss, rot) for rot in rotamer_names(resn)}
        for ss in SS_CLASSES
    }


__all__ = [
    "RESIDUES", "SS_CLASSES",
    "IC_PATH", "MAINCHAIN_ATOMS",
    "build_template",
    "rotamer_names", "ss_torsions",
    "templates", "get_template",
    "SS_TOKEN", "ROTAMER_TOKEN", "rotamer_token", "template_token",
]

# Promote every residue module's literal `IC_*` named templates into the package
# namespace, so `from ...template import IC_Gly_HelixAlpha` works.  Each residue
# module (gly.py, ser.py, ...) carries its own `IC_<Res>_<SS>[_<rotamer>]`
# constants and an `__all__` (generated by
# scripts/generate_internal_coord_template_named.py); here we only bind them.
for _resn in RESIDUES:
    _mod = importlib.import_module(f"{__name__}.{_resn.lower()}")
    _names = list(getattr(_mod, "__all__", ()))
    for _n in _names:
        globals()[_n] = getattr(_mod, _n)
    __all__ = list(__all__) + _names
