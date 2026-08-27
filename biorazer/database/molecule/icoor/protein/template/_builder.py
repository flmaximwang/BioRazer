# -*- coding: utf-8 -*-
"""Build per-residue, per-secondary-structure, per-rotamer ideal
``InternalCoord`` templates from the ``biorazer.database`` torsion and geometry
tables.  This is the **write / generative** path (the read path is
``InternalCoord.from_atomarray``; its side-chain topology
``IC_PATH`` is reused here as the residue atom order).

Design (user-defined, 2026)
---------------------------
Each amino-acid template is an :class:`~biorazer.structure.objects.InternalCoord`
whose heavy atoms are ``N, CA, C, O, CB, <side chain>`` (Gly: no ``CB``/side
chain), anchored at ``{N, CA, C}`` --- the residue's three backbone atoms.  The
carbonyl ``O`` and the whole side chain (grown along ``IC_PATH``)
are reconstructible from the anchor.  A template is a **per-conformer ideal
snapshot**:

* backbone: canonical Engh-Huber / Rosetta ideal; the secondary-structure
  ``phi``/``psi``/``omega`` means (from ``torsion_angle.backbone``) are carried
  on the template instance (``phi``/``psi``/``omega`` attrs) --- a single
  residue's own atoms cannot encode them (they need the neighbor residues).
* side chain: built to a chosen **rotamer**; the rotatable chi torsions are set
  exactly (subtree rigid rotation), the rest ideal.

**Direct table fill (no Cartesian round-trip).**  A template is filled straight
from the database tables: the anchor frame ``{N, CA, C}`` is the only piece of
absolute geometry (placed identically to the old coordinate build), and every
bond length / bond angle / dihedral is assigned directly from
``bond.length`` / ``bond.angle`` / ``bond.dihedral``.  Because a rotamer is a
*rigid* subtree rotation about a chi bond:

* bond lengths and bond angles are invariant under the rotation -> taken
  straight from the tables;
* every non-chi grow dihedral is also invariant (the dihedral either lies
  inside the rigid group or its reference frame sits outside the rotated
  subtree) -> taken straight from the tables;
* the chi dihedrals themselves equal the rotamer's bin targets ->
  overridden from ``rotamer_targets()``.

This is verified exhaustively: for all 20 residues x 12 SS classes x common
rotamers, the table-filled template reproduces the old
measure-from-coordinates geometry bit-for-bit (bond/angle/dihedral all equal).
"""

from __future__ import annotations

import numpy as np

from biorazer.structure.objects.internal_coords import InternalCoord, InternalCoordAtom
from biorazer.database.molecule.icoor.protein.topology import BACKBONE_IC_PATH, IC_PATH
from biorazer.database.molecule.bond.length.protein import AMINO_ACID_SIDECHAIN_BOND, AMINO_ACID_BOND_LENGTH
from biorazer.database.molecule.bond.angle.generic import AMINO_ACID_BOND_ANGLE
from biorazer.database.molecule.bond.angle.protein import AMINO_ACID_SIDECHAIN_BOND_ANGLE
from biorazer.database.molecule.bond.dihedral.protein import (
    SS_BB_TORSION_ANGLE,
    SIDECHAIN_CHI,
    SIDECHAIN_IC_DIHEDRAL,
    ROTAMER_BIN,
)


# --------------------------------------------------------------------------- #
# rotamer scheme
# --------------------------------------------------------------------------- #
# chi bin centers live in torsion_angle.sidechain (single source of truth).
# g-/t/g+ = -60/180/+60 (standard Dunbrack rotamer definitions).
_ROT_BIN = ROTAMER_BIN

#: Per amino acid: number of freely rotatable chi used to define its common
#: rotamers.  ``0`` -> single canonical template; ``1`` -> 3 chi1 rotamers;
#: ``2`` -> 9 chi1 x chi2 rotamers (chi3+ kept canonical).  Pro keeps a single
#: canonical ring-puckering template (its chi are ring-constrained).
_COMMON_CHI_AXES = {
    "ALA": 0, "GLY": 0, "PRO": 0,
    "CYS": 1, "SER": 1, "THR": 1, "VAL": 1,
    "ASN": 2, "ASP": 2, "GLN": 2, "GLU": 2, "HIS": 2, "ILE": 2, "LEU": 2,
    "LYS": 2, "MET": 2, "PHE": 2, "TRP": 2, "TYR": 2, "ARG": 2,
}


def rotamer_names(resn):
    """List of common rotamer names for a residue, each with a ``"canonical"``
    representative first (chi is a first-class axis where the residue has
    rotatable chi)."""
    n = _COMMON_CHI_AXES.get(resn, 0)
    if n == 0:
        return ["canonical"]
    if n == 1:
        return ["canonical"] + list(_ROT_BIN)
    return ["canonical"] + [f"{a}/{b}" for a in _ROT_BIN for b in _ROT_BIN]


def rotamer_targets(resn, rotamer):
    """Map a rotamer name -> {rotator atom: chi target (deg)}.  For the
    ``2``-axis residues the chi3/chi4 rotators are kept at their canonical
    (ideal) value (measured from the canonical build unless already set)."""
    chis = SIDECHAIN_CHI.get(resn, [])
    rotators = [q[3] for q in chis]
    if _COMMON_CHI_AXES.get(resn, 0) == 0:
        return {}
    if _COMMON_CHI_AXES.get(resn, 0) == 1:
        bins = [rotamer]
    else:
        bins = rotamer.split("/")
    targets = {}
    for i, rot in enumerate(rotators):
        if i < len(bins) and bins[i] in _ROT_BIN:
            targets[rot] = _ROT_BIN[bins[i]]["mean"]
    return targets


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _element(atom_name):
    """Heavy-atom element from a PDB atom name (first letter for N/O/S, else C)."""
    return atom_name[0] if atom_name[0] in ("N", "O", "S") else "C"


# --------------------------------------------------------------------------- #
# public API
# --------------------------------------------------------------------------- #
def ss_torsions(ss):
    """``{phi, psi, omega}`` mean values for a secondary-structure class, from
    :data:`biorazer.database.molecule.bond.dihedral.protein.by_ss.SS_BB_TORSION_ANGLE`."""
    v = SS_BB_TORSION_ANGLE[ss]
    return {"phi": v["phi"]["mean"], "psi": v["psi"]["mean"], "omega": v["omega"]["mean"]}


def build_template(resn, ss, rotamer="canonical"):
    """Build one ``InternalCoord`` template for ``resn`` at ``ss`` and ``rotamer``
    by **direct table fill** -- no all-atom Cartesian placement, no
    measure-from-coordinates round-trip.

    The anchor frame ``{N, CA, C}`` is the only absolute geometry (N at origin,
    CA along +x, C from the N-CA-C bond angle), placed identically to the old
    coordinate build.  Every bond length / bond angle / dihedral is assigned
    directly from the database tables; the chi dihedrals of a rotamer are
    overridden with the rotamer's bin targets (rigid subtree rotation: all other
    values are invariant, see the module docstring).

    Parameters
    ----------
    resn : str
        Three-letter residue name (upper-case).
    ss : str
        Secondary-structure class key of ``SS_BB_TORSION_ANGLE``.
    rotamer : str
        A name from :func:`rotamer_names` (default ``"canonical"``).
    """
    bl = AMINO_ACID_BOND_LENGTH
    ba = AMINO_ACID_BOND_ANGLE

    order = ["N", "CA", "C", "O"]
    seen = set(order)
    for quad in IC_PATH[resn]:
        l = quad[3]
        if l not in seen:
            order.append(l)
            seen.add(l)
    atoms = [InternalCoordAtom(res_name=resn, name=nm, element=_element(nm),
                        chain_id="A", res_id=1) for nm in order]
    ic = InternalCoord(atoms=atoms)
    idx = {nm: n for n, nm in enumerate(order)}

    # anchor frame: N at origin, CA along +x, C from CA-C length + N-CA-C angle
    ic.anchor[idx["N"]] = (0.0, 0.0, 0.0)
    ca = np.array([bl[("N", "CA")]["mean"], 0.0, 0.0], float)
    ic.anchor[idx["CA"]] = tuple(ca)
    ang = np.radians(180.0 - ba[("N", "CA", "C")]["mean"])
    c = ca + bl[("CA", "C")]["mean"] * np.array([np.cos(ang), np.sin(ang), 0.0])
    ic.anchor[idx["C"]] = tuple(c)

    # anchor triple is a fully-specified rigid body: its own bonds and angle
    ic.bond_distances[(idx["N"], idx["CA"])] = bl[("N", "CA")]["mean"]
    ic.bond_distances[(idx["CA"], idx["C"])] = bl[("CA", "C")]["mean"]
    ic.bond_angles[(idx["N"], idx["CA"], idx["C"])] = ba[("N", "CA", "C")]["mean"]

    # carbonyl O branch -- the "intra" backbone grow quads (only quads whose
    # (k, l) bond geometry exists in the generic tables are grown; template
    # residues have no OXT so only the O quad lands here)
    for i, j, k, l in BACKBONE_IC_PATH["intra"]:
        if (k, l) not in bl:
            continue
        ic.bond_distances[(idx[k], idx[l])] = bl[(k, l)]["mean"]
        ic.bond_angles[(idx[j], idx[k], idx[l])] = ba[(j, k, l)]["mean"]
        ic.dihedra[(idx[i], idx[j], idx[k], idx[l])] = 180.0

    # side chain straight from the tables; chi quads overridden by the rotamer
    sc_bond = AMINO_ACID_SIDECHAIN_BOND[resn]
    sc_angle = AMINO_ACID_SIDECHAIN_BOND_ANGLE[resn]
    targets = rotamer_targets(resn, rotamer)
    for quad in IC_PATH[resn]:
        i, j, k, l = quad
        ic.bond_distances[(idx[k], idx[l])] = sc_bond[(k, l)]["mean"]
        ic.bond_angles[(idx[j], idx[k], idx[l])] = sc_angle[(j, k, l)]["mean"]
        ic.dihedra[(idx[i], idx[j], idx[k], idx[l])] = (
            targets[l] if l in targets else SIDECHAIN_IC_DIHEDRAL[resn][quad]["mean"]
        )

    t = ss_torsions(ss)
    ic.ss = ss
    ic.phi = float(t["phi"])
    ic.psi = float(t["psi"])
    ic.omega = float(t["omega"])
    ic.rotamer = rotamer
    return ic


def build_template_direct(resn, ss, rotamer="canonical"):
    """Direct-fill alias of :func:`build_template`.

    Kept for backward compatibility: ``build_template`` is now the direct
    table-fill path (no Cartesian round-trip), so both names are equivalent.
    """
    return build_template(resn, ss, rotamer)


def make_residue_templates(resn):
    """``{ss: {rotamer: InternalCoord}}`` for all SS classes x common rotamers."""
    out = {}
    for ss in SS_BB_TORSION_ANGLE:
        out[ss] = {r: build_template(resn, ss, r) for r in rotamer_names(resn)}
    return out
