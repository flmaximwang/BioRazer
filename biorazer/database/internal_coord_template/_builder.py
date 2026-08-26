# -*- coding: utf-8 -*-
"""Build per-residue, per-secondary-structure, per-rotamer ideal
``InternalCoord`` templates from the ``biorazer.database`` torsion and geometry
tables.  This is the **write / generative** path (the read path is
``InternalCoord.from_atomarray``; its side-chain topology
``SIDE_CHAIN_IC_PATH`` is reused here as the residue atom order).

Design (user-defined, 2026)
---------------------------
Each amino-acid template is an :class:`~biorazer.structure.objects.InternalCoord`
whose heavy atoms are ``N, CA, C, O, CB, <side chain>`` (Gly: no ``CB``/side
chain), anchored at ``{N, CA, C}`` --- the residue's three backbone atoms.  The
carbonyl ``O`` and the whole side chain (grown along ``SIDE_CHAIN_IC_PATH``)
are reconstructible from the anchor.  A template is a **per-conformer ideal
snapshot**:

* backbone: canonical Engh-Huber / Rosetta ideal; the secondary-structure
  ``phi``/``psi``/``omega`` means (from ``torsion_angle.backbone``) are carried
  on the template instance (``phi``/``psi``/``omega`` attrs) --- a single
  residue's own atoms cannot encode them (they need the neighbor residues).
* side chain: built to a chosen **rotamer**; the rotatable chi torsions are set
  exactly (subtree rigid rotation), the rest ideal.  Because the
  ``SIDE_CHAIN_IC_PATH`` quads reference backbone atoms, the side chain
  ``dihedra``/``bond_angles`` stored in the template are measured *at that
  conformer*, so ``to_coords()`` reproduces it exactly (round-trip invariant;
  this is the Rosetta residue-params model).
"""

from __future__ import annotations

import copy

import numpy as np

from biorazer.structure.objects.internal_coords import InternalCoord, AtomRecord, _place, _dihedral
from biorazer.database.internal_coord_template._topology import SIDE_CHAIN_IC_PATH
from biorazer.database.bond.backbone import AMINO_ACID_BOND_LENGTH, AMINO_ACID_BOND_ANGLE
from biorazer.database.bond.sidechain import AMINO_ACID_SIDECHAIN_BOND
from biorazer.database.torsion_angle.backbone import SS_BB_TORSION_ANGLE
from biorazer.database.torsion_angle.sidechain import (
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
            targets[rot] = _ROT_BIN[bins[i]]
    return targets


# --------------------------------------------------------------------------- #
# canonical ideal-residue coordinates (built from the database tables, so the
# runtime has no Rosetta / external dependency).
# --------------------------------------------------------------------------- #
def _build_coords(resn):
    """Place N, CA, C (canonical backbone), O, then every side-chain atom along
    ``SIDE_CHAIN_IC_PATH`` using the database bond/angle/dihedral values.
    Returns ``{atom_name: coord}`` for all heavy atoms."""
    N = np.array([0.0, 0.0, 0.0])
    CA = N + AMINO_ACID_BOND_LENGTH[("N", "CA")]["mean"] * np.array([1.0, 0.0, 0.0])
    C = (CA + AMINO_ACID_BOND_LENGTH[("CA", "C")]["mean"]
         * np.array([np.cos(np.radians(180 - AMINO_ACID_BOND_ANGLE[("N","CA","C")]["mean"])),
                     np.sin(np.radians(180 - AMINO_ACID_BOND_ANGLE[("N","CA","C")]["mean"])), 0.0]))
    coord = {"N": N, "CA": CA, "C": C}
    blen = AMINO_ACID_BOND_LENGTH[("C", "O")]["mean"]
    bang = AMINO_ACID_BOND_ANGLE[("CA", "C", "O")]["mean"]
    coord["O"] = _place(N, CA, C, blen, bang, 180.0)
    # side chain, in SIDE_CHAIN_IC_PATH topological order
    for quad in SIDE_CHAIN_IC_PATH[resn]:
        i, j, k, l = quad
        if not all(n in coord for n in (i, j, k)):
            continue
        entry = AMINO_ACID_SIDECHAIN_BOND[resn][quad]
        blen = entry["mean"]
        bang = entry["angle"]
        dih = SIDECHAIN_IC_DIHEDRAL[resn][quad]
        coord[l] = _place(coord[i], coord[j], coord[k], blen, bang, dih)
    return coord


# --------------------------------------------------------------------------- #
# chi-set machinery (exact rotamer application via subtree rigid rotation)
# --------------------------------------------------------------------------- #
def _parent_map(resn):
    pm = {l: k for _, _, k, l in SIDE_CHAIN_IC_PATH[resn]}
    pm.setdefault("CA", "N")
    pm.setdefault("C", "CA")
    pm.setdefault("O", "CA")
    pm["N"] = None
    return pm


def _subtree(resn, root):
    pm = _parent_map(resn)
    res = {root}
    changed = True
    while changed:
        changed = False
        for atom, par in pm.items():
            if par in res and atom not in res:
                res.add(atom)
                changed = True
    return res


def _rotate_about(coords, axis_pt, axis_dir, deg):
    axis_dir = axis_dir / np.linalg.norm(axis_dir)
    rad = np.radians(deg)
    c, s = np.cos(rad), np.sin(rad)

    def rot(v):
        return (c * v + s * np.cross(axis_dir, v)
                + (1 - c) * np.dot(v, axis_dir) * axis_dir)
    return {k: axis_pt + rot(p - axis_pt) for k, p in coords.items()}


def _set_chi(resn, coord, chi_quad, target):
    i, j, k, l = chi_quad
    axis_dir = coord[k] - coord[j]
    axis_pt = coord[k]
    group = {a: c for a, c in coord.items() if a in _subtree(resn, k)}
    others = {a: c for a, c in coord.items() if a not in _subtree(resn, k)}
    current = _dihedral(coord[i], coord[j], coord[k], coord[l])
    return {**others, **_rotate_about(group, axis_pt, axis_dir, target - current)}


def _apply_rotamer(resn, coord, targets):
    for quad in SIDECHAIN_CHI[resn]:
        if quad[3] in targets:
            coord = _set_chi(resn, coord, quad, targets[quad[3]])
    return coord


# --------------------------------------------------------------------------- #
# measurement -> InternalCoord
# --------------------------------------------------------------------------- #
def _element(atom_name):
    """Heavy-atom element from a PDB atom name (first letter for N/O/S, else C)."""
    return atom_name[0] if atom_name[0] in ("N", "O", "S") else "C"


def _measure_ic(resn, coord, ss, rotamer, phi, psi, omega):
    order = ["N", "CA", "C", "O"]
    seen = set(order)
    for quad in SIDE_CHAIN_IC_PATH[resn]:
        l = quad[3]
        if l in coord and l not in seen:
            order.append(l)
            seen.add(l)
    atoms = [AtomRecord(res_name=resn, name=nm, element=_element(nm),
                        chain_id="A", res_id=1) for nm in order]
    ic = InternalCoord(atoms=atoms)
    idx = {nm: n for n, nm in enumerate(order)}
    ic.anchor = {idx["N"]: tuple(np.asarray(coord["N"], float)),
                 idx["CA"]: tuple(np.asarray(coord["CA"], float)),
                 idx["C"]: tuple(np.asarray(coord["C"], float))}

    def record(quad):
        i, j, k, l = quad
        ci, cj, ck, cl = (coord[nm] for nm in quad)
        blen = float(np.linalg.norm(cl - ck))
        cos = np.dot(cj - ck, cl - ck) / (np.linalg.norm(cj - ck) * np.linalg.norm(cl - ck))
        bang = float(np.degrees(np.arccos(np.clip(cos, -1, 1))))
        ic.bond_distances[(idx[k], idx[l])] = blen
        ic.bond_angles[(idx[j], idx[k], idx[l])] = bang
        ic.dihedra[(idx[i], idx[j], idx[k], idx[l])] = _dihedral(ci, cj, ck, cl)

    # backbone O branch
    record(("N", "CA", "C", "O"))
    # side chain growth
    for quad in SIDE_CHAIN_IC_PATH[resn]:
        if all(n in coord for n in quad):
            record(quad)
    ic.ss = ss
    ic.phi = float(phi)
    ic.psi = float(psi)
    ic.omega = float(omega)
    ic.rotamer = rotamer
    return ic


# --------------------------------------------------------------------------- #
# public API
# --------------------------------------------------------------------------- #
def ss_torsions(ss):
    """``{phi, psi, omega}`` mean values for a secondary-structure class, from
    :data:`torsion_angle.backbone.SS_BB_TORSION_ANGLE`."""
    v = SS_BB_TORSION_ANGLE[ss]
    return {"phi": v["phi"]["mean"], "psi": v["psi"]["mean"], "omega": v["omega"]["mean"]}


def build_template(resn, ss, rotamer="canonical"):
    """Build one ``InternalCoord`` template for ``resn`` at ``ss`` and ``rotamer``.

    Parameters
    ----------
    resn : str
        Three-letter residue name (upper-case).
    ss : str
        Secondary-structure class key of ``SS_BB_TORSION_ANGLE``.
    rotamer : str
        A name from :func:`rotamer_names` (default ``"canonical"``).
    """
    coord = _build_coords(resn)
    coord = _apply_rotamer(resn, coord, rotamer_targets(resn, rotamer))
    t = ss_torsions(ss)
    return _measure_ic(resn, coord, ss, rotamer, t["phi"], t["psi"], t["omega"])


def build_template_direct(resn, ss, rotamer="canonical"):
    """Direct-fill an ``InternalCoord`` from the database tables -- no all-atom
    Cartesian placement (the user's "just fill the parameters" model).

    Only the anchor frame ``N/CA/C`` is placed -- the anchor is required by
    :class:`~biorazer.structure.objects.InternalCoord` to locate the residue --
    and the carbonyl ``O`` plus every side-chain atom are described purely by
    their ``(bond, bond-angle, dihedral)`` entries taken straight from
    ``bond.backbone`` / ``bond.sidechain`` / ``torsion_angle.sidechain``.

    .. warning::
        This path is **exact only for the canonical conformer**.  A rotamer is
        a *rigid subtree rotation* about a chi bond, and the chi torsion
        (``SIDECHAIN_CHI``, e.g. ``(N,CA,CB,OG)``) is **not a key of the stored
        grow-dihedrals** (``SIDECHAIN_IC_DIHEDRAL`` keys them by the
        ``SIDE_CHAIN_IC_PATH`` grow quad, e.g. ``(CA,N,CB,OG)``).  Setting chi
        therefore shifts the stored grow-dihedral geometrically -- a value that
        is not in any table -- so a pure table fill cannot reproduce a rotamer
        conformer.  Hence this function only implements ``rotamer="canonical"``.
    """
    if rotamer != "canonical":
        raise NotImplementedError(
            f"build_template_direct supports only rotamer='canonical' "
            f"(a rotamer is a rigid rotation; its stored grow-dihedrals must be "
            f"measured from geometry -- use build_template).")

    bl = AMINO_ACID_BOND_LENGTH
    ba = AMINO_ACID_BOND_ANGLE

    order = ["N", "CA", "C", "O"]
    seen = set(order)
    for quad in SIDE_CHAIN_IC_PATH[resn]:
        if quad[3] not in seen:
            order.append(quad[3])
            seen.add(quad[3])
    atoms = [AtomRecord(res_name=resn, name=nm, element=_element(nm),
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

    # carbonyl O branch (N, CA, C, O): trans
    ic.bond_distances[(idx["C"], idx["O"])] = bl[("C", "O")]["mean"]
    ic.bond_angles[(idx["CA"], idx["C"], idx["O"])] = ba[("CA", "C", "O")]["mean"]
    ic.dihedra[(idx["N"], idx["CA"], idx["C"], idx["O"])] = 180.0

    # side chain straight from the tables
    sc_bond = AMINO_ACID_SIDECHAIN_BOND[resn]
    for quad in SIDE_CHAIN_IC_PATH[resn]:
        entry = sc_bond[quad]
        i, j, k, l = (idx[nm] for nm in quad)
        ic.bond_distances[(k, l)] = entry["mean"]
        ic.bond_angles[(j, k, l)] = entry["angle"]
        ic.dihedra[(i, j, k, l)] = SIDECHAIN_IC_DIHEDRAL[resn][quad]

    t = ss_torsions(ss)
    ic.ss = ss
    ic.phi = float(t["phi"])
    ic.psi = float(t["psi"])
    ic.omega = float(t["omega"])
    ic.rotamer = "canonical"
    return ic


def make_residue_templates(resn):
    """``{ss: {rotamer: InternalCoord}}`` for all SS classes x common rotamers."""
    out = {}
    for ss in SS_BB_TORSION_ANGLE:
        out[ss] = {r: build_template(resn, ss, r) for r in rotamer_names(resn)}
    return out