"""Converters between biotite's ``AtomArray`` and biorazer's ``InternalCoord``.

:class:`AtomArray_InternalCoord` (read path) and :class:`InternalCoord_AtomArray`
(write path) are biorazer's ``A_B`` in-memory bridge converters -- ``A`` is the
source object, ``B`` the target -- so they implement the parent's
:meth:`~biorazer.io.Converter.convert` (one transform, no ``read()`` /
``write()``): the source object is held in ``input_io`` and ``convert()`` returns
the result.

The conversion lives here, not on the objects: ``InternalCoord`` carries no
``from_atomarray`` / ``to_atomarray`` methods, so the objects stay free of
knowledge about each other (see
:mod:`biorazer.structure.objects`).  The two directions are **not** symmetric and
should not be expected to round-trip:

* the read path records every bond / angle / dihedral **measured from the input**
  array, so ``InternalCoord -> AtomArray -> InternalCoord`` on the same object
  reproduces the coordinates it was built from (``~1e-14`` A), and
* the write path *reconstructs* the coordinates from the internal coordinates by
  frame growth (:meth:`~biorazer.structure.objects.InternalCoord.to_coords`),
  then packs the atom records back into a new array.  Annotation categories that
  ``InternalCoordAtom`` does not carry are therefore not restored -- notably
  ``hetero``, which comes back ``False`` (biotite's ``AtomArray`` constructor
  default) as it does for the input array.
"""

from __future__ import annotations

import numpy as np

from biorazer.io import Converter
from biorazer.structure.objects import AtomArray, InternalCoord, InternalCoordAtom
from biorazer.structure.objects.internal_coords import NULL_ALT, dihedral


def _altloc(atom_array, index) -> str:
    """Read one ``AtomArray`` row's altLoc label, ``""`` when it has none.

    The ``altloc_id`` category only exists on an array read with
    ``altloc="all"`` (biotite 1.6 does not create it otherwise), and even then
    it spells "no alternate conformation" as a blank column, ``"."`` or
    ``"?"`` -- all normalised to ``""`` here, so a record has exactly one
    spelling for "none" (see :data:`~biorazer.structure.objects.internal_coords.NULL_ALT`).
    """
    labels = getattr(atom_array, "altloc_id", None)
    if labels is None:
        return ""
    label = str(labels[index])
    return "" if label in NULL_ALT else label


def _atom_record(atom_array, index) -> InternalCoordAtom:
    """Read one ``AtomArray`` row into an :class:`InternalCoordAtom`.

    Parameters
    ----------
    atom_array : AtomArray
        The source array.
    index : int
        Row (atom) index to read.

    Returns
    -------
    InternalCoordAtom
    """
    arr = atom_array
    return InternalCoordAtom(ins_code=str(arr.ins_code[index]),
                             chain_id=str(arr.chain_id[index]),
                             res_name=str(arr.res_name[index]),
                             res_id=int(arr.res_id[index]),
                             name=str(arr.atom_name[index]),
                             element=str(arr.element[index]),
                             altloc_id=_altloc(arr, index))


class AtomArray_InternalCoord(Converter):
    """Builds an :class:`InternalCoord` from a biotite ``AtomArray``.

    Protein-aware construction (**the default**, ``quads=None``) runs a
    two-pass build:

    * **Main-chain pass** (uniform, same for every residue): the backbone
      ``N -> CA -> C -> O`` is walked residue by residue, linking residue
      ``i``'s ``C`` to residue ``i+1``'s ``N`` (peptide bond).  It records
      the cross-residue quads ``(N_i, CA_i, C_i, N_{i+1})``,
      ``(CA_i, C_i, N_{i+1}, CA_{i+1})``, ``(C_i, N_{i+1}, CA_{i+1}, C_{i+1})``
      and the per-residue carbonyl branch ``(N, CA, C, O)`` -- the exact
      quads of :data:`~biorazer.database.molecule.icoor.protein.topology.BACKBONE_IC_PATH`
      (``"peptide"`` / ``"intra"`` groups).  **Every** quad, the carbonyl
      ``O`` branch included, stores the value **measured from the input**
      (``record`` reads it straight off ``arr``), which is what makes
      :meth:`~biorazer.structure.objects.InternalCoord.to_coords` reproduce the
      input to ``~1e-14`` A.  The constraint the
      measured ``O`` value satisfies in a real structure (sp2 coplanarity of
      the carbonyl carbon -> ``dihedral(N, CA, C, O) = psi - 180``, see
      :func:`~biorazer.database.molecule.icoor.protein.topology.carbonyl_o_dihedral`)
      is therefore *not* applied here -- the **write** paths (template /
      ``build_side_chain``) use it, because they have no ``O`` to read.
    * **Side-chain pass** (per residue): each standard amino acid's side
      chain is grown off the already-placed backbone using its per-residue
      grow-path table ``IC_PATH`` (chi rotamers; see
      ``biorazer.database.molecule.icoor.protein.topology``).  Non-standard /
      non-protein atoms (water, ligands, hydrogens) are not covered.

    Every atom of the input gets a record -- including the ``altloc_id``
    label, so a selector can filter on alternate conformations.  The **grow
    tree**, however, keys atoms by ``(chain, res_id, ins_code, name)``: e.g.
    per name the last copy in array order wins, and the peptide quads are
    only recorded when that picked ``C_i - N_{i+1}`` distance is within the
    C-N bond-length bound.  Picking a copy per atom is therefore a lottery on
    a multi-conformer structure, and losing it once orphans everything
    downstream -- measured on 2VB1, whose residue 5 ``N`` has an altloc A
    copy 1.375 A from ``C4`` (just over the 1.371 A bound) and a B copy at
    1.241 A: array order keeps B at 4->5 but B again at 43->44, where the A
    copy was the bonded one (1.313 A vs 1.454 A), so the chain breaks at 44
    and 654 of 2900 records become unplaceable
    (:class:`InternalCoord_AtomArray` then raises ``Unreachable atoms``).

    Reading with ``altloc="first"`` is **not** the fix: it keeps per-atom
    first copies, which on the same file break earlier still (961 records
    unplaceable, from residue 5).  An ``InternalCoord`` round trip needs a
    conformer set that is geometrically self-consistent, which is the
    caller's decision, not one this builder can make.

    Anchors default to the first three backbone atoms ``N, CA, C`` of every
    chain (one connected-component root per chain).

    Anchor-frame geometry is recorded so the anchor is a fully-specified
    rigid body: the ``N-CA`` and ``CA-C`` bonds of the anchor triple go into
    ``bond_distances`` and the ``N-CA-C`` bond angle into ``bond_angles``.
    (The peptide ``C_i - N_{i+1}`` bond is recorded by the cross-residue
    quads as usual; only the anchor triple itself has no dihedral, which is
    fine -- a dihedral needs four atoms, and the anchor is a rigid frame
    with no parent.)  This keeps ``anchor`` self-describing: a ``to_coords``
    round-trip on the anchor atoms alone needs no extra bookkeeping, and
    downstream code that modifies anchor positions can always recover the
    pair distances from the bond map.  These two records are set by
    :func:`record` for every quad that grows one of the anchor atoms (the
    carbonyl ``O`` branch and the peptide link), and any remaining missing
    bond of the anchor triple itself (``N-CA`` or ``CA-C`` of a terminal
    residue, or of a chain that never grows) is filled at the end of the
    per-chain loop.

    For a **general graph** (ligands, rings, arbitrary connectivity) pass
    explicit ``quads`` (a list of ``(i, j, k, l)`` atom-index quadruples) --
    the generic path; its anchors default to ``{0: first atom}``.
    """

    def convert(self, quads=None, anchor=None) -> InternalCoord:
        """Build the :class:`InternalCoord` from ``self.input_io``.

        Parameters
        ----------
        quads : list[tuple[int,int,int,int]] or None
            Explicit (i,j,k,l) atom-index quadruples for the general-graph
            path.  If ``None`` (default) the protein-aware two-pass build
            described in the class docstring is used.
        anchor : dict[int, tuple[float,float,float]] or None
            Absolute coordinates for anchor atoms given as a mapping.  If
            ``None`` (default) anchors are auto-detected (per chain's first
            ``N, CA, C`` for the protein path; ``{0: first atom}`` for the
            general ``quads`` path).

        Returns
        -------
        InternalCoord

        Notes
        -----
        Bond lengths and angles are derived from the input array for the
        parent/child pairs of each quad, so they are exact (not idealised).
        """
        from biorazer.database.molecule.icoor.protein.topology import (
            BACKBONE_IC_PATH,
            IC_PATH,
        )
        from biorazer.database.molecule.bond.length.protein import AMINO_ACID_BOND_LENGTH

        arr = self.input_io
        n = len(arr)
        atoms = [_atom_record(arr, i) for i in range(n)]
        ic = InternalCoord(atoms=atoms)

        def record(quad):
            """Fill bond/angle/dihedra for one quad from ``arr`` (exact)."""
            i, j, k, l = quad
            c0 = np.asarray(arr.coord[i], float)
            c1 = np.asarray(arr.coord[j], float)
            c2 = np.asarray(arr.coord[k], float)
            c3 = np.asarray(arr.coord[l], float)
            ic.bond_distances.setdefault((k, l),
                                         float(np.linalg.norm(c3 - c2)))
            v1 = c1 - c2
            v2 = c3 - c2
            cos = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
            ic.bond_angles.setdefault((j, k, l),
                                      float(np.degrees(np.arccos(
                                          np.clip(cos, -1, 1)))))
            ic.dihedra[quad] = dihedral(c0, c1, c2, c3)

        if quads is not None:
            # generic graph: explicit quads (legacy behaviour)
            for quad in quads:
                record(quad)
            if anchor is not None:
                ic.anchor = dict(anchor)
            else:
                ic.anchor = {0: tuple(np.asarray(arr.coord[0], float))}
            return ic

        def fill_anchor_geometry():
            """Record the anchor triple's own bonds (N-CA, CA-C) and its bond
            angle (N-CA-C) if not already covered by a grow quad.

            ``record`` stores ``bond_distances[(k, l)]`` and
            ``bond_angles[(j, k, l)]`` for every quad, so once the anchor atoms
            participate in any quad as parents these entries exist.  This fills
            the remaining holes so ``anchor`` is a fully-specified rigid body:
            the two anchor bonds and the one anchor angle are always queryable
            from the maps.
            """
            for i, j in ((nN, nCA), (nCA, nC)):
                if (i, j) not in ic.bond_distances:
                    ic.bond_distances[(i, j)] = float(np.linalg.norm(
                        np.asarray(arr.coord[j], float)
                        - np.asarray(arr.coord[i], float)))
            if (nN, nCA, nC) not in ic.bond_angles:
                v1 = np.asarray(arr.coord[nN], float) - np.asarray(arr.coord[nCA], float)
                v2 = np.asarray(arr.coord[nC], float) - np.asarray(arr.coord[nCA], float)
                cos = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
                ic.bond_angles[(nN, nCA, nC)] = float(np.degrees(
                    np.arccos(np.clip(cos, -1, 1))))

        # ---- protein-aware two-pass build --------------------------------
        # Group atoms into residues (chain_id, res_id, ins_code), preserving
        # atom order, and map atom name -> atom index within each residue.
        residues = {}      # key -> {"res_name": str, "atoms": {name: idx}}
        chain_keys = {}    # chain_id -> [keys in file order]
        for i in range(n):
            a = atoms[i]
            key = (a.chain_id, a.res_id, a.ins_code)
            if key not in residues:
                residues[key] = {"res_name": a.res_name.upper(), "atoms": {}}
                chain_keys.setdefault(a.chain_id, []).append(key)
            residues[key]["atoms"][a.name] = i

        auto_anchor = anchor is None
        ic.anchor = {} if auto_anchor else dict(anchor)

        for ckeys in chain_keys.values():
            for r_i, key in enumerate(ckeys):
                res = residues[key]["atoms"]
                name = residues[key]["res_name"]
                if not all(nm in res for nm in ("N", "CA", "C")):
                    continue          # incomplete residue: cannot extend chain
                nN, nCA, nC = (res[nm] for nm in ("N", "CA", "C"))

                # root frame of this chain = its first three backbone atoms
                if auto_anchor and r_i == 0:
                    ic.anchor.setdefault(
                        nN, tuple(np.asarray(arr.coord[nN], float)))
                    ic.anchor.setdefault(
                        nCA, tuple(np.asarray(arr.coord[nCA], float)))
                    ic.anchor.setdefault(
                        nC, tuple(np.asarray(arr.coord[nC], float)))

                # anchor triple must be a fully-specified rigid body: its own
                # N-CA / CA-C bonds and N-CA-C angle.  ``record`` already
                # covers them whenever a quad grows one of these atoms (the
                # carbonyl O branch or the peptide link); this fills any hole
                # (terminal residue / chain with no growth).
                fill_anchor_geometry()

                # carbonyl O (and C-terminal OXT) as branches off C --
                # the "intra" backbone grow quads, collected here so they can be
                # recorded with the rest of the backbone below.
                o_quads = []
                for spec in BACKBONE_IC_PATH["intra"]:
                    if all(nm in res for nm in spec):
                        o_quads.append(tuple(res[nm] for nm in spec))

                # peptide link to the next residue in the same chain; only
                # connect when the C_i - N_{i+1} distance is chemically
                # plausible (within the C-N bond-length upper bound), otherwise
                # the chain is broken here and we must not keep growing.
                c_n_ub = AMINO_ACID_BOND_LENGTH[("C", "N")]["up"]
                if r_i + 1 < len(ckeys):
                    nxt = residues[ckeys[r_i + 1]]["atoms"]
                    if all(nm in nxt for nm in ("N", "CA", "C")):
                        mN, mCA, mC = (nxt[nm] for nm in ("N", "CA", "C"))
                        c_n_dist = float(np.linalg.norm(
                            np.asarray(arr.coord[mN], float)
                            - np.asarray(arr.coord[nC], float)))
                        if c_n_dist <= c_n_ub:
                            # the "peptide" backbone grow quads: each grows
                            # one atom of residue i+1 from the frame spanning
                            # the peptide bond (see BACKBONE_IC_PATH)
                            def _bb(name):
                                if name.endswith("_i"):
                                    return res[name[:-2]]
                                if name.endswith("_{i+1}"):
                                    return nxt[name[:-6]]
                                return res[name]

                            for spec in BACKBONE_IC_PATH["peptide"]:
                                record(tuple(_bb(nm) for nm in spec))

                # carbonyl O (and C-terminal OXT) branch quads.  The read path
                # records the **measured** value for every quad, the O one
                # included: its job is to reproduce the input coordinates
                # exactly (``to_coords`` round-trip ~1e-14 A), and a real O sits
                # a few 0.01 A off the ideal peptide plane.  The constraint that
                # measured value satisfies in a real structure -- sp2
                # coplanarity of C, i.e. dihedral(N, CA, C, O) = psi - 180 -- is
                # defined once in topology.carbonyl_o_dihedral and used by the
                # *write* paths (template / build_side_chain), which have no O
                # to read.
                for quad in o_quads:
                    record(quad)

                # side chain: per-residue grow path (chi rotamers)
                for spec in IC_PATH.get(name, ()):
                    if all(nm in res for nm in spec):
                        record(tuple(res[nm] for nm in spec))
        return ic


class InternalCoord_AtomArray(Converter):
    """Rebuilds a biotite ``AtomArray`` (coordinates + annotations) from an
    :class:`InternalCoord`.

    The coordinates are grown from the anchor by
    :meth:`~biorazer.structure.objects.InternalCoord.to_coords`; the atom
    records supply the annotations.  Categories ``InternalCoordAtom`` does not
    carry (``hetero``, ``b_factor``, ``charge``, ...) come back at biotite's
    constructor defaults.

    The ``altloc_id`` annotation is restored **only if some record carries a
    label**; an all-empty one is left out, the way a hand-built ``AtomArray``
    has no such category.  Note that the label then still does not reach the
    *file*: biotite 1.6's PDB writer ignores the annotation (measured -- the
    altLoc column comes out blank), so it survives ``InternalCoord`` ->
    ``AtomArray`` only in memory.
    """

    def convert(self, tol=1e-6) -> AtomArray:
        """Rebuild the array from ``self.input_io``.

        Parameters
        ----------
        tol : float
            Tolerance in Angstrom for the redundant placements of a ring /
            any atom reachable through more than one dihedral (see
            :meth:`~biorazer.structure.objects.InternalCoord.to_coords`).

        Returns
        -------
        AtomArray
        """
        ic = self.input_io
        coords = ic.to_coords(tol=tol)
        n = len(ic.atoms)
        aa = AtomArray(n)
        aa.coord = np.array([coords[i] for i in range(n)], float)
        aa.chain_id = np.array([a.chain_id for a in ic.atoms])
        aa.res_name = np.array([a.res_name for a in ic.atoms], dtype="U3")
        aa.res_id = np.array([a.res_id for a in ic.atoms], dtype=np.int32)
        aa.atom_name = np.array([a.name for a in ic.atoms], dtype="U4")
        aa.element = np.array([a.element for a in ic.atoms])
        aa.ins_code = np.array([a.ins_code for a in ic.atoms])
        altloc = np.array([a.altloc_id for a in ic.atoms])
        if altloc.any():          # 全空就别建这个 category (与手工 AtomArray 一致)
            aa.set_annotation("altloc_id", altloc)
        return aa
