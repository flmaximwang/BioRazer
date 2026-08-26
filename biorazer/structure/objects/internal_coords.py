# -*- coding: utf-8 -*-
"""Internal-coordinate representation of a biomolecular structure.

Design (user-defined, 2026)
---------------------------
This is a *generative*, internal-coordinate (Z-matrix-like) description of a
structure: a list of atoms with their annotations, an *anchor* set of atoms
that carry absolute coordinates, and the bond lengths / bond angles /
dihedrals that let any atom be reconstructed from a few already-located
"parent" atoms.

Attributes
----------
atoms : list[AtomRecord]
    One entry per atom, carrying PDB-style annotations but no coordinates:
    ``(ins_code, chain_id, res_name, res_id, name, element)``.
anchor : dict[int, tuple[float,float,float]]
    ``{atom_index: (x,y,z)}``.  The anchor can be any atoms, but per the
    user's design it must be **3 consecutive atoms of a single dihedral**
    (positions 0,1,2 or 1,2,3 of that dihedral -- not a non-adjacent triple
    like 0,1,4) so that exactly one atom can be grown from them immediately.
bond_distances : dict[tuple[int,int], float]
    ``{(i,j): distance in Angstrom}`` -- bond length between atoms i and j.
bond_angles : dict[tuple[int,int,int], float]
    ``{(i,j,k): angle in degree}`` -- angle at atom j between i and k.
dihedra : dict[tuple[int,int,int,int], float]
    ``{(i,j,k,l): angle in degree}`` -- the dihedral of 4 ordered atoms.

All angles (``bond_angles`` and ``dihedra``) are in **degree**; only
``bond_distances`` is in Angstrom.

The reconstruction ("grow") rule for one target atom is:

* an atom :math:`l` whose dihedral ``(i,j,k,l)`` is present and whose parents
  ``i,j,k`` already have coordinates is placed using
  ``bond_distances[(k,l)]``, ``bond_angles[(j,k,l)]`` and ``dihedra[(i,j,k,l)]``
  as the Z-matrix triple (bond / angle at ``k`` / dihedral ``i-j-k-l``).

Because every dihedral contributes a directed dependency ``{i,j,k} -> l``,
the resulting graph must be a DAG for a unique reconstruction.  Rings / other
cycles make some atom reachable through more than one dihedral; in that case
the *independent* placements are computed and compared, and disagreement
raises ``ValueError`` (an inconsistent cycle), agreement passes (the ring
closes consistently).
"""

from __future__ import annotations

import numpy as np

from biorazer.structure.objects.external import AtomArray


def _place(B, C, D, blen, bang, dih_deg):
    """Place a new atom A given already-placed parents B, C, D.

    The frame is built **consistently** with :func:`_dihedral`, so that
    ``dihedral(B, C, D, A) == dih_deg`` exactly.  Frame vectors (k-centered,
    ``k == D``, i.e. B, C, D are the parent quads' i, j, k):

    * ``z = (C-D)`` direction ``k -> j``
    * ``x`` in the ``(B - C)`` direction, orthogonalised against z
    * ``y = cross(z, x)``

    Parameters
    ----------
    B, C, D : array_like (3,)
        Coordinates of the three parent atoms (the i, j, k of the dihedral).
    blen : float
        Bond length |A - D|.
    bang : float
        Bond angle at D between A and D and C (**degree**).
    dih_deg : float
        Dihedral angle ``(B, C, D, A)`` (**degree**), same convention as
        :func:`_dihedral` returns.

    Returns
    -------
    numpy.ndarray
        Coordinate of A.
    """
    B = np.asarray(B, float)
    C = np.asarray(C, float)
    D = np.asarray(D, float)
    z = D - C                       # direction k -> j
    z = z / np.linalg.norm(z)
    t = B - C
    t = t - np.dot(t, z) * z        # project out the z component
    if np.linalg.norm(t) < 1e-12:   # degenerate: pick any in-plane axis
        t = np.array([1.0, 0.0, 0.0], float)
        t = t - np.dot(t, z) * z
    bang = np.radians(bang)
    chi = np.radians(dih_deg)
    x = t / np.linalg.norm(t)
    y = np.cross(z, x)
    base = (-np.cos(bang) * z
            + np.sin(bang) * (np.cos(chi) * x + np.sin(chi) * y))
    return D + blen * base


def _dihedral(p0, p1, p2, p3):
    """Signed dihedral (degree) of 4 points p0..p3 (~N-CA-C-N)."""
    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return float(np.degrees(np.arctan2(y, x)))


class InternalCoordAtom:
    """A single atom's PDB-style annotation (no coordinates).

    Attribute access mirrors an ``AtomArray`` row: ``rec.chain_id = "B"`` etc.
    """

    __slots__ = ("ins_code", "chain_id", "res_name", "res_id", "name",
                 "element")

    def __init__(self, ins_code="", chain_id="A", res_name="GLY", res_id=1,
                 name="N", element="N"):
        self.ins_code = ins_code
        self.chain_id = chain_id
        self.res_name = res_name
        self.res_id = res_id
        self.name = name
        self.element = element

    @classmethod
    def from_atom(cls, atom_array, index):
        arr = atom_array
        return cls(ins_code=str(arr.ins_code[index]),
                   chain_id=str(arr.chain_id[index]),
                   res_name=str(arr.res_name[index]),
                   res_id=int(arr.res_id[index]),
                   name=str(arr.atom_name[index]),
                   element=str(arr.element[index]))

    def __repr__(self):
        return f"AtomRecord({self.chain_id}:{self.res_id}:{self.res_name}:{self.name})"


def _annotation_accessor(rec_attr, name, dtype, cast, doc):
    """Build a property+setter mapping an atom annotation to a numpy array view.

    Mirrors ``AtomArray.<name>`` so filtering on an ``InternalCoord`` reads
    naturally (``mask = ic.chain_id == "A"`` yields a boolean array over all
    atoms).  The getter returns a numpy array; the setter accepts a scalar
    (broadcast to every atom) or a length-``len(ic)`` sequence/array.

    Parameters
    ----------
    rec_attr : str
        Attribute name on each :class:`InternalCoordAtom`.
    name : str
        Public accessor name (for error messages).
    dtype : numpy dtype or None
        Array dtype for the getter.  ``None`` -> a per-atom-width ``U`` array.
    cast : callable
        Applied to each value on both get and set (e.g. ``int`` for ``res_id``).
    doc : str
        Docstring for the property.
    """

    def getter(self):
        if not self.atoms:
            return np.array([], dtype="U1" if dtype is None else dtype)
        if dtype is None:
            width = max(len(str(getattr(a, rec_attr))) for a in self.atoms)
            return np.array([cast(getattr(a, rec_attr)) for a in self.atoms],
                            dtype=f"U{width}")
        return np.array([cast(getattr(a, rec_attr)) for a in self.atoms],
                        dtype=dtype)

    def setter(self, value):
        n = len(self.atoms)
        if isinstance(value, (str, bytes)):
            vals = [value] * n
        else:
            vals = list(value)
            if len(vals) != n:
                raise ValueError(
                    f"{name} must be a scalar or a sequence of length {n}, "
                    f"got {len(vals)}")
        for a, v in zip(self.atoms, vals):
            setattr(a, rec_attr, cast(v))

    return property(getter, setter, doc=doc)


class InternalCoord:
    """Internal-coordinate (generative) description of a structure.

    See the module docstring for the exact schema and the grow rule.
    """

    def __init__(self, atoms=None, anchor=None, bond_distances=None,
                 bond_angles=None, dihedra=None):
        self.atoms = list(atoms) if atoms is not None else []
        self.anchor = anchor if anchor is not None else {}
        self.bond_distances = bond_distances if bond_distances is not None else {}
        self.bond_angles = bond_angles if bond_angles is not None else {}
        self.dihedra = dihedra if dihedra is not None else {}

    def __len__(self):
        return len(self.atoms)

    def __getitem__(self, key):
        """Select atoms by integer, slice, or boolean / index array.

        * ``ic[i]`` (int) -> the ``AtomRecord`` at rank ``i`` (as before).
        * ``ic[mask]`` (boolean array of length ``len(ic)``), ``ic[idx]``
          (integer index array), or ``ic[slice]`` -> **a new ``InternalCoord``**
          restricted to the selected atoms.  Its connectivity maps
          (``dihedra``/``bond_angles``/``bond_distances``) and ``anchor`` are
          reindexed to the new 0-based numbering; any map entries touching
          atoms outside the selection are dropped.
        """
        if isinstance(key, (int, np.integer)):
            return self.atoms[int(key)]
        if isinstance(key, slice):
            sel = list(range(*key.indices(len(self.atoms))))
            return self._subset(sel)
        key = np.asarray(key)
        if key.dtype == bool:
            if key.ndim != 1 or len(key) != len(self.atoms):
                raise ValueError(
                    f"boolean mask must be 1-D of length {len(self.atoms)}, "
                    f"got shape {key.shape}")
            sel = np.nonzero(key)[0].tolist()
            return self._subset(sel)
        if key.dtype.kind in "iu":
            return self._subset([int(i) for i in key])
        raise TypeError(
            f"index must be int, slice, boolean mask, or integer array; "
            f"got {type(key).__name__}")

    def _subset(self, sel):
        """New :class:`InternalCoord` restricted to atom ranks ``sel``.

        Connectivity maps and the anchor are reindexed to the new 0-based
        numbering; entries touching atoms outside ``sel`` are dropped.
        """
        remap = {old: new for new, old in enumerate(sel)}
        new_ic = type(self)(atoms=[self.atoms[i] for i in sel])
        for (i, j), d in self.bond_distances.items():
            if i in remap and j in remap:
                new_ic.bond_distances[(remap[i], remap[j])] = d
        for (i, j, k), ang in self.bond_angles.items():
            if i in remap and j in remap and k in remap:
                new_ic.bond_angles[(remap[i], remap[j], remap[k])] = ang
        for (i, j, k, l), dih in self.dihedra.items():
            if all(x in remap for x in (i, j, k, l)):
                new_ic.dihedra[(remap[i], remap[j], remap[k], remap[l])] = dih
        new_ic.anchor = {remap[i]: c for i, c in self.anchor.items()
                         if i in remap}
        return new_ic

    # ------------------------------------------------------------------ #
    #  connectivity (runtime, derived from dihedra)
    # ------------------------------------------------------------------ #
    @property
    def fragments(self):
        """Split all atoms into connected components (subgraphs) of ``dihedra``.

        Two atoms are in the same fragment iff they are connected through a
        path of ``dihedra`` edges (each dihedral ``(i, j, k, l)`` contributes
        the edges ``i-j``, ``j-k``, ``k-l``).  This is computed at runtime from
        ``self.dihedra`` only; atoms that appear in no dihedral form their own
        single-atom fragments.

        Returns
        -------
        list[list[int]]
            The connected components, each a sorted list of atom indices.
        """
        n = len(self.atoms)

        adj = [[] for _ in range(n)]
        for (i, j, k, l) in self.dihedra:
            adj[i].append(j)
            adj[j].append(i)
            adj[j].append(k)
            adj[k].append(j)
            adj[k].append(l)
            adj[l].append(k)

        seen = [False] * n
        fragments = []
        for start in range(n):
            if seen[start]:
                continue
            comp = []
            stack = [start]
            seen[start] = True
            while stack:
                x = stack.pop()
                comp.append(x)
                for y in adj[x]:
                    if not seen[y]:
                        seen[y] = True
                        stack.append(y)
            fragments.append(sorted(comp))
        return fragments

    # ------------------------------------------------------------------ #
    #  atom representation + pandas-table views
    # ------------------------------------------------------------------ #
    def atom_repr(self, i):
        """Human/PDB-style tag for atom ``i`` (no coordinates).

        Format: ``{chain_id}:{res_id}:{res_name}:{name}``, e.g. ``A:1:SER:N``.
        """
        a = self.atoms[i]
        assert isinstance(a, InternalCoordAtom)
        return f"{a.chain_id}:{a.res_id}:{a.res_name}:{a.name}"

    # Per-atom annotation views, mirroring ``AtomArray``.  Each is a
    # property+setter: ``ic.chain_id`` returns a numpy array over all atoms
    # (so ``ic.chain_id == "A"`` yields a mask), and assigning broadcasts a
    # scalar or takes a length-``len(ic)`` sequence.
    chain_id = _annotation_accessor(
        "chain_id", "chain_id", None, str,
        'Chain ID of every atom as a numpy array; `ic.chain_id == "A"` '
        "yields a mask.  Set with a scalar or a length-`len(ic)` sequence.")
    res_id = _annotation_accessor(
        "res_id", "res_id", np.int32, int,
        "Residue ID of every atom as an int32 numpy array.  Set with a scalar "
        "or a length-``len(ic)`` sequence.")
    res_name = _annotation_accessor(
        "res_name", "res_name", "U3", str,
        "Residue name of every atom as a ``U3`` numpy array.  Set with a "
        "scalar or a length-``len(ic)`` sequence.")
    atom_name = _annotation_accessor(
        "name", "atom_name", "U4", str,
        "Atom name of every atom as a ``U4`` numpy array (mirrors "
        "``AtomArray.atom_name``).  Set with a scalar or a length-``len(ic)`` "
        "sequence.")
    element = _annotation_accessor(
        "element", "element", None, str,
        "Element of every atom as a numpy array.  Set with a scalar or a "
        "length-``len(ic)`` sequence.")
    ins_code = _annotation_accessor(
        "ins_code", "ins_code", None, str,
        "Insertion code of every atom as a numpy array (mirrors "
        "``AtomArray.ins_code``).  Set with a scalar or a length-``len(ic)`` "
        "sequence.")

    def dihedra_pd(self):
        """Dihedrals as a pandas table (easy filtering).

        Columns: ``i, j, k, l, dihedral, type`` -- the four atoms (as repr
        tags), the dihedral angle in **degree**, and the torsion ``type``
        annotated from the official definitions:

        * backbone ``phi`` / ``psi`` / ``omega`` --
          :data:`~biorazer.database.molecule.bond.dihedral.protein.MAINCHAIN_TORSION_DEFINITIONS`
          (IUPAC: ``phi = C_{i-1}-N_i-CA_i-C_i``, ``psi = N_i-CA_i-C_i-N_{i+1}``,
          ``omega = CA_i-C_i-N_{i+1}-CA_{i+1}``);
        * side chain ``chi1``..``chi4`` --
          :data:`~biorazer.database.molecule.bond.dihedral.protein.SIDECHAIN_CHI`
          (official Rosetta ``CHI`` rows, per residue).

        ``from_atomarray`` records backbone quads in the official atom order,
        so ``phi``/``psi``/``omega`` annotate directly; side-chain quads are
        stored in the official ICOOR order (bonded parent in slot ``k``), so
        the first quads equal the official chi definitions.  Dihedrals that
        match no official torsion (e.g. the carbonyl ``O`` branch
        ``(N, CA, C, O)``) get an empty ``type``.
        """
        import pandas as pd

        rows = [
            (self.atom_repr(i), self.atom_repr(j), self.atom_repr(k),
             self.atom_repr(l), ang, self._torsion_type(i, j, k, l))
            for (i, j, k, l), ang in self.dihedra.items()
        ]
        return pd.DataFrame(rows, columns=["i", "j", "k", "l", "dihedral",
                                           "type"])

    def _torsion_type(self, i, j, k, l):
        """Official torsion name for the dihedral ``(i, j, k, l)``.

        ``phi`` / ``psi`` / ``omega`` / ``chi1``..``chi4`` by **exact
        atom-name sequence** against the official definitions
        (``MAINCHAIN_TORSION_DEFINITIONS`` / ``SIDECHAIN_CHI``); ``""`` when
        no official torsion matches.
        """
        from biorazer.database.molecule.bond.dihedral.protein import (
            MAINCHAIN_TORSION_DEFINITIONS,
            SIDECHAIN_CHI,
        )

        names = tuple(self.atoms[x].name for x in (i, j, k, l))
        for ttype, def_names in MAINCHAIN_TORSION_DEFINITIONS.items():
            if names == def_names:
                return ttype
        resn = self.atoms[i].res_name.upper()
        for n, chi_quad in enumerate(SIDECHAIN_CHI.get(resn, [])):
            if names == chi_quad:
                return f"chi{n + 1}"
        return ""

    def bond_distances_pd(self):
        """Bond lengths as a pandas table (easy filtering).

        Columns: ``i, j, distance`` -- the two bonded atoms (repr tags) and
        the distance in Angstrom.
        """
        import pandas as pd

        rows = [(self.atom_repr(i), self.atom_repr(j), d)
                for (i, j), d in self.bond_distances.items()]
        return pd.DataFrame(rows, columns=["i", "j", "distance"])

    def bond_angles_pd(self):
        """Bond angles as a pandas table (easy filtering).

        Columns: ``i, j, k, angle`` -- the three atoms (repr tags) and the
        angle at atom ``j`` in **degree**.
        """
        import pandas as pd

        rows = [(self.atom_repr(i), self.atom_repr(j), self.atom_repr(k), ang)
                for (i, j, k), ang in self.bond_angles.items()]
        return pd.DataFrame(rows, columns=["i", "j", "k", "angle"])

    # ------------------------------------------------------------------ #
    #  growth machinery
    # ------------------------------------------------------------------ #
    def _deps(self):
        """Map each target atom l -> list of (parents, quad) it depends on.

        Every dihedral ``(i,j,k,l)`` adds a directed dependency
        ``{i,j,k} -> l``.  A ring makes one ``l`` appear from several quads.
        """
        deps = {}
        for quad in self.dihedra:
            i, j, k, l = quad
            entry = (tuple(quad[:3]), quad)
            deps.setdefault(l, []).append(entry)
        return deps

    def to_coords(self, tol=1e-6):
        """Reconstruct all coordinates from the anchor by graph growth.

        Returns a dict ``{atom_index: (x,y,z)}``.

        Raises
        ------
        ValueError
            * If some atoms are left unreachable (parents never located), or
            * on an inconsistent cycle: an atom reached from two dihedrals
              that disagree beyond ``tol``.
        """
        coords = dict(self.anchor)
        placed = set(coords)
        deps = self._deps()
        remain = dict(deps)
        progress = True
        while remain and progress:
            progress = False
            for l, parents_list in list(remain.items()):
                if not all(all(p in placed for p in parents)
                           for parents, _ in parents_list):
                    continue
                newcoord = None
                for parents, quad in parents_list:
                    i, j, k, l = quad
                    blen = self.bond_distances[(k, l)]
                    bang = self.bond_angles[(j, k, l)]
                    dih = self.dihedra[quad]
                    pos = _place(coords[i], coords[j], coords[k], blen, bang, dih)
                    if newcoord is None:
                        newcoord = pos
                    elif np.linalg.norm(pos - newcoord) > tol:
                        raise ValueError(
                            f"Inconsistent cycle at atom {l}: {quad} gives "
                            f"{np.round(newcoord, 4)} vs {np.round(pos, 4)}")
                if l in placed:
                    if np.linalg.norm(newcoord - coords[l]) > tol:
                        raise ValueError(f"Inconsistent coordinate for atom {l}")
                else:
                    coords[l] = newcoord
                    placed.add(l)
                del remain[l]
                progress = True
        if remain:
            raise ValueError(
                f"Unreachable atoms (parents never located): {sorted(remain)}")
        return coords

    def to_atomarray(self, tol=1e-6):
        """Rebuild to a biotite ``AtomArray`` (coordinates + annotations)."""
        coords = self.to_coords(tol=tol)
        n = len(self.atoms)
        aa = AtomArray(n)
        aa.coord = np.array([coords[i] for i in range(n)], float)
        aa.chain_id = np.array([a.chain_id for a in self.atoms])
        aa.res_name = np.array([a.res_name for a in self.atoms], dtype="U3")
        aa.res_id = np.array([a.res_id for a in self.atoms], dtype=np.int32)
        aa.atom_name = np.array([a.name for a in self.atoms], dtype="U4")
        aa.element = np.array([a.element for a in self.atoms])
        aa.ins_code = np.array([a.ins_code for a in self.atoms])
        return aa

    # ------------------------------------------------------------------ #
    #  construction from an AtomArray
    # ------------------------------------------------------------------ #
    @classmethod
    def from_atomarray(cls, arr, quads=None, anchor=None):
        """Build an :class:`InternalCoord` from a biotite ``AtomArray``.

        Protein-aware construction (**the default**, ``quads=None``) runs a
        two-pass build:

        * **Main-chain pass** (uniform, same for every residue): the backbone
          ``N -> CA -> C -> O`` is walked residue by residue, linking residue
          ``i``'s ``C`` to residue ``i+1``'s ``N`` (peptide bond).  It records
          the cross-residue quads ``(N_i, CA_i, C_i, N_{i+1})``,
          ``(CA_i, C_i, N_{i+1}, CA_{i+1})``, ``(C_i, N_{i+1}, CA_{i+1}, C_{i+1})``
          and the per-residue carbonyl branch ``(N, CA, C, O)``.
        * **Side-chain pass** (per residue): each standard amino acid's side
          chain is grown off the already-placed backbone using its per-residue
          grow-path table ``IC_PATH`` (chi rotamers; see
          ``biorazer.database.molecule.icoor.protein.topology``).  Non-standard /
          non-protein atoms (water, ligands, hydrogens) are not covered.

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

        Parameters
        ----------
        arr : AtomArray
            Input structure.
        quads : list[tuple[int,int,int,int]] or None
            Explicit (i,j,k,l) atom-index quadruples for the general-graph
            path.  If ``None`` (default) the protein-aware two-pass build above
            is used.
        anchor : dict[int, tuple[float,float,float]] or None
            Absolute coordinates for anchor atoms given as a mapping.  If
            ``None`` (default) anchors are auto-detected (per chain's first
            ``N, CA, C`` for the protein path; ``{0: first atom}`` for the
            general ``quads`` path).

        Notes
        -----
        Bond lengths and angles are derived from the input ``arr`` for the
        parent/child pairs of each quad, so they are exact (not idealised).
        """
        from biorazer.database.molecule.icoor.protein.topology import IC_PATH
        from biorazer.database.molecule.bond.length.generic import AMINO_ACID_BOND_LENGTH

        n = len(arr)
        atoms = [InternalCoordAtom.from_atom(arr, i) for i in range(n)]
        ic = cls(atoms=atoms)

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
            ic.dihedra[quad] = _dihedral(c0, c1, c2, c3)

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

                # carbonyl O (and C-terminal OXT) as branches off C
                if "O" in res:
                    record((nN, nCA, nC, res["O"]))
                if "OXT" in res:
                    record((nN, nCA, nC, res["OXT"]))

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
                            record((nN, nCA, nC, mN))    # C_i - N_{i+1} peptide
                            record((nCA, nC, mN, mCA))   # N_{i+1} - CA_{i+1}
                            record((nC, mN, mCA, mC))    # CA_{i+1} - C_{i+1}

                # side chain: per-residue grow path (chi rotamers)
                for spec in IC_PATH.get(name, ()):
                    if all(nm in res for nm in spec):
                        record(tuple(res[nm] for nm in spec))
        return ic