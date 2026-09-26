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
atoms : list[InternalCoordAtom]
    One entry per atom, carrying PDB-style annotations but no coordinates:
    ``(ins_code, chain_id, res_name, res_id, name, element, altloc_id)``.
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

A ``build_template`` template is a **single-residue** coordinate set.  The
parameters it was built from (``ss``, ``rotamer``, and the ss-class mean
``phi``/``psi``/``omega``) are deliberately **not** attributes of this
container: they describe a *residue build*, not a coordinate set, and a
multi-residue instance (built by
:class:`~biorazer.structure.bridge.atom_array.AtomArray_InternalCoord`, or by
``connect_internal_coords``) has no single value for them.
:func:`~biorazer.database.molecule.icoor.protein.template.build_template`
returns them separately as a
:class:`~biorazer.database.molecule.icoor.protein.template.TemplateSpec`.

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

This module holds the object and the operations that act on it alone.  Turning
a biotite ``AtomArray`` into an ``InternalCoord`` -- and back -- is a conversion
between two objects, so it lives in the bridge:
:class:`~biorazer.structure.bridge.atom_array.AtomArray_InternalCoord` (read
path) and
:class:`~biorazer.structure.bridge.atom_array.InternalCoord_AtomArray` (write
path).  Neither class carries a ``from_atomarray`` / ``to_atomarray`` method.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


def _place(B, C, D, blen, bang, dih_deg):
    """Place a new atom A given already-placed parents B, C, D.

    The frame is built **consistently** with :func:`dihedral`, so that
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
        :func:`dihedral` returns.

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


def dihedral(p0, p1, p2, p3):
    """Signed dihedral (degree) of 4 points p0..p3 (~N-CA-C-N).

    This is the repository's **single source of truth** for the sign
    convention of a dihedral angle: the returned value is the IUPAC-signed
    torsion about the ``p1-p2`` axis (the sign Dunbrack / Rosetta chi values
    use).  Every other module must call this function rather than carrying a
    private copy of the formula.

    It is exactly the angle :func:`_place` consumes/produces, i.e.
    ``dihedral(B, C, D, _place(B, C, D, ...))`` round-trips.
    """
    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return float(np.degrees(np.arctan2(y, x)))


#: Backwards-compatible alias (the function used to be private).
_dihedral = dihedral


@dataclass(repr=False, eq=False, slots=True)
class InternalCoordAtom:
    """A single atom's PDB-style annotation (no coordinates).

    Attribute access mirrors an ``AtomArray`` row: ``rec.chain_id = "B"`` etc.

    Fields
    ------
    ins_code : str
        PDB insertion code (``""`` when absent).
    chain_id : str
        Chain identifier.
    res_name : str
        Three-letter residue name, upper-case (e.g. ``"GLY"``).
    res_id : int
        Residue sequence number.
    name : str
        Atom name (``"N"``, ``"CA"``, ``"CB"``, ...).
    element : str | None
        Element symbol.  ``None`` (the default) derives it from ``name``:
        its first character when that is ``N``/``O``/``S``, otherwise ``"C"``.
    altloc_id : str
        Alternate-conformation label (the PDB altLoc column), empty when the
        atom has no alternate conformation.  Mirrors ``AtomArray.altloc_id``,
        including its biotite quirk: only an array read with
        ``altloc="all"`` carries that category at all (see
        :data:`NULL_ALT`).

    ``repr=False`` keeps the custom ``__repr__`` below; ``eq=False`` keeps
    **identity** equality and therefore hashability, because the record is
    *mutable* (the annotation setters of :class:`InternalCoord` write into
    it in place) -- a field-based ``__eq__`` would come with
    ``__hash__ = None`` and break any set/dict use.
    """

    ins_code: str = ""
    chain_id: str = "A"
    res_name: str = "GLY"
    res_id: int = 1
    name: str = "N"
    element: str | None = None
    altloc_id: str = ""

    def __post_init__(self):
        if self.element is None:
            self.element = (self.name[0] if self.name[0] in ("N", "O", "S")
                            else "C")

    def __repr__(self):
        return f"AtomRecord({self.chain_id}:{self.res_id}:{self.res_name}:{self.name})"


#: biotite spells "no alternate conformation" three ways depending on the
#: source: an empty PDB altLoc column is ``" "``, an mmCIF ``"."``, and a
#: missing value ``"?"``.  :class:`InternalCoordAtom` (and the selector's
#: field view) normalise all three to the empty string.
NULL_ALT = ("", " ", ".", "?")


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


@dataclass(repr=False, eq=False, slots=True)
class InternalCoord:
    """Internal-coordinate (generative) description of a structure.

    See the module docstring for the exact schema, the grow rule, and the
    meaning of each field.

    Fields
    ------
    atoms : list[InternalCoordAtom]
        The atom records (annotations only), in the numbering all the maps
        below are keyed by.
    anchor : dict[int, tuple[float, float, float]]
        ``{atom_index: (x, y, z)}`` -- the atoms that carry absolute
        coordinates (the growth roots).
    bond_distances : dict[tuple[int, int], float]
        ``{(i, j): Angstrom}``.
    bond_angles : dict[tuple[int, int, int], float]
        ``{(i, j, k): degree}`` -- angle at ``j``.
    dihedra : dict[tuple[int, int, int, int], float]
        ``{(i, j, k, l): degree}`` -- the generative map: ``l`` is grown from
        the parents ``(i, j, k)``.

    The container holds **geometry only** -- no build parameters.  The
    ``ss`` / ``rotamer`` / ss-mean ``phi``/``psi``/``omega`` a
    ``build_template`` call used are returned beside the template as a
    :class:`~biorazer.database.molecule.icoor.protein.template.TemplateSpec`,
    not attached here (see the module docstring).

    ``repr=False`` (a generated repr would dump every atom and map) and
    ``eq=False`` (the maps and atom records are mutated in place, so
    identity comparison plus hashability is kept) are deliberate.

    The per-atom annotation views below (``chain_id``, ``res_id``,
    ``res_name``, ``atom_name``, ``element``, ``ins_code``) are
    *properties*, not fields: they are numpy views over ``atoms`` for
    filtering (``ic.chain_id == "A"`` yields a mask).
    """

    atoms: list[InternalCoordAtom] = field(default_factory=list)
    anchor: dict[int, tuple[float, float, float]] = field(default_factory=dict)
    bond_distances: dict[tuple[int, int], float] = field(default_factory=dict)
    bond_angles: dict[tuple[int, int, int], float] = field(default_factory=dict)
    dihedra: dict[tuple[int, int, int, int], float] = field(
        default_factory=dict)

    def __post_init__(self):
        """Keep the pre-dataclass container contract.

        ``None`` means "empty" for every map (the old ``__init__`` accepted
        ``None`` for all five), and ``atoms`` is **copied** so a caller's
        list is never aliased by the instance.
        """
        self.atoms = [] if self.atoms is None else list(self.atoms)
        if self.anchor is None:
            self.anchor = {}
        if self.bond_distances is None:
            self.bond_distances = {}
        if self.bond_angles is None:
            self.bond_angles = {}
        if self.dihedra is None:
            self.dihedra = {}

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
        An alternate-conformation label is appended in parentheses:
        ``A:1:SER:CB(A)``.
        """
        a = self.atoms[i]
        assert isinstance(a, InternalCoordAtom)
        alt = f"({a.altloc_id})" if a.altloc_id else ""
        return f"{a.chain_id}:{a.res_id}:{a.res_name}:{a.name}{alt}"

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
    altloc_id = _annotation_accessor(
        "altloc_id", "altloc_id", None, str,
        "Alternate-conformation label of every atom as a numpy array (mirrors "
        "``AtomArray.altloc_id``; empty = no alternate conformation).  Set "
        "with a scalar or a length-``len(ic)`` sequence.")

    def dihedra_pd(self):
        """Dihedrals as a pandas table (easy filtering).

        Columns: ``i, j, k, l, dihedral, type`` -- the four atoms (as repr
        tags), the dihedral angle in **degree**, and the torsion ``type``
        annotated from the official definitions:

        * backbone ``phi`` / ``psi`` / ``omega`` --
          :data:`~biorazer.database.molecule.bond.dihedral.protein.ALIAS_QUAD`
          (IUPAC: ``phi = C_{i-1}-N_i-CA_i-C_i``, ``psi = N_i-CA_i-C_i-N_{i+1}``,
          ``omega = CA_i-C_i-N_{i+1}-CA_{i+1}``);
        * side chain ``chi1``..``chi4`` --
          :data:`~biorazer.database.molecule.bond.dihedral.protein.SIDECHAIN_CHI`
          (official Rosetta ``CHI`` rows, per residue).

        The read path
        (:class:`~biorazer.structure.bridge.atom_array.AtomArray_InternalCoord`)
        records backbone quads in the official atom order,
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
        (``ALIAS_QUAD`` / ``SIDECHAIN_CHI``); ``""`` when
        no official torsion matches.
        """
        from biorazer.database.molecule.bond.dihedral.protein import (
            ALIAS_QUAD,
            SIDECHAIN_CHI,
        )

        names = tuple(self.atoms[x].name for x in (i, j, k, l))
        for ttype, def_names in ALIAS_QUAD.items():
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
                assert newcoord is not None     # parents_list is never empty
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
