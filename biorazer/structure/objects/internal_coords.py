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
    ``(insert_code, chain_id, res_name, res_id, name, element)``.
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


class AtomRecord:
    """A single atom's PDB-style annotation (no coordinates).

    Attribute access mirrors an ``AtomArray`` row: ``rec.chain_id = "B"`` etc.
    """

    __slots__ = ("insert_code", "chain_id", "res_name", "res_id", "name",
                 "element")

    def __init__(self, insert_code="", chain_id="A", res_name="GLY", res_id=1,
                 name="N", element="N"):
        self.insert_code = insert_code
        self.chain_id = chain_id
        self.res_name = res_name
        self.res_id = res_id
        self.name = name
        self.element = element

    @classmethod
    def from_atom(cls, atom_array, index):
        arr = atom_array
        return cls(insert_code=str(arr.ins_code[index]),
                   chain_id=str(arr.chain_id[index]),
                   res_name=str(arr.res_name[index]),
                   res_id=int(arr.res_id[index]),
                   name=str(arr.atom_name[index]),
                   element=str(arr.element[index]))

    def __repr__(self):
        return f"AtomRecord({self.chain_id}:{self.res_id}:{self.res_name}:{self.name})"


class InternalCoord:
    """Internal-coordinate (generative) description of a structure.

    See the module docstring for the exact schema and the grow rule.
    """

    def __init__(self, atoms=None, anchor=None, bond_distances=None,
                 bond_angles=None, dihedra=None):
        self.atoms = atoms if atoms is not None else []
        self.anchor = anchor if anchor is not None else {}
        self.bond_distances = bond_distances if bond_distances is not None else {}
        self.bond_angles = bond_angles if bond_angles is not None else {}
        self.dihedra = dihedra if dihedra is not None else {}

    def __len__(self):
        return len(self.atoms)

    def __getitem__(self, i):
        return self.atoms[i]

    # ------------------------------------------------------------------ #
    #  atom representation + pandas-table views
    # ------------------------------------------------------------------ #
    def atom_repr(self, i):
        """Human/PDB-style tag for atom ``i`` (no coordinates).

        Format: ``{chain_id}:{res_id}:{res_name}:{name}``, e.g. ``A:1:SER:N``.
        """
        a = self.atoms[i]
        return f"{a.chain_id}:{a.res_id}:{a.res_name}:{a.name}"

    def dihedra_pd(self):
        """Dihedrals as a pandas table (easy filtering).

        Columns: ``i, j, k, l, dihedral`` -- the four atoms (as repr tags) and
        the dihedral angle in **degree**.
        """
        import pandas as pd

        rows = [
            (self.atom_repr(i), self.atom_repr(j), self.atom_repr(k),
             self.atom_repr(l), ang)
            for (i, j, k, l), ang in self.dihedra.items()
        ]
        return pd.DataFrame(rows, columns=["i", "j", "k", "l", "dihedral"])

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
        aa.ins_code = np.array([a.insert_code for a in self.atoms])
        return aa

    # ------------------------------------------------------------------ #
    #  construction from an AtomArray
    # ------------------------------------------------------------------ #
    @classmethod
    def from_atomarray(cls, arr, quads=None, anchor=None):
        """Build an :class:`InternalCoord` from a biotite ``AtomArray``.

        Parameters
        ----------
        arr : AtomArray
            Input structure.
        quads : list[tuple[int,int,int,int]] or None
            The (i,j,k,l) atom-index quadruples to record as dihedrals.  If
            ``None``, consecutive quadruples ``(0,1,2,3),(1,2,3,4),...`` are
            used (assumes a linear backbone already in polymer order).
        anchor : dict[int, tuple[float,float,float]] or None
            Absolute coordinates for the anchor atoms given as a mapping.
            Defaults to ``{0: (x,y,z)}`` (just the first atom).

        Notes
        -----
        Bond lengths and angles are derived from the input ``arr`` for the
        parent/child pairs of each quad, so they are exact (not idealised).
        """
        n = len(arr)
        atoms = [AtomRecord.from_atom(arr, i) for i in range(n)]
        ic = cls(atoms=atoms)
        if anchor is not None:
            ic.anchor = anchor
        else:
            ic.anchor = {0: tuple(np.asarray(arr.coord[0], float))}
        if quads is None:
            quads = [tuple(range(i, i + 4)) for i in range(max(0, n - 3))]
        for quad in quads:
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
        return ic