"""Manipulation of :class:`~biorazer.structure.objects.internal_coords.InternalCoord`.

The operand of this module is biorazer's own internal-coordinate container
(:mod:`biorazer.structure.objects.internal_coords`), not biotite's
``AtomArray`` -- the AtomArray operations live in
:mod:`biorazer.structure.manipulation.atom_array`.
"""

import numpy as np

from biorazer.database.molecule.bond.angle.generic import (
    AMINO_ACID_BACKBONE_BOND_ANGLE,
)
from biorazer.database.molecule.bond.dihedral.protein import OMEGA_TRANS
from biorazer.database.molecule.bond.length.protein import AMINO_ACID_BOND_LENGTH

from ...objects import InternalCoord


def _residue_ca_index(ic, index, expected_name, param_name):
    """Index of the CA atom in the residue containing ``index``.

    ``index`` must be a valid atom index whose atom name equals
    ``expected_name`` (so the seam geometry is built around the right
    backbone atoms, catching off-by-N index bugs early).  Returns the index
    of the ``CA`` atom of that same residue.
    """
    if not isinstance(index, (int, np.integer)):
        raise TypeError(
            f"{param_name} must be an integer atom index, "
            f"got {type(index).__name__}")
    index = int(index)
    n = len(ic.atoms)
    if not 0 <= index < n:
        raise ValueError(
            f"{param_name} must be a valid atom index in [0, {n - 1}], "
            f"got {index}")
    target = ic.atoms[index]
    if target.name != expected_name:
        raise ValueError(
            f"atom at {param_name}={index} is {target.name!r} "
            f"({ic.atom_repr(index)}); expected {expected_name!r}")
    for i, a in enumerate(ic.atoms):
        if (a.chain_id == target.chain_id
                and a.res_id == target.res_id
                and a.ins_code == target.ins_code
                and a.name == "CA"):
            return i
    raise ValueError(
        f"residue containing {ic.atom_repr(index)} has no CA atom; "
        f"cannot define the peptide-bond angles / omega")


def connect_internal_coords(
    N_terminal_ic: InternalCoord,
    C_terminal_ic: InternalCoord,
    C_index: int,
    N_index: int,
    bond_length: float | None = None,
    angle_CA_C_N: float | None = None,
    angle_C_N_CA: float | None = None,
    omega: float | None = None,
) -> InternalCoord:
    """Connect two :class:`InternalCoord` fragments with a new peptide bond.

    Joins the C-terminal end of ``N_terminal_ic`` to the N-terminal end of
    ``C_terminal_ic`` by a **new peptide bond**: the carbonyl carbon at
    ``C_index`` (of the N-terminal fragment, which is fragment's ``C``) is
    bonded to the amide nitrogen at ``N_index`` (of the C-terminal fragment,
    its ``N``).  The two ``CA`` atoms of the terminal residues are located
    automatically from the residue context of ``C_index`` / ``N_index``.

    The result is a **new** ``InternalCoord`` that is the concatenation of
    the two fragments (``N_terminal_ic`` first, then ``C_terminal_ic``), with
    the new peptide-bond geometry recorded at the seam:

    * ``bond_distances[(C_index, N_index)]`` = the ``C--N`` bond length,
    * ``bond_angles[(CA_n, C_index, N_index)]`` = the angle at ``C`` between
      the N-terminal residue's ``CA`` and the new ``N``,
    * ``bond_angles[(C_index, N_index, CA_c)]`` = the angle at ``N`` between
      the new ``C`` and the C-terminal residue's ``CA``,
    * ``dihedra[(CA_n, C_index, N_index, CA_c)]`` = ``omega`` across the
      peptide plane (``CA_n``/``CA_c`` = the two terminal residues' ``CA``).

    Both fragments' ``anchor`` dictionaries are preserved and merged, so the
    returned object reconstructs both fragments in place (a single connected
    component via the seam dihedral).  Because the C-terminal fragment keeps
    its own anchors, :meth:`InternalCoord.to_coords` checks the seam against
    the requested geometry: if the C-terminal fragment was not pre-positioned
    to sit at ``omega``/the two angles it raises ``ValueError``
    ("Inconsistent coordinate").  Callers that want to *place* the C-terminal
    fragment at an ideal seam must first orient it (e.g. rotate about the new
    ``C--N`` bond to set ``omega``) before building the fragments.

    Units follow :mod:`biorazer.structure.objects.internal_coords`:
    ``bond_length`` in Angstrom; ``angle_CA_C_N``, ``angle_C_N_CA`` and
    ``omega`` in **degree**.

    Parameters
    ----------
    N_terminal_ic : InternalCoord
        N-terminal fragment; provides the carbonyl carbon ``C_index``.
    C_terminal_ic : InternalCoord
        C-terminal fragment; provides the amide nitrogen ``N_index``.
    C_index : int
        Atom index (in ``N_terminal_ic``) of the carbonyl carbon ``C`` to
        join.
    N_index : int
        Atom index (in ``C_terminal_ic``) of the amide nitrogen ``N`` to join.
    bond_length : float, optional
        New ``C--N`` peptide bond length in Angstrom.  Default is the Engh &
        Huber mean ``AMINO_ACID_BOND_LENGTH[("C", "N")]["mean"]`` (1.329 A).
    angle_CA_C_N : float, optional
        Angle at ``C`` between the N-terminal residue's ``CA`` and the new
        ``N``, in degree.  Default ``AMINO_ACID_BOND_ANGLE[("CA","C","N")]``
        mean (116.2).
    angle_C_N_CA : float, optional
        Angle at ``N`` between the new ``C`` and the C-terminal residue's
        ``CA``, in degree.  Default ``AMINO_ACID_BOND_ANGLE[("C","N","CA")]``
        mean (121.7).
    omega : float, optional
        Dihedral ``(CA, C, N, CA)`` across the new peptide bond, in degree.
        Default ``OMEGA_TRANS["mean"]`` (180.0, trans).

    Returns
    -------
    InternalCoord
        A new merged ``InternalCoord`` with the recorded peptide seam.

    Raises
    ------
    ValueError
        If an index is out of range or its atom is not the expected ``C`` /
        ``N``, or if the terminal residue has no ``CA`` atom.
    """
    if bond_length is None:
        bond_length = AMINO_ACID_BOND_LENGTH[("C", "N")]["mean"]
    if angle_CA_C_N is None:
        angle_CA_C_N = AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "N")]["mean"]
    if angle_C_N_CA is None:
        angle_C_N_CA = AMINO_ACID_BACKBONE_BOND_ANGLE[("C", "N", "CA")]["mean"]
    if omega is None:
        omega = OMEGA_TRANS["mean"]

    # Terminal-residue CA atoms, used by both angles and omega.
    ca_nt = _residue_ca_index(N_terminal_ic, C_index, "C", "C_index")
    ca_ct = _residue_ca_index(C_terminal_ic, N_index, "N", "N_index")

    offset = len(N_terminal_ic)

    out = InternalCoord()
    out.atoms = list(N_terminal_ic.atoms) + list(C_terminal_ic.atoms)

    # Merge connectivity maps, reindexing the C-terminal fragment's entries.
    out.bond_distances = dict(N_terminal_ic.bond_distances)
    for (i, j), d in C_terminal_ic.bond_distances.items():
        out.bond_distances[(i + offset, j + offset)] = d
    out.bond_angles = dict(N_terminal_ic.bond_angles)
    for (i, j, k), a in C_terminal_ic.bond_angles.items():
        out.bond_angles[(i + offset, j + offset, k + offset)] = a
    out.dihedra = dict(N_terminal_ic.dihedra)
    for (i, j, k, l), d in C_terminal_ic.dihedra.items():
        out.dihedra[(i + offset, j + offset, k + offset, l + offset)] = d

    # Preserve both fragments' anchors so each reconstructs in place.
    out.anchor = dict(N_terminal_ic.anchor)
    for i, (x, y, z) in C_terminal_ic.anchor.items():
        out.anchor[i + offset] = (x, y, z)

    # The new peptide seam.
    c = C_index                   # carbonyl C of the N-terminal fragment
    n = N_index + offset          # amide N of the C-terminal fragment

    # The seam dihedral (CA, C, N, CA) grows the C-terminal residue's CA,
    # which needs that residue's N-CA bond.  The read path
    # (``AtomArray_InternalCoord``) omits the first residue's N-CA (it seeds
    # the anchor), so fill it from the anchor coordinates when missing.
    if (n, ca_ct + offset) not in out.bond_distances:
        n_xyz = out.anchor.get(n)
        ca_xyz = out.anchor.get(ca_ct + offset)
        if n_xyz is not None and ca_xyz is not None:
            out.bond_distances[(n, ca_ct + offset)] = float(
                np.linalg.norm(np.asarray(n_xyz, float)
                               - np.asarray(ca_xyz, float)))
        else:
            raise ValueError(
                "C-terminal fragment is missing the N-CA bond geometry "
                "of its N-terminal residue and its N/CA are not anchored; "
                "cannot build a reconstructible seam")

    out.bond_distances[(c, n)] = bond_length
    out.bond_angles[(ca_nt, c, n)] = angle_CA_C_N
    out.bond_angles[(c, n, ca_ct + offset)] = angle_C_N_CA
    out.dihedra[(ca_nt, c, n, ca_ct + offset)] = omega
    return out
