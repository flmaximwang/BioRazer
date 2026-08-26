"""Helpers for writing PDB record lines (LINK / SSBOND / SEGID).

These are the low-level formatters used by :class:`AtomArray_Pdb`
(``biorazer.structure.io.protein``) to emit intermolecular LINK records,
disulfide SSBOND records and the SEGID annotation columns that biotite's own
PDB writer does not produce.
"""

import numpy as np
from biotite.structure import AtomArray
import biotite.structure as bio_struc
from biotite.structure.io.pdb.hybrid36 import encode_hybrid36


# Same limit as biotite's PDB writer (_PDB_MAX_RESIDUES), so that the
# residue sequence numbers in LINK/SSBOND match the ATOM records.
_PDB_MAX_RESIDUES = 9999


def _format_pdb_atom_name(atom_name: str, element: str) -> str:
    """
    Format a biotite atom name into the 4-column PDB atom-name field
    (names of atoms with a one-letter element are shifted right by one).
    """
    name = atom_name.strip()
    if len(element) == 1 and len(name) < 4:
        name = " " + name
    return name.ljust(4)


def _format_pdb_res_id(res_id: int, hybrid36: bool) -> str:
    """Format a residue sequence number exactly like the ATOM records."""
    if hybrid36:
        return encode_hybrid36(res_id, 4)
    if res_id > 0:
        res_id = ((res_id - 1) % _PDB_MAX_RESIDUES) + 1
    return f"{res_id:>4}"


def _format_link_record(array: AtomArray, protein_i: int, ligand_j: int, hybrid36: bool) -> str:
    """Build a LINK record for a covalent bond between a protein atom and a ligand atom."""
    dist = np.linalg.norm(array.coord[protein_i] - array.coord[ligand_j])
    return (
        "LINK  "
        + " " * 6
        + _format_pdb_atom_name(array.atom_name[protein_i], array.element[protein_i])
        + " "
        + f"{array.res_name[protein_i]:>3}"
        + " "
        + f"{array.chain_id[protein_i]:<1}"
        + _format_pdb_res_id(array.res_id[protein_i], hybrid36)
        + f"{array.ins_code[protein_i]:<1}"
        + " " * 15
        + _format_pdb_atom_name(array.atom_name[ligand_j], array.element[ligand_j])
        + " "
        + f"{array.res_name[ligand_j]:>3}"
        + " "
        + f"{array.chain_id[ligand_j]:<1}"
        + _format_pdb_res_id(array.res_id[ligand_j], hybrid36)
        + f"{array.ins_code[ligand_j]:<1}"
        + "  "
        + f"{'1555':>6}"
        + " "
        + f"{'1555':>6}"
        + " "
        + f"{dist:>5.2f}"
        + "  "
    )


def _format_ssbond_record(
    array: AtomArray, cys_i: int, cys_j: int, serial: int, hybrid36: bool
) -> str:
    """Build an SSBOND record for a CYS SG - CYS SG disulfide bond."""
    dist = np.linalg.norm(array.coord[cys_i] - array.coord[cys_j])
    return (
        "SSBOND"
        + " "
        + f"{serial:>3}"
        + " "
        + "CYS"
        + " "
        + f"{array.chain_id[cys_i]:<1}"
        + " "
        + _format_pdb_res_id(array.res_id[cys_i], hybrid36)
        + f"{array.ins_code[cys_i]:<1}"
        + "   "
        + "CYS"
        + " "
        + f"{array.chain_id[cys_j]:<1}"
        + " "
        + _format_pdb_res_id(array.res_id[cys_j], hybrid36)
        + f"{array.ins_code[cys_j]:<1}"
        + " " * 23
        + f"{'1555':>6}"
        + " "
        + f"{'1555':>6}"
        + " "
        + f"{dist:>5.2f}"
        + "  "
    )


def _format_link_records(array: AtomArray, hybrid36: bool) -> list[str]:
    """Find covalent bonds between a canonical amino acid and a small molecule."""
    protein_mask = bio_struc.filter_canonical_amino_acids(array)
    ligand_mask = ~protein_mask & ~bio_struc.filter_solvent(array)
    records = []
    for atom_i, atom_j, _ in array.bonds.as_array():
        if protein_mask[atom_i] and ligand_mask[atom_j]:
            records.append(_format_link_record(array, atom_i, atom_j, hybrid36))
        elif protein_mask[atom_j] and ligand_mask[atom_i]:
            records.append(_format_link_record(array, atom_j, atom_i, hybrid36))
    return records


def _format_ssbond_records(array: AtomArray, hybrid36: bool) -> list[str]:
    """Find disulfide bonds (SG of two CYS residues) and format them as SSBOND."""
    cys_sg = (array.res_name == "CYS") & (np.char.strip(array.atom_name) == "SG")
    records = []
    for atom_i, atom_j, _ in array.bonds.as_array():
        if atom_i != atom_j and cys_sg[atom_i] and cys_sg[atom_j]:
            records.append(
                _format_ssbond_record(array, atom_i, atom_j, len(records) + 1, hybrid36)
            )
    return records


def _inject_seg_ids(lines: list[str], array: AtomArray) -> None:
    """
    Write the ``seg_id`` annotation into the SEGID columns (73-76) of the
    ATOM/HETATM lines.

    Why this is needed: biotite's PDB writer has no ``seg_id`` support (no
    such annotation category; columns 73-76 are always blank), but SEGID is
    the standard field for grouping atoms by segment (e.g. repeat units of
    a fiber). Only lines in ``lines`` that start with ATOM/HETATM are
    touched; the k-th such line corresponds to the k-th atom of ``array``
    (biotite writes atoms in array order). Values are left-justified in
    the 4 columns, per PDB convention.

    Parameters
    ----------
    lines : list[str]
        The PDB record lines to modify in place.
    array : AtomArray
        The structure being written; the ``seg_id`` annotation must exist
        when ``lines`` contains ATOM/HETATM records.

    Returns
    -------
    None
    """
    if "seg_id" not in array.get_annotation_categories():
        return
    seg_ids = np.char.ljust(np.asarray(array.seg_id, dtype=str), 4)
    atom_i = 0
    for idx, line in enumerate(lines):
        if line.startswith(("ATOM", "HETATM")):
            lines[idx] = line[:72] + seg_ids[atom_i] + line[76:]
            atom_i += 1
