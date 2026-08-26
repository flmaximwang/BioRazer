import io
import os
import tempfile
from pathlib import Path
import numpy as np
from biotite.structure.io import pdb, pdbx
from biotite.structure import AtomArray
import biotite.structure as bio_struc
import biotite.sequence as bio_seq
from biorazer.database.amino_acid import AMINO_ACIDS_3TO1_UPPER
from biotite.structure.io.pdb.hybrid36 import encode_hybrid36
from biorazer.io import Converter
from biorazer.sequence.io import StrDict_Fasta


class Pdb_AtomArray(Converter):
    def read(self, **kwargs) -> AtomArray:
        return pdb.get_structure(pdb.PDBFile.read(self.input_io), **kwargs)[0]


class Cif_AtomArray(Converter):
    def read(self, **kwargs) -> AtomArray:
        return pdbx.get_structure(pdbx.CIFFile.read(self.input_io), **kwargs)[0]


class AtomArray_Cif(Converter):
    def write(self, tmp, **kwargs):
        output_file_obj = pdbx.CIFFile()
        pdbx.set_structure(output_file_obj, tmp, **kwargs)
        output_file_obj.write(self.output_io)


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


class AtomArray_Pdb(Converter):
    def write(self, tmp: AtomArray, hybrid36: bool = False):
        """
        Write a PDB file. Intermolecular covalent bonds between a protein
        and a small molecule are additionally written as LINK records, and
        disulfide bonds as SSBOND records, appended at the end of the file.
        If the array carries a ``seg_id`` annotation, it is written into
        the SEGID columns (73-76) of the ATOM/HETATM records.
        """
        output_file_obj = pdb.PDBFile()
        pdb.set_structure(output_file_obj, tmp, hybrid36=hybrid36)
        if tmp.bonds is not None and tmp.coord.ndim == 2:
            output_file_obj.lines.extend(
                _format_link_records(tmp, hybrid36)
                + _format_ssbond_records(tmp, hybrid36)
            )
        _inject_seg_ids(output_file_obj.lines, tmp)
        output_file_obj.write(self.output_io)


class Pdb_StrDict(Converter):
    """
    Converts a PDB file to a sequence dictionary.
    """

    def read(self, remove_gaps=False, **kwargs) -> dict:
        structure = Pdb_AtomArray(self.input_io, self.output_io).read(**kwargs)
        chain_ids = bio_struc.get_chains(structure)
        res = {}
        for chain_id in chain_ids:
            chain_structure = structure[structure.chain_id == chain_id]
            res_ids, res_names = bio_struc.get_residues(chain_structure)
            res_ids = list(res_ids)
            res_names = list(res_names)
            one_char_res_names = []
            for i in range(min(res_ids), max(res_ids) + 1):
                try:
                    idx = res_ids.index(i)
                except ValueError:
                    # Missing residue
                    one_char_res_names.append("-") # - means gaps
                    continue
                res_name = res_names[idx]
                if len(res_name) != 3:
                    # Nucleotides
                    break
                try:
                    one_char_res_names.append(AMINO_ACIDS_3TO1_UPPER[res_name])
                except KeyError:
                    # Non-standard amino acid or ligand
                    one_char_res_names.append("X") # X means unrecognized residues

            sequence = "".join(one_char_res_names)
            if remove_gaps:
                res[chain_id] = sequence.replace("-", "")

        return res


# ---------------------------------------------------------------------------
# Optional-dependency converters: biopython internal coords / PyRosetta Pose
#
# These converters build representations that require extra packages.
# biopython (Bio.PDB.internal_coords) is a CORE dependency of biorazer, so the
# internal-coordinate converters are always available (the lazy import below is
# just defensive). PyRosetta, however, is NOT a core dependency (it is not on
# PyPI), so all PyRosetta imports are LAZY: importing this module never fails
# when pyrosetta is missing -- the ImportError is raised only when a Pose
# converter is actually used. Declare pyrosetta via the ``pyrosetta`` extra in
# pyproject.toml.
# ---------------------------------------------------------------------------


def _import_internal_coords():
    """Import the Bio.PDB machinery used by the internal-coordinate converters."""
    try:
        from Bio.PDB import MMCIFIO, MMCIFParser, PDBIO, PDBParser
        from Bio.PDB import Model as BioModel, Structure as BioStructure
        from Bio.PDB.internal_coords import IC_Chain
    except ImportError as e:  # pragma: no cover - biopython is a core dependency
        raise ImportError(
            "Biopython is required for the internal-coordinate converters "
            "(Pdb_InternalCoord / Cif_InternalCoord / InternalCoord_Pdb / "
            "InternalCoord_Cif); it is a core dependency of biorazer."
        ) from e
    return PDBParser, MMCIFParser, PDBIO, MMCIFIO, BioStructure, BioModel, IC_Chain


def _ic_chains(structure, IC_Chain):
    """Build one ``IC_Chain`` per chain (first model), converting XYZ -> internal coords."""
    chains = []
    model = next(iter(structure))
    for chain in model:
        ic = IC_Chain(chain)
        ic.atom_to_internal_coordinates()
        chains.append(ic)
    return chains


def _rebuild_structure(ic_chains, BioStructure, BioModel):
    """
    Regenerate Cartesian coords from a chain's internal coords and reassemble
    a ``Bio.PDB.Structure`` holding all the (rebuilt) chains.
    """
    if not isinstance(ic_chains, (list, tuple)):
        ic_chains = [ic_chains]
    for ic in ic_chains:
        ic.internal_to_atom_coordinates()
    structure = BioStructure.Structure("biorazer")
    model = BioModel.Model(0)
    for ic in ic_chains:
        chain = ic.chain
        if chain.parent is not None:
            chain.parent.detach_child(chain.id)
        model.add(chain)
    structure.add(model)
    return structure


def _written_text(output_io):
    """Return the written text when the target is an ``io.StringIO``, else None."""
    if isinstance(output_io, io.StringIO):
        return output_io.getvalue()
    return None


def _io_target(output_io):
    """Bio.PDB PDBIO/MMCIFIO accept a str filename or a file object (not a Path)."""
    return str(output_io) if isinstance(output_io, Path) else output_io


class Pdb_InternalCoord(Converter):
    """
    Converts a PDB file to biopython's internal-coordinate representation.

    Reads the file with ``Bio.PDB.PDBParser`` and returns a list of
    :class:`Bio.PDB.internal_coords.IC_Chain` (one per chain of the first
    model), each with its Cartesian coordinates already converted to
    internal/torsion coordinates via ``atom_to_internal_coordinates``.
    """

    def read(self, **kwargs):
        PDBParser, _, _, _, _, _, IC_Chain = _import_internal_coords()
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("biorazer", self.input_io)
        return _ic_chains(structure, IC_Chain)


class Cif_InternalCoord(Converter):
    """
    Converts an mmCIF file to biopython's internal-coordinate representation.

    Reads the file with ``Bio.PDB.MMCIFParser`` and returns a list of
    :class:`Bio.PDB.internal_coords.IC_Chain` (one per chain of the first
    model), each with its Cartesian coordinates already converted to
    internal/torsion coordinates.
    """

    def read(self, **kwargs):
        _, MMCIFParser, _, _, _, _, IC_Chain = _import_internal_coords()
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("biorazer", self.input_io)
        return _ic_chains(structure, IC_Chain)


class InternalCoord_Pdb(Converter):
    """
    Regenerates a PDB file from biopython's internal-coordinate representation.

    ``tmp`` is an :class:`Bio.PDB.internal_coords.IC_Chain` or a list of them
    (as returned by :class:`Pdb_InternalCoord` / :class:`Cif_InternalCoord`).
    Each chain is rebuilt back to Cartesian coordinates from its
    internal/torsion coordinates via ``internal_to_atom_coordinates`` and the
    chains are written as a single PDB file.
    """

    def write(self, tmp, **kwargs):
        _, _, PDBIO, _, BioStructure, BioModel, _ = _import_internal_coords()
        structure = _rebuild_structure(tmp, BioStructure, BioModel)
        writer = PDBIO()
        writer.set_structure(structure)
        writer.save(_io_target(self.output_io))
        return _written_text(self.output_io)


class InternalCoord_Cif(Converter):
    """
    Regenerates an mmCIF file from biopython's internal-coordinate
    representation (the inverse of :class:`Cif_InternalCoord`).
    """

    def write(self, tmp, **kwargs):
        _, _, _, MMCIFIO, BioStructure, BioModel, _ = _import_internal_coords()
        structure = _rebuild_structure(tmp, BioStructure, BioModel)
        writer = MMCIFIO()
        writer.set_structure(structure)
        writer.save(_io_target(self.output_io))
        return _written_text(self.output_io)


_PYROSETTA_INITIALIZED = False


def _import_pyrosetta():
    """Lazily import PyRosetta, raising an informative ImportError if missing."""
    try:
        import pyrosetta
    except ImportError as e:
        raise ImportError(
            "PyRosetta is required for the Pose converters (Pdb_Pose / Cif_Pose / "
            "Pose_Pdb / Pose_Cif). PyRosetta is not on PyPI; install it, e.g. via "
            "conda (-c https://conda.rosettacommons.org -c conda-forge pyrosetta) or "
            "from a pyrosetta wheel, then `pip install .[pyrosetta]`."
        ) from e
    return pyrosetta


def _ensure_pyrosetta_init():
    """Call ``pyrosetta.init`` once per process (it is process-global state)."""
    global _PYROSETTA_INITIALIZED
    if _PYROSETTA_INITIALIZED:
        return
    _import_pyrosetta().init("-mute all", silent=True)
    _PYROSETTA_INITIALIZED = True


def _pose_from_io(input_io, suffix):
    """Read a PyRosetta Pose from a str/Path (auto-detected by extension) or StringIO."""
    pyrosetta = _import_pyrosetta()
    _ensure_pyrosetta_init()
    if isinstance(input_io, (str, Path)):
        return pyrosetta.pose_from_file(str(input_io))
    # io.StringIO: PyRosetta needs a real file, so stage the text with a
    # format suffix so pose_from_file can auto-detect the format.
    with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as fh:
        fh.write(input_io.getvalue().encode("utf-8"))
        name = fh.name
    try:
        return pyrosetta.pose_from_file(name)
    finally:
        os.remove(name)


def _dump_pose(pose, output_io, suffix):
    """Write a PyRosetta Pose to a str/Path or StringIO, in PDB or mmCIF format."""
    if isinstance(output_io, (str, Path)):
        if suffix == ".pdb":
            pose.dump_pdb(str(output_io))
        else:
            pose.dump_cif(str(output_io))
        return None
    # io.StringIO: dump to a temp file, then read the text back.
    with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as fh:
        name = fh.name
    try:
        if suffix == ".pdb":
            pose.dump_pdb(name)
        else:
            pose.dump_cif(name)
        with open(name, "r") as fh:
            output_io.write(fh.read())
    finally:
        os.remove(name)
    return output_io.getvalue()


class Pdb_Pose(Converter):
    """
    Converts a PDB file to a PyRosetta :class:`rosetta.core.pose.Pose`.

    The input may be a file path (format detected by extension) or an
    ``io.StringIO``. ``pyrosetta.init`` is called automatically (once per
    process) if it has not been already.
    """

    def read(self, **kwargs):
        return _pose_from_io(self.input_io, ".pdb")


class Cif_Pose(Converter):
    """
    Converts an mmCIF file to a PyRosetta :class:`rosetta.core.pose.Pose`.

    The input may be a file path (format detected by extension) or an
    ``io.StringIO``. ``pyrosetta.init`` is called automatically (once per
    process) if it has not been already.
    """

    def read(self, **kwargs):
        return _pose_from_io(self.input_io, ".cif")


class Pose_Pdb(Converter):
    """
    Writes a PyRosetta :class:`rosetta.core.pose.Pose` as a PDB file
    (the inverse of :class:`Pdb_Pose`).
    """

    def write(self, tmp, **kwargs):
        return _dump_pose(tmp, self.output_io, ".pdb")


class Pose_Cif(Converter):
    """
    Writes a PyRosetta :class:`rosetta.core.pose.Pose` as an mmCIF file
    (the inverse of :class:`Cif_Pose`).
    """

    def write(self, tmp, **kwargs):
        return _dump_pose(tmp, self.output_io, ".cif")

