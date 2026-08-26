"""Converters between structure file formats and in-memory representations.

This module is a package: the public :class:`Converter` subclasses live here
(``__init__.py``), while their helper functions are split across sibling
modules by concern:

- :mod:`._pdb_records` — PDB LINK/SSBOND/SEGID record writers (used by
  :class:`AtomArray_Pdb`).
- :mod:`._io` — generic ``str``/``Path``/``io.StringIO`` target helpers.
- :mod:`._internal_coords` — biopython ``Bio.PDB.internal_coords`` helpers.
- :mod:`._pose` — PyRosetta Pose helpers (imports stay lazy: PyRosetta is an
  optional dependency).
"""

from biotite.structure.io import pdb, pdbx
from biotite.structure import AtomArray
import biotite.structure as bio_struc
from biorazer.database.amino_acid import AMINO_ACIDS_3TO1_UPPER
from biorazer.io import Converter

from ._pdb_records import _format_link_records, _format_ssbond_records, _inject_seg_ids
from ._io import _io_target, _written_text
from ._internal_coords import (
    MMCIFIO,
    MMCIFParser,
    PDBIO,
    PDBParser,
    InternalCoord_SMCRA,
    SMCRA_InternalCoord,
)
from ._pose import _dump_pose, _pose_from_io


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
# biopython (Bio.PDB.internal_coords) is a CORE dependency of biorazer, so the
# internal-coordinate converters are always available. PyRosetta, however, is
# NOT a core dependency (it is not on PyPI), so all PyRosetta imports in
# ``._pose`` are LAZY: importing this package never fails when pyrosetta is
# missing -- the ImportError is raised only when a Pose converter is actually
# used. Declare pyrosetta via the ``pyrosetta`` extra in pyproject.toml.
# ---------------------------------------------------------------------------


class Pdb_InternalCoord(Converter):
    """
    Converts a PDB file to biopython's internal-coordinate representation.

    Reads the file with ``Bio.PDB.PDBParser`` and returns a list of
    :class:`Bio.PDB.internal_coords.IC_Chain` (one per chain of the first
    model), each with its Cartesian coordinates already converted to
    internal/torsion coordinates via ``atom_to_internal_coordinates``.
    """

    def read(self, **kwargs):
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("biorazer", self.input_io)
        return SMCRA_InternalCoord(input_io=structure).convert()


class Cif_InternalCoord(Converter):
    """
    Converts an mmCIF file to biopython's internal-coordinate representation.

    Reads the file with ``Bio.PDB.MMCIFParser`` and returns a list of
    :class:`Bio.PDB.internal_coords.IC_Chain` (one per chain of the first
    model), each with its Cartesian coordinates already converted to
    internal/torsion coordinates.
    """

    def read(self, **kwargs):
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("biorazer", self.input_io)
        return SMCRA_InternalCoord(input_io=structure).convert()


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
        structure = InternalCoord_SMCRA(input_io=tmp).convert()
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
        structure = InternalCoord_SMCRA(input_io=tmp).convert()
        writer = MMCIFIO()
        writer.set_structure(structure)
        writer.save(_io_target(self.output_io))
        return _written_text(self.output_io)


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
