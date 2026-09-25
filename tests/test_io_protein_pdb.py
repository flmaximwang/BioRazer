"""Tests for :class:`biorazer.structure.io.protein.AtomArray_Pdb`.

biotite writes the ATOM/HETATM records; biorazer adds the LINK / SSBOND
records and the SEGID columns on top of that.  Both are written by editing
the serialised file, which biotite 1.7 silently broke -- its ``PDBFile`` is
now backed by the Rust implementation, where ``lines`` is a read-only copy,
so anything that post-edited it disappeared without an exception.  These
tests therefore pin the produced *text*, not the in-memory ``PDBFile``.
"""

import io

import numpy as np

from biorazer.structure.io.protein import AtomArray_Pdb
from biorazer.structure.objects import AtomArray, BondList


def _cys_cys_hem(res_id: int = 10) -> AtomArray:
    """Two CYS residues bonded SG-SG, plus a HEM bonded to the first SG."""
    array = AtomArray(11)
    array.atom_name = np.array(["N", "CA", "C", "O", "SG"] * 2 + ["FE"])
    array.res_name = np.array(["CYS"] * 5 + ["CYS"] * 5 + ["HEM"])
    array.chain_id = np.array(["A"] * 10 + ["L"])
    array.res_id = np.array([res_id] * 5 + [res_id + 1] * 5 + [1])
    array.ins_code = np.array([""] * 11)
    array.element = np.array(["N", "C", "C", "O", "S"] * 2 + ["FE"])
    array.coord = np.zeros((11, 3))
    # (SG_i, SG_j) disulfide bond;  (SG_i, FE) protein - ligand bond
    array.bonds = BondList(11, np.array([[4, 9, 1], [4, 10, 1]]))
    array.set_annotation("seg_id", np.array(["FIB1"] * 10 + ["LIG1"]))
    return array


def _write(array: AtomArray, hybrid36: bool = False) -> str:
    buffer = io.StringIO()
    AtomArray_Pdb(output_io=buffer).write(array, hybrid36=hybrid36)
    return buffer.getvalue()


def _records(text: str, prefix: str) -> list[str]:
    return [line for line in text.splitlines() if line.startswith(prefix)]


class TestAtomArrayPdbRecords:
    def test_writes_one_link_record_for_the_protein_ligand_bond(self):
        link = _records(_write(_cys_cys_hem()), "LINK")
        assert len(link) == 1
        assert "SG" in link[0] and "CYS" in link[0]      # protein side
        assert "FE" in link[0] and "HEM" in link[0]      # ligand side

    def test_writes_one_ssbond_record_for_the_disulfide(self):
        ssbond = _records(_write(_cys_cys_hem()), "SSBOND")
        assert len(ssbond) == 1
        fields = ssbond[0].split()
        assert fields[0:2] == ["SSBOND", "1"]
        # both bonded residues are the two CYS of chain A, in order
        assert [fields[i] for i in (2, 5)] == ["CYS", "CYS"]
        assert [fields[i] for i in (3, 6)] == ["A", "A"]
        assert [fields[i] for i in (4, 7)] == ["10", "11"]

    def test_writes_no_records_without_bonds(self):
        array = _cys_cys_hem()
        array.bonds = None
        text = _write(array)
        assert _records(text, "LINK") == []
        assert _records(text, "SSBOND") == []

    def test_writes_link_and_ssbond_with_hybrid36(self):
        """The combination that biotite 1.7 dropped silently."""
        for res_id in (10, 10000):
            text = _write(_cys_cys_hem(res_id), hybrid36=True)
            assert len(_records(text, "LINK")) == 1
            assert len(_records(text, "SSBOND")) == 1

    def test_writes_the_seg_id_into_the_segid_columns(self):
        """Columns 73-76 of every ATOM/HETATM record, in atom order."""
        array = _cys_cys_hem()
        text = _write(array)
        atom_lines = _records(text, "ATOM") + _records(text, "HETATM")
        assert len(atom_lines) == len(array)
        for line, seg_id in zip(atom_lines, array.seg_id):
            assert line[72:76] == seg_id
        # the injected columns must not have moved anything else
        assert all(len(line) == 80 for line in atom_lines)

    def test_ends_with_exactly_one_newline(self):
        """biotite 1.6 terminates the last record, 1.7 does not."""
        text = _write(_cys_cys_hem())
        assert text.endswith("\n") and not text.endswith("\n\n")

    def test_keeps_the_record_columns_when_hybrid36_is_on(self):
        """resSeq is right-justified in a 4-column field, as in ATOM records.

        ``encode_hybrid36`` returns an unpadded string ('10', not '  10'), so
        the LINK/SSBOND records used to be shorter than the ATOM records and
        every column after the resSeq field was shifted left.
        """
        plain = _write(_cys_cys_hem(), hybrid36=False)
        hybrid = _write(_cys_cys_hem(), hybrid36=True)
        for prefix in ("LINK", "SSBOND"):
            assert len(_records(hybrid, prefix)[0]) == len(_records(plain, prefix)[0])
        # protein side of LINK, then the ligand side
        assert _records(hybrid, "LINK")[0][22:26] == "  10"
        assert _records(hybrid, "LINK")[0][52:56] == "   1"
        # both cysteines of the SSBOND
        assert _records(hybrid, "SSBOND")[0][17:21] == "  10"
        assert _records(hybrid, "SSBOND")[0][31:35] == "  11"

    def test_keeps_the_record_columns_in_the_hybrid_letter_range(self):
        """Above 9999 hybrid-36 goes to letters: 'A000', 'A001', ..."""
        hybrid = _write(_cys_cys_hem(res_id=10000), hybrid36=True)
        assert _records(hybrid, "LINK")[0][22:26] == "A000"
        assert _records(hybrid, "LINK")[0][52:56] == "   1"
        assert _records(hybrid, "SSBOND")[0][17:21] == "A000"
        assert _records(hybrid, "SSBOND")[0][31:35] == "A001"
