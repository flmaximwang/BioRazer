"""Tests for :mod:`biorazer.structure.bridge.sequence`."""

import numpy as np
import pytest

from biorazer.structure.bridge import AtomArray_ProteinSequence
from biorazer.structure.objects import AtomArray


def _array(residue_names, chain_ids=None, res_ids=None):
    """A minimal CA-only ``AtomArray`` carrying the given residue names."""
    size = len(residue_names)
    array = AtomArray(size)
    array.coord = np.zeros((size, 3))
    array.chain_id = np.array(chain_ids if chain_ids else ["A"] * size)
    array.res_id = np.array(res_ids if res_ids else list(range(1, size + 1)))
    array.res_name = np.array(residue_names)
    array.atom_name = np.array(["CA"] * size)
    array.element = np.array(["C"] * size)
    array.hetero = np.zeros(size, bool)
    return array


class TestAtomArrayProteinSequence:
    def test_reads_one_sequence_per_chain(self):
        array = _array(["ALA", "GLY", "SER"], chain_ids=["A", "A", "B"])
        sequences = AtomArray_ProteinSequence(input_io=array).convert()
        assert {str(chain) for chain in sequences} == {"A", "B"}
        assert str(sequences[np.str_("A")]) == "AG"
        assert str(sequences[np.str_("B")]) == "S"

    def test_reads_multiple_atoms_per_residue_once(self):
        array = _array(["ALA", "ALA", "GLY"], res_ids=[1, 1, 2])
        sequences = AtomArray_ProteinSequence(input_io=array).convert()
        assert str(sequences[np.str_("A")]) == "AG"

    def test_maps_non_standard_residues(self):
        """MSE reads as M and UNK as X, via biotite's 3-to-1 table."""
        array = _array(["ALA", "MSE", "UNK"], chain_ids=["A", "A", "B"])
        sequences = AtomArray_ProteinSequence(input_io=array).convert()
        assert str(sequences[np.str_("A")]) == "AM"
        assert str(sequences[np.str_("B")]) == "X"

    def test_rejects_a_residue_the_table_does_not_know(self):
        """A ligand/solvent residue is named instead of being skipped."""
        array = _array(["ALA", "GLY", "HOH"])
        with pytest.raises(ValueError) as error:
            AtomArray_ProteinSequence(input_io=array).convert()
        assert "'HOH'" in str(error.value)
        assert "'A'" in str(error.value)

    def test_rejects_an_object_that_is_not_an_atom_array(self):
        with pytest.raises(TypeError):
            AtomArray_ProteinSequence(input_io="not an array").convert()
