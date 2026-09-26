"""Converters between biotite's ``AtomArray`` and biotite's ``ProteinSequence``.

The structure -> sequence direction lives here rather than on the object:
per the package convention, object -> object conversions live in ``bridge``
(see :mod:`biorazer.structure.bridge.atom_array` and
:mod:`biorazer.structure.bridge.icchain`).

Scope note: ``ProteinSequence`` belongs to biotite's *sequence* layer, so it
comes from :mod:`biorazer.sequence.objects` (the sequence catalogue), not from
biotite directly -- the same "one home per object name" rule that governs
:mod:`biorazer.structure.objects`.  This module only *bridges* the two layers;
it defines no sequence object of its own.
"""

import biotite.structure as bt_struct

from biorazer.io import Converter
from biorazer.sequence.objects import ProteinSequence
from biorazer.structure.objects import AtomArray


def _residue_letter(chain_id, residue_name):
    """3-letter residue name -> 1-letter protein letter, with a clear error.

    ``ProteinSequence.convert_letter_3to1`` raises a bare ``KeyError`` for any
    name it does not know (``HOH``, ``EDO``, an unknown ligand, a nucleotide);
    the chain and residue are added here so the caller learns *what* in the
    structure was not protein.
    """
    try:
        return ProteinSequence.convert_letter_3to1(residue_name)
    except KeyError:
        raise ValueError(
            f"chain {str(chain_id)!r}: residue {str(residue_name)!r} has no "
            "3-to-1 protein letter -- a chain carrying solvent, ligands or "
            "nucleotides cannot be read as a ProteinSequence"
        ) from None


class AtomArray_ProteinSequence(Converter):
    """
    Reads the protein sequence of every chain of an ``AtomArray``.

    ``input_io`` holds a biotite ``AtomArray``; :meth:`convert` returns a
    ``{chain_id: ProteinSequence}`` dict with one entry per chain id in the
    array, in the order the chains appear in it.

    Residues are read in array order through ``biotite.structure.get_residues``
    and mapped with the standard 3-to-1 table (so ``MSE`` reads as ``M`` and
    ``UNK`` as ``X``).  A residue the table does not know raises ``ValueError``
    naming the chain and residue, instead of silently skipping it or producing
    a sequence shorter than the chain.
    """

    def convert(self, **kwargs):
        array = self.input_io
        if not isinstance(array, AtomArray):
            raise TypeError(
                f"expected an AtomArray, got {type(array).__name__}")
        sequences = {}
        for chain_id in bt_struct.get_chains(array):
            chain = array[array.chain_id == chain_id]
            _, residue_names = bt_struct.get_residues(chain)
            letters = [_residue_letter(chain_id, name) for name in residue_names]
            sequences[chain_id] = ProteinSequence("".join(letters))
        return sequences
