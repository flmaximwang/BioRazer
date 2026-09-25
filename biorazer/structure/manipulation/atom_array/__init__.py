"""Manipulation of biotite's ``AtomArray``: annotations, side chains, sequences.

The operand of every function here is an ``AtomArray``
(:mod:`biorazer.structure.objects.bt_atom_array`) -- residue ids and
small-molecule bonds are read or written on it, parts of it are removed,
taken from a second array, or rebuilt.  The operations on biorazer's own
internal-coordinate container live in the sibling package
:mod:`biorazer.structure.manipulation.internal_coord`.

- :mod:`~biorazer.structure.manipulation.atom_array.annotation` -- residue ids,
  small-molecule bonds.
- :mod:`~biorazer.structure.manipulation.atom_array.modification` -- side
  chains, hydrogens, backbone-only scaffolds.
- :mod:`~biorazer.structure.manipulation.atom_array.mutation` -- side chains
  rebuilt from a rotamer library, and point mutation.

The private helpers of :mod:`~biorazer.structure.manipulation.atom_array.util`
are not re-exported: they are shared by the modules above, not a public API.
"""

from biorazer.structure.manipulation.atom_array.annotation import (
    add_bonds_to_organic,
    get_renumbered_res_ids,
)
from biorazer.structure.manipulation.atom_array.modification import (
    add_hydrogens,
    mutate_without_side_chains,
    remove_side_chains,
    replace_side_chains,
)
from biorazer.structure.manipulation.atom_array.mutation import (
    build_side_chain,
    dihedral,
    mutate,
    rotamer_candidates,
    ss_from_phi_psi,
)

__all__ = [
    # annotation
    "get_renumbered_res_ids",
    "add_bonds_to_organic",
    # modification
    "add_hydrogens",
    "remove_side_chains",
    "replace_side_chains",
    "mutate_without_side_chains",
    # mutation
    "dihedral",
    "ss_from_phi_psi",
    "rotamer_candidates",
    "build_side_chain",
    "mutate",
]
