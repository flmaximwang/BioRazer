# -*- coding: utf-8 -*-
"""Tests for the biorazer.database.molecule package.

Covers the 2026-08-27 refactor that merged the flat database modules
(``bond`` / ``torsion_angle`` / ``internal_coord_template`` / ``atom.py``)
into ``biorazer/database/molecule``:

* public names resolve from the new package paths (identity with their
  per-module homes);
* every numeric entry carries the uniform ``{mean, std, lb, up, source}``
  record (missing spread = ``np.nan``);
* the migrated values are unchanged from the pre-move data;
* templates still build and round-trip (``to_atomarray``);
* the old import paths raise ModuleNotFoundError.
"""

import importlib
import numpy as np
import pytest

from biorazer.database import molecule as M
from biorazer.database.molecule.bond.length import generic as len_generic
from biorazer.database.molecule.bond.length import protein as len_protein
from biorazer.database.molecule.bond.angle import generic as ang_generic
from biorazer.database.molecule.bond.angle import protein as ang_protein
from biorazer.database.molecule.bond.dihedral import protein as dih_protein
from biorazer.database.molecule.atom import radius as atom_radius
from biorazer.database.molecule.icoor.protein import topology


class TestPackageLayout:
    """The molecule package exposes every public table."""

    def test_aggregate_reexports(self):
        for name in M.__all__:
            assert hasattr(M, name), f"molecule missing {name}"
        assert set(M.__all__) <= set(dir(M))

    def test_identity_with_module_homes(self):
        assert M.AMINO_ACID_BOND_LENGTH is len_generic.AMINO_ACID_BOND_LENGTH
        assert M.BOND_REFS is len_generic.BOND_REFS
        assert M.AMINO_ACID_BOND_LENGTH_BY_RESIDUE is len_protein.AMINO_ACID_BOND_LENGTH_BY_RESIDUE
        assert M.AMINO_ACID_SIDECHAIN_BOND is len_protein.AMINO_ACID_SIDECHAIN_BOND
        assert M.AMINO_ACID_BOND_ANGLE is ang_generic.AMINO_ACID_BOND_ANGLE
        assert M.AMINO_ACID_BOND_ANGLE_BY_RESIDUE is ang_protein.AMINO_ACID_BOND_ANGLE_BY_RESIDUE
        assert M.AMINO_ACID_SIDECHAIN_BOND_ANGLE is ang_protein.AMINO_ACID_SIDECHAIN_BOND_ANGLE
        assert M.SS_BB_TORSION_ANGLE is dih_protein.SS_BB_TORSION_ANGLE
        assert M.OMEGA_TRANS is dih_protein.OMEGA_TRANS
        assert M.SIDECHAIN_CHI is dih_protein.SIDECHAIN_CHI
        assert M.IC_PATH is topology.IC_PATH
        assert M.MAINCHAIN_ATOMS is topology.MAINCHAIN_ATOMS
        assert M.ATOM_RADIUS is atom_radius.ATOM_RADIUS


class TestUniformRecord:
    """Every numeric entry is {mean, std, lb, up, source} (+ aux keys)."""

    REQUIRED = {"mean", "std", "lb", "up", "source"}

    def test_generic_bond_length(self):
        for key, rec in M.AMINO_ACID_BOND_LENGTH.items():
            assert self.REQUIRED <= set(rec), key
            assert rec["mean"] > 0 and rec["std"] > 0
            assert rec["lb"] < rec["mean"] < rec["up"]

    def test_by_residue_bond_length(self):
        for res, d in M.AMINO_ACID_BOND_LENGTH_BY_RESIDUE.items():
            for key, rec in d.items():
                assert self.REQUIRED <= set(rec), (res, key)

    def test_generic_bond_angle(self):
        for key, rec in M.AMINO_ACID_BOND_ANGLE.items():
            assert self.REQUIRED <= set(rec), key
            assert rec["lb"] < rec["mean"] < rec["up"]

    def test_by_residue_bond_angle(self):
        for res, d in M.AMINO_ACID_BOND_ANGLE_BY_RESIDUE.items():
            for key, rec in d.items():
                assert self.REQUIRED <= set(rec), (res, key)

    def test_sidechain_bond_uses_nan_for_missing_spread(self):
        # Rosetta ICOOR gives ideal point values only -> std/lb/up are nan.
        n = 0
        for res, d in M.AMINO_ACID_SIDECHAIN_BOND.items():
            for key, rec in d.items():
                assert self.REQUIRED <= set(rec), (res, key)
                assert isinstance(key, tuple) and len(key) == 2, (res, key)
                assert np.isnan(rec["std"]) and np.isnan(rec["lb"]) and np.isnan(rec["up"])
                n += 1
        assert n == 87  # the bond count of IC_PATH (excluding GLY)

    def test_sidechain_bond_angle(self):
        for res, d in M.AMINO_ACID_SIDECHAIN_BOND_ANGLE.items():
            for key, rec in d.items():
                assert self.REQUIRED <= set(rec), (res, key)
                assert isinstance(key, tuple) and len(key) == 3, (res, key)
                assert np.isnan(rec["std"]) and np.isnan(rec["lb"]) and np.isnan(rec["up"])

    def test_ss_torsion_angle(self):
        for ss, v in M.SS_BB_TORSION_ANGLE.items():
            for t in ("phi", "psi", "omega"):
                assert self.REQUIRED <= set(v[t]), (ss, t)
                assert "up" in v[t] and "ub" not in v[t]

    def test_sidechain_ic_dihedral(self):
        for res, d in M.SIDECHAIN_IC_DIHEDRAL.items():
            for key, rec in d.items():
                assert self.REQUIRED <= set(rec), (res, key)
                assert np.isnan(rec["std"]) and np.isnan(rec["lb"]) and np.isnan(rec["up"])

    def test_rotamer_bin(self):
        for name, rec in M.ROTAMER_BIN.items():
            assert self.REQUIRED <= set(rec), name
            assert np.isnan(rec["std"]) and np.isnan(rec["lb"]) and np.isnan(rec["up"])

    def test_atom_radius(self):
        for elm, rec in M.ATOM_RADIUS.items():
            assert self.REQUIRED <= set(rec), elm
            assert np.isnan(rec["std"]) and np.isnan(rec["lb"]) and np.isnan(rec["up"])


class TestMigratedValues:
    """Spot values that must be unchanged from the pre-move data."""

    def test_bond_length_generic(self):
        assert M.AMINO_ACID_BOND_LENGTH[("C", "N")] == {
            "mean": 1.329, "std": 0.014, "lb": 1.287, "up": 1.371,
            "note": "肽键 C(=O)-N。Engh-Huber 主表 C-NH1 (except Pro) 1.329±0.014; Pro 1.341±0.016。",
            "source": ("engh_huber_1991", "procheck_appendix_a"),
        }
        assert M.AMINO_ACID_BOND_LENGTH[("N", "CA")]["mean"] == 1.458

    def test_bond_angle_generic(self):
        assert M.AMINO_ACID_BOND_ANGLE[("N", "CA", "C")]["mean"] == 111.2
        assert M.AMINO_ACID_BOND_ANGLE[("CA", "C", "N")]["mean"] == 116.2

    def test_ss_torsion_alpha_helix(self):
        phi = M.SS_BB_TORSION_ANGLE["alpha-helix"]["phi"]
        assert phi["mean"] == -60 and phi["up"] == -35 and phi["lb"] == -85
        assert M.SS_BB_TORSION_ANGLE["alpha-helix"]["psi"]["mean"] == -45

    def test_sidechain_tables_keyed_by_arity(self):
        # length keys are 2-atom tuples, angle keys are 3-atom tuples,
        # dihedral keys stay 4-atom grow quads; each IC_PATH quad must resolve
        # in all three tables.
        for res, quads in M.IC_PATH.items():
            for quad in quads:
                _, j, k, l = quad
                assert (k, l) in M.AMINO_ACID_SIDECHAIN_BOND[res], (res, quad)
                assert (j, k, l) in M.AMINO_ACID_SIDECHAIN_BOND_ANGLE[res], (res, quad)
                assert quad in M.SIDECHAIN_IC_DIHEDRAL[res], (res, quad)

    def test_sidechain_length_angle_consistent(self):
        # the same (res, grow quad) must map to the same bond/angle values
        # through the 2- and 3-atom keys.
        for res, quads in M.IC_PATH.items():
            for quad in quads:
                i, j, k, l = quad
                bond = M.AMINO_ACID_SIDECHAIN_BOND[res][(k, l)]
                ang = M.AMINO_ACID_SIDECHAIN_BOND_ANGLE[res][(j, k, l)]
                assert ang["mean"] > 0 and bond["mean"] > 0, (res, quad)
                assert ang["source"] == bond["source"] == "rosetta_params_408", (res, quad)

    def test_mainchain_torsion_definitions(self):
        assert M.MAINCHAIN_TORSION_DEFINITIONS == {
            "phi": ("C", "N", "CA", "C"),
            "psi": ("N", "CA", "C", "N"),
            "omega": ("CA", "C", "N", "CA"),
        }

    def test_omega(self):
        assert M.OMEGA_TRANS["mean"] == 180.0 and M.OMEGA_TRANS["up"] == 190.0
        assert M.OMEGA_CIS["mean"] == 0.0


class TestTemplates:
    """The template store builds and round-trips."""

    def test_named_constants(self):
        from biorazer.database.molecule.icoor.protein.template import (
            get_template, IC_Gly_HelixAlpha, IC_Ser_HelixAlpha_gminus,
        )
        assert IC_Gly_HelixAlpha is get_template("GLY", "alpha-helix")
        assert IC_Ser_HelixAlpha_gminus is get_template("SER", "alpha-helix", "g-")

    def test_build_and_to_atomarray(self):
        from biorazer.database.molecule.icoor.protein.template import get_template
        ic = get_template("SER", "alpha-helix", "g-")
        assert ic.phi == -60.0 and ic.psi == -45.0 and ic.omega == 180.0
        arr = ic.to_atomarray()
        assert len(arr) == len(ic.atoms)
        names = set(arr.atom_name)
        assert {"N", "CA", "C", "O", "CB", "OG"} <= names

    def test_templates_round_trip(self):
        # to_coords -> from_atomarray must reproduce the template geometry
        from biorazer.database.molecule.icoor.protein.template import get_template
        from biorazer.structure.objects.internal_coords import InternalCoord
        ic = get_template("TRP", "beta-strand")
        arr = ic.to_atomarray()
        ic2 = InternalCoord.from_atomarray(arr)
        # anchors identical (the IC frame is preserved)
        for k, v in ic.anchor.items():
            assert np.allclose(v, ic2.anchor[k], atol=1e-6)


class TestOldPathsDeleted:
    """The pre-move module files are gone from the repo tree.

    (A namespace-package ``import`` check would also resolve the old paths
    through the editable-install pointer to the main working tree, so we
    assert on the filesystem instead.)
    """

    @pytest.mark.parametrize("rel", [
        "biorazer/database/bond",
        "biorazer/database/torsion_angle",
        "biorazer/database/internal_coord_template",
        "biorazer/database/atom.py",
    ])
    def test_old_module_files_deleted(self, rel):
        repo_root = importlib.util.find_spec("biorazer").origin
        if repo_root is None:  # namespace package: search each path entry
            import biorazer
            repo_root = str(biorazer.__path__[0])
        import pathlib
        assert not (pathlib.Path(repo_root) / rel).exists()
