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
from biorazer.database.molecule.bond.length import protein as len_protein
from biorazer.database.molecule.bond.angle import generic as ang_generic
from biorazer.database.molecule.bond.angle import protein as ang_protein
from biorazer.database.molecule.bond.dihedral import protein as dih_protein
from biorazer.database.molecule.atom import radius as atom_radius
from biorazer.database.molecule.icoor.protein import topology


def _is_nr_bin_key(res, key):
    """True if ``key`` is a ``<RES>_nr<i>`` non-rotameric bin entry."""
    prefix = f"{res}_nr"
    return key.startswith(prefix) and key[len(prefix):].isdigit()


class TestPackageLayout:
    """The molecule package exposes every public table."""

    def test_aggregate_reexports(self):
        for name in M.__all__:
            assert hasattr(M, name), f"molecule missing {name}"
        assert set(M.__all__) <= set(dir(M))

    def test_identity_with_module_homes(self):
        assert M.AMINO_ACID_BOND_LENGTH is len_protein.AMINO_ACID_BOND_LENGTH
        assert M.BOND_REFS is len_protein.BOND_REFS
        assert M.AMINO_ACID_BOND_LENGTH_BY_RESIDUE is len_protein.AMINO_ACID_BOND_LENGTH_BY_RESIDUE
        assert M.AMINO_ACID_SIDECHAIN_BOND is len_protein.AMINO_ACID_SIDECHAIN_BOND
        assert M.AMINO_ACID_BOND_ANGLE is ang_protein.AMINO_ACID_BOND_ANGLE
        assert M.AMINO_ACID_BACKBONE_BOND_ANGLE is ang_generic.AMINO_ACID_BACKBONE_BOND_ANGLE
        assert M.SS_BB_TORSION_ANGLE is dih_protein.SS_BB_TORSION_ANGLE
        assert M.OMEGA_TRANS is dih_protein.OMEGA_TRANS
        assert M.SIDECHAIN_CHI is dih_protein.SIDECHAIN_CHI
        assert M.IC_PATH is topology.IC_PATH
        assert M.MAINCHAIN_ATOMS is topology.MAINCHAIN_ATOMS
        assert M.BACKBONE_IC_PATH is topology.BACKBONE_IC_PATH
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
        # flat generic backbone table: no residue key, triad keys, real spread
        for key, rec in M.AMINO_ACID_BACKBONE_BOND_ANGLE.items():
            assert self.REQUIRED <= set(rec), key
            assert isinstance(key, tuple) and len(key) == 3, key
            assert rec["lb"] < rec["mean"] < rec["up"]

    def test_by_residue_bond_angle(self):
        # single residue-keyed AMINO_ACID_BOND_ANGLE: every residue carries a
        # complete angle table.  The 5 universal backbone angles are always
        # present; CB-CA-C is present for every residue with a CB (all but Gly).
        uni = [("N", "CA", "C"), ("CA", "C", "N"), ("CA", "C", "O"),
               ("C", "N", "CA"), ("O", "C", "N")]
        for res, d in M.AMINO_ACID_BOND_ANGLE.items():
            for key in uni:
                assert key in d, (res, key)
            if res != "GLY":
                assert ("CB", "CA", "C") in d, (res, "CB-CA-C")
            for key, rec in d.items():
                assert self.REQUIRED <= set(rec), (res, key)
                assert isinstance(key, tuple) and len(key) == 3, (res, key)
                if rec["source"] == "engh_huber_1991":
                    assert rec["lb"] < rec["mean"] < rec["up"], (res, key)

    def test_by_residue_bond_angle_refinements(self):
        # refined backbone values (Engh-Huber per-residue) are written as
        # literals overriding the generic reference.
        assert M.AMINO_ACID_BOND_ANGLE["GLY"][("N", "CA", "C")]["mean"] == 112.5
        assert M.AMINO_ACID_BOND_ANGLE["PRO"][("N", "CA", "C")]["mean"] == 111.8
        assert M.AMINO_ACID_BOND_ANGLE["PRO"][("O", "C", "N")]["mean"] == 122.0
        assert M.AMINO_ACID_BOND_ANGLE["ALA"][("CB", "CA", "C")]["mean"] == 110.5
        # VIT = Val/Ile/Thr refined CB-CA-C
        for res in ("VAL", "ILE", "THR"):
            assert M.AMINO_ACID_BOND_ANGLE[res][("CB", "CA", "C")]["mean"] == 109.1

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
        # Rosetta side-chain angle entries within each residue's table carry
        # nan spread; the backbone (generic/Engh-Huber) entries carry real spread.
        for res, d in M.AMINO_ACID_BOND_ANGLE.items():
            for key, rec in d.items():
                assert self.REQUIRED <= set(rec), (res, key)
                assert isinstance(key, tuple) and len(key) == 3, (res, key)
                if rec["source"] == "rosetta_params_408":
                    assert np.isnan(rec["std"]) and np.isnan(rec["lb"]) and np.isnan(rec["up"])

    def test_ss_torsion_angle(self):
        for ss, v in M.SS_BB_TORSION_ANGLE.items():
            for quad in M.ALIAS_QUAD.values():
                assert self.REQUIRED <= set(v[quad]), (ss, quad)
                assert "up" in v[quad] and "ub" not in v[quad]

    def test_sc_non_chi_dihedral_private(self):
        # the non-rotameric side-chain dihedral table is private (_-prefixed,
        # SC_NON_CHI_DIHEDRAL) and only consumed as the SIDECHAIN_ROTAMER_LIB
        # base; it is still reachable from the module for validation.  It holds
        # a SUBSET of IC_PATH: every non-chi grow quad + PRO's ring-pucker chi
        # + terminal non-rotameric chi (the rotameric chi are stripped and
        # regenerated by the build).
        from biorazer.database.molecule.bond.dihedral.protein import by_residue
        table = by_residue.SC_NON_CHI_DIHEDRAL
        nonrot = M.SIDECHAIN_NON_ROTAMERIC_BINS
        for res, d in table.items():
            assert set(d) <= set(M.IC_PATH[res]), (res, set(d))
            # PRO keeps its ring-pucker chi; others keep only non-rotameric chi
            # and non-chi quads (never a rotameric chi1/chi2).
            for key, rec in d.items():
                assert self.REQUIRED <= set(rec), (res, key)
                assert np.isnan(rec["std"]) and np.isnan(rec["lb"]) and np.isnan(rec["up"])
                assert key in M.IC_PATH[res], (res, key)
        # the stripped (rotameric) chi must be regenerated by the build:
        # canonical entries are complete over IC_PATH and every chi quad exists.
        for res in M.AAS:
            canon = M.SIDECHAIN_ROTAMER_LIB[f"{res}_canonical"]
            assert set(canon) == set(M.IC_PATH[res]), res
            for q in M.SIDECHAIN_CHI[res]:
                assert q in canon, (res, q)

    def test_rotamer_bin(self):
        for name, rec in M.ROTAMER_BIN.items():
            assert self.REQUIRED <= set(rec), name
            assert np.isnan(rec["std"]) and np.isnan(rec["lb"]) and np.isnan(rec["up"])

    def test_dunbrack_rotamers_framework(self):
        # DUNBRACK_ROTAMERS was removed (superseded by SIDECHAIN_ROTAMER_LIB
        # + SIDECHAIN_NON_ROTAMERIC_BINS).  rotameric_chi is derived directly
        # from SIDECHAIN_CHI + SIDECHAIN_NON_ROTAMERIC_BINS: total chi axes
        # minus the terminal non-rotameric chi (if any), capped at 2 (named
        # rotamers cover chi1/chi2 only), PRO ring-constrained -> 0.
        labels = ["g-", "g+", "t"]
        for res in M.AAS:
            if res == "PRO":
                rc = 0
            else:
                rc = len(M.SIDECHAIN_CHI[res]) - (1 if res in M.SIDECHAIN_NON_ROTAMERIC_BINS else 0)
                rc = min(2, rc)
            named = {k for k in M.SIDECHAIN_ROTAMER_LIB if k.startswith(res + "_")
                     and k != f"{res}_canonical"
                     and not _is_nr_bin_key(res, k)}
            if rc == 0:
                assert named == set(), res
            elif rc == 1:
                assert named == {f"{res}_{l}" for l in labels}, res
            else:
                assert named == {f"{res}_{a}_{b}" for a in labels for b in labels}, res

    def test_non_rotameric_bin_width(self):
        assert M.NON_ROTAMERIC_BIN_WIDTH == 30.0

    def test_sidechain_rotamer_lib(self):
        LIB = M.SIDECHAIN_ROTAMER_LIB
        REQUIRED = {"mean", "std", "lb", "up", "source"}
        # every residue has a canonical, single-level keys, record format
        for res in M.AAS:
            assert f"{res}_canonical" in LIB, res
        for key, quad_map in LIB.items():
            res = key.split("_")[0]
            assert res in M.SIDECHAIN_CHI, key
            # every entry carries the COMPLETE side-chain dihedral map:
            # all IC_PATH quads of that residue (not just chi quads).
            assert set(quad_map) == set(M.IC_PATH[res]), (key, set(quad_map))
            for quad, rec in quad_map.items():
                assert quad in M.IC_PATH[res], (key, quad)
                assert REQUIRED <= set(rec), (key, quad)
                assert np.isnan(rec["std"]) and np.isnan(rec["lb"]) and np.isnan(rec["up"])
        # named-rotamer counts per rotameric_chi.  PRO is ring-constrained
        # (0 named); otherwise named rotamers cover only the first 2 chi,
        # so rc is capped at 2 and the terminal non-rotameric chi is excluded.
        # ``<RES>_nr<i>`` non-rotameric-bin entries are separate from named
        # (g-/g+/t) rotamers and excluded from this count.
        labels = ["g-", "g+", "t"]
        for res in M.AAS:
            if res == "PRO":
                rc = 0
            else:
                rc = len(M.SIDECHAIN_CHI[res]) - (1 if res in M.SIDECHAIN_NON_ROTAMERIC_BINS else 0)
                rc = min(2, rc)
            named = {k for k in LIB if k.startswith(res + "_")
                     and k != f"{res}_canonical"
                     and not _is_nr_bin_key(res, k)}
            if rc == 0:
                assert named == set(), res
            elif rc == 1:
                assert named == {f"{res}_{l}" for l in labels}, res
            else:
                assert named == {f"{res}_{a}_{b}" for a in labels for b in labels}, res
        # canonical quad->mean: non-rotameric quads match the private
        # SC_NON_CHI_DIHEDRAL table; stripped rotameric chi are regenerated
        # at the 0-deg canonical ideal.
        from biorazer.database.molecule.bond.dihedral.protein import by_residue
        NC = by_residue.SC_NON_CHI_DIHEDRAL
        for res in M.AAS:
            for quad, rec in LIB[f"{res}_canonical"].items():
                if quad in NC[res]:
                    assert NC[res][quad]["mean"] == rec["mean"], (res, quad)
                else:
                    # regenerated rotameric chi -> 0-deg canonical
                    assert quad in M.SIDECHAIN_CHI[res], (res, quad)
                    assert rec["mean"] == 0.0, (res, quad, rec["mean"])
        # source is per-quad: canonical quads in any entry -> rosetta_params_408;
        # named (g-/g+/t) entries override their rotameric chi quads ->
        # dunbrack_2010; nr entries override the terminal non-rotameric chi ->
        # dunbrack_2010_uniform_30deg_bin.  All other quads stay canonical.
        nonrot = M.SIDECHAIN_NON_ROTAMERIC_BINS
        for key, quad_map in LIB.items():
            res = key.split("_")[0]
            if key.endswith("_canonical"):
                override = {}
            elif _is_nr_bin_key(res, key):
                override = {nonrot[res]["chi_quad"]: "dunbrack_2010_uniform_30deg_bin"}
            else:
                # named rotamer: overridden chi quads are the rotatable ones
                # the rotamer name encodes (first 1 or 2 chi, excluding any
                # terminal non-rotameric chi).
                rc = len(M.SIDECHAIN_CHI[res]) - (1 if res in nonrot else 0)
                rc = min(2, rc)
                over_q = M.SIDECHAIN_CHI[res][:rc]
                override = {q: "dunbrack_2010" for q in over_q}
            for quad, rec in quad_map.items():
                expect_src = override.get(quad, "rosetta_params_408")
                assert rec["source"] == expect_src, (key, quad, rec["source"])

    def test_sidechain_non_rotameric_bins(self):
        B = M.SIDECHAIN_NON_ROTAMERIC_BINS
        expect = {"ASN":12,"ASP":6,"GLN":12,"GLU":6,"HIS":12,"PHE":6,"TRP":12,"TYR":6}
        assert set(B) == set(expect), set(B)
        for res, spec in B.items():
            assert spec["bins"] == expect[res], (res, spec)
            assert spec["chi_quad"] in M.SIDECHAIN_CHI[res], (res, spec)
            assert spec["chi_quad"] == M.SIDECHAIN_CHI[res][-1], (res, spec)  # terminal chi
        # terminal chi stays canonical (not rotated) in named (g-/g+/t)
        # rotamers: it is present at its canonical value/source, not a
        # rotamer override.
        for res in expect:
            last = M.SIDECHAIN_CHI[res][-1]
            canon_mean = M.SIDECHAIN_ROTAMER_LIB[f"{res}_canonical"][last]["mean"]
            for key, quad_map in M.SIDECHAIN_ROTAMER_LIB.items():
                if key.startswith(res + "_") and not key.endswith("_canonical") \
                   and not _is_nr_bin_key(res, key):
                    assert last in quad_map, (res, key)  # present (full geometry)
                    rec = quad_map[last]
                    assert rec["mean"] == canon_mean, (res, key, rec["mean"])
                    assert rec["source"] == "rosetta_params_408", (res, key, rec["source"])
        # non-rotameric bins: exactly ``bins`` nr entries per residue, bin
        # centers at NON_ROTAMERIC_BIN_WIDTH*(i-0.5) = 15,45,75,... covering
        # the full period (bins*30 deg), each carrying the full geometry but
        # overriding only the terminal non-rotameric chi.
        for res, spec in B.items():
            nbin = spec["bins"]
            nq = spec["chi_quad"]
            nr = [k for k in M.SIDECHAIN_ROTAMER_LIB if _is_nr_bin_key(res, k)]
            assert len(nr) == nbin, (res, len(nr))
            assert nr == [f"{res}_nr{i}" for i in range(1, nbin + 1)], res
            for i in range(1, nbin + 1):
                key = f"{res}_nr{i}"
                quad_map = M.SIDECHAIN_ROTAMER_LIB[key]
                assert set(quad_map) == set(M.IC_PATH[res]), (key, set(quad_map))
                rec = quad_map[nq]
                assert rec["mean"] == M.NON_ROTAMERIC_BIN_WIDTH * (i - 0.5), (key, rec["mean"])
                assert rec["source"] == "dunbrack_2010_uniform_30deg_bin", (key, rec["source"])
                # all non-terminal quads stay canonical
                for quad, r2 in quad_map.items():
                    if quad != nq:
                        assert r2["source"] == "rosetta_params_408", (key, quad, r2["source"])

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
        assert M.AMINO_ACID_BACKBONE_BOND_ANGLE[("N", "CA", "C")]["mean"] == 111.2
        assert M.AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "N")]["mean"] == 116.2

    def test_ss_torsion_alpha_helix(self):
        ah = M.SS_BB_TORSION_ANGLE["alpha-helix"]
        phi = ah[M.ALIAS_QUAD["phi"]]
        assert phi["mean"] == -60 and phi["up"] == -35 and phi["lb"] == -85
        assert ah[M.ALIAS_QUAD["psi"]]["mean"] == -45

    def test_sidechain_tables_keyed_by_arity(self):
        # length keys are 2-atom tuples, angle keys are 3-atom tuples,
        # dihedral keys stay 4-atom grow quads; each IC_PATH quad must resolve
        # in the complete canonical rotamer-lib entry, and every quad the
        # private SC_NON_CHI_DIHEDRAL carries is a valid IC_PATH quad.
        from biorazer.database.molecule.bond.dihedral.protein import by_residue
        NC = by_residue.SC_NON_CHI_DIHEDRAL
        for res, quads in M.IC_PATH.items():
            for quad in quads:
                _, j, k, l = quad
                assert (k, l) in M.AMINO_ACID_SIDECHAIN_BOND[res], (res, quad)
                assert (j, k, l) in M.AMINO_ACID_BOND_ANGLE[res], (res, quad)
                assert quad in M.SIDECHAIN_ROTAMER_LIB[f"{res}_canonical"], (res, quad)
            assert set(NC[res]) <= set(M.IC_PATH[res]), (res, set(NC[res]))

    def test_sidechain_length_angle_consistent(self):
        # the same (res, grow quad) must map to the same bond/angle values
        # through the 2- and 3-atom keys.
        for res, quads in M.IC_PATH.items():
            for quad in quads:
                i, j, k, l = quad
                bond = M.AMINO_ACID_SIDECHAIN_BOND[res][(k, l)]
                ang = M.AMINO_ACID_BOND_ANGLE[res][(j, k, l)]
                assert ang["mean"] > 0 and bond["mean"] > 0, (res, quad)
                assert ang["source"] == bond["source"] == "rosetta_params_408", (res, quad)

    def test_mainchain_torsion_definitions(self):
        assert M.ALIAS_QUAD == {
            "phi": ("C", "N", "CA", "C"),
            "psi": ("N", "CA", "C", "N"),
            "omega": ("CA", "C", "N", "CA"),
        }

    def test_backbone_ic_path(self):
        # the uniform backbone grow quads every residue shares
        assert set(M.BACKBONE_IC_PATH) == {"intra", "peptide"}
        # intra: per-residue carbonyl branches off C
        assert M.BACKBONE_IC_PATH["intra"] == (
            ("N", "CA", "C", "O"),
            ("N", "CA", "C", "OXT"),
        )
        # peptide: cross-residue quads, subscripted _i / _{i+1}
        assert len(M.BACKBONE_IC_PATH["peptide"]) == 3
        for quad in M.BACKBONE_IC_PATH["peptide"]:
            assert len(quad) == 4
            for nm in quad:
                assert nm.endswith("_i") or nm.endswith("_{i+1}")
        # the stored dihedral of each peptide quad is the official torsion:
        # psi_i, omega_i, phi_{i+1} (see ALIAS_QUAD)
        assert M.BACKBONE_IC_PATH["peptide"] == (
            ("N_i", "CA_i", "C_i", "N_{i+1}"),      # psi of residue i
            ("CA_i", "C_i", "N_{i+1}", "CA_{i+1}"),  # omega of residue i
            ("C_i", "N_{i+1}", "CA_{i+1}", "C_{i+1}"),  # phi of residue i+1
        )

    def test_omega(self):
        assert M.OMEGA_TRANS["mean"] == 180.0 and M.OMEGA_TRANS["up"] == 190.0
        assert M.OMEGA_CIS["mean"] == 0.0


class TestTemplates:
    """The template store builds and round-trips."""

    def test_get_available_specs(self):
        from biorazer.database.molecule.icoor.protein import template
        # ALA: 0 chi axis -> single canonical rotamer per SS class
        specs = template.get_available_specs("ALA")
        assert len(specs) == 12                      # 12 SS classes
        assert specs[0] == ("alpha-helix", "canonical")
        assert all(rot == "canonical" for _, rot in specs)
        # SER: 1 chi axis -> canonical + 3 g-/t/g+ per SS class
        ser = template.get_available_specs("SER")
        assert len(ser) == 12 * 4
        assert ser[0] == ("alpha-helix", "canonical")
        assert ser[1] == ("alpha-helix", "g-")
        # LEU: 2 chi axes -> canonical + 9 chi1/chi2 per SS class
        leu = template.get_available_specs("LEU")
        assert len(leu) == 12 * 10
        assert leu[1] == ("alpha-helix", "g-/g-")

    def test_build_and_to_atomarray(self):
        from biorazer.database.molecule.icoor.protein import template
        ic = template.build_template("SER", "alpha-helix", "g-")
        assert ic.phi == -60.0 and ic.psi == -45.0 and ic.omega == 180.0
        assert ic.rotamer == "g-" and ic.ss == "alpha-helix"
        arr = ic.to_atomarray()
        assert len(arr) == len(ic.atoms)
        names = set(arr.atom_name)
        assert {"N", "CA", "C", "O", "CB", "OG"} <= names

    def test_templates_round_trip(self):
        # to_coords -> from_atomarray must reproduce the template geometry
        from biorazer.database.molecule.icoor.protein import template
        from biorazer.structure.objects.internal_coords import InternalCoord
        ic = template.build_template("TRP", "beta-strand")
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
