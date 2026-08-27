# -*- coding: utf-8 -*-
"""Emit the protein side-chain bond-length / bond-angle dicts from Rosetta
fa_standard params.

Targets:
* ``AMINO_ACID_SIDECHAIN_BOND`` -> biorazer/database/molecule/bond/length/
  protein/by_residue.py  (the merged length module; paste after the
  ``AMINO_ACID_BOND_LENGTH_BY_RESIDUE`` dict)
* ``AMINO_ACID_SIDECHAIN_BOND_ANGLE`` -> the **side-chain section** of
  ``AMINO_ACID_BOND_ANGLE`` in biorazer/database/molecule/bond/angle/
  protein.py (the merged residue-keyed angle module)

Run (from the repo root)::

    /opt/envs/BioRazer/bin/python scripts/generate_database_bond_sidechain.py > /tmp/sc.py

The emitted text contains only the ``AMINO_ACID_SIDECHAIN_BOND`` (keyed
by 2-atom ``(k, l)`` bond) and ``AMINO_ACID_SIDECHAIN_BOND_ANGLE`` (keyed
by 3-atom ``(j, k, l)`` angle) dicts (the ``{mean, std, lb, up, source,
...}`` uniform record; std/lb/up = np.nan); paste the side-chain angle dict
as the side-chain section of the merged ``AMINO_ACID_BOND_ANGLE`` (the
backbone per-residue refinements Gly/Pro/Ala/VIT sit above it).
"""
import numpy as np
import sys
import pathlib

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from biorazer.structure.objects.internal_coords import _place, _dihedral
from biorazer.database.molecule.icoor.protein.topology import IC_PATH

PARAMS = "/opt/Rosetta/rosetta.source.release-408/main/database/chemical/residue_type_sets/fa_standard/residue_types/l-caa/"
HEAVY = set("CNOS")
AAS = ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
       "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]


def parse_icoor(aa):
    icoor = {}
    for line in open(PARAMS + aa + ".params", encoding="utf-8"):
        if line.startswith("ICOOR_INTERNAL"):
            p = line.split()
            icoor[p[1]] = (float(p[2]), float(p[3]), float(p[4]), p[5], p[6], p[7])
    return icoor


def build_canonical(aa):
    icoor = parse_icoor(aa)
    N = np.array([0.0, 0.0, 0.0])
    CA = N + 1.458 * np.array([1.0, 0.0, 0.0])
    C = CA + 1.523 * np.array([0.361, 0.932, 0.0])
    coord = {"N": N, "CA": CA, "C": C}
    coord["O"] = _place(N, CA, C, 1.231, 120.8, 180.0)
    while True:
        grew = False
        for nm, (dih, theta, ln, p1, p2, p3) in icoor.items():
            if nm in coord or nm[0] not in HEAVY:
                continue
            if p1 in coord and p2 in coord and p3 in coord:
                coord[nm] = _place(coord[p3], coord[p2], coord[p1], ln, 180.0 - theta, dih)
                grew = True
        if not grew:
            break
    return coord, icoor


def quad_geom(aa, coord, quad):
    i, j, k, l = quad
    ci, cj, ck, cl = (coord[n] for n in quad)
    blen = float(np.linalg.norm(cl - ck))
    cos = np.dot(cj - ck, cl - ck) / (np.linalg.norm(cj - ck) * np.linalg.norm(cl - ck))
    bang = float(np.degrees(np.arccos(np.clip(cos, -1, 1))))
    dih = _dihedral(ci, cj, ck, cl)
    return blen, bang, dih


def main():
    lines = []
    w = lines.append
    w('import numpy as np')
    w('')
    w('')
    w('#: 每种氨基酸侧链的规范键长 (k,l), key 为 2 原子元组')
    w('#: {(res_name): {(k, l): {...}}}')
    w('#: 含第一条 CB 键 (CA,CB) (其键长亦见 generic 主表)。')
    w('#: std/lb/up 为 np.nan (Rosetta ICOOR 只给理想点值, 无 CSD sigma)。')
    w('AMINO_ACID_SIDECHAIN_BOND = {')
    for aa in AAS:
        coord, _ = build_canonical(aa)
        w(f'    "{aa}": {{')
        for quad in IC_PATH[aa]:
            if not all(n in coord for n in quad):
                continue
            blen, _, _ = quad_geom(aa, coord, quad)
            i, j, k, l = quad
            w('        (%r, %r): {"mean": %.4f, "std": np.nan, '
              '"lb": np.nan, "up": np.nan, '
              '"note": "Rosetta fa_standard %s.params ICOOR ideal", '
              '"source": "rosetta_params_408"},'
              % (k, l, blen, aa))
        w('    },')
    w('}')
    w('')
    w('')
    w('#: 每种氨基酸侧链的规范键角 (j,k,l), key 为 3 原子元组, 顶点在 k')
    w('#: {(res_name): {(j, k, l): {mean, std, lb, up, source}}}')
    w('#: std/lb/up 为 np.nan (Rosetta ICOOR 只给理想点值, 无 CSD sigma)。')
    w('#: 这是合并后 AMINO_ACID_BOND_ANGLE 的「侧链」段 (粘贴到骨架细分段之后)。')
    w('AMINO_ACID_SIDECHAIN_BOND_ANGLE = {')
    for aa in AAS:
        coord, _ = build_canonical(aa)
        w(f'    "{aa}": {{')
        for quad in IC_PATH[aa]:
            if not all(n in coord for n in quad):
                continue
            _, bang, _ = quad_geom(aa, coord, quad)
            i, j, k, l = quad
            w('        (%r, %r, %r): {"mean": %.2f, '
              '"std": np.nan, "lb": np.nan, "up": np.nan, '
              '"note": "Rosetta fa_standard %s.params ICOOR ideal", '
              '"source": "rosetta_params_408"},'
              % (j, k, l, bang, aa))
        w('    },')
    w('}')
    w('')
    w('BOND_SIDECHAIN_REFS = {')
    w('    "rosetta_params_408": "Rosetta 408 main/database/chemical/'
      'residue_type_sets/fa_standard/residue_types/l-caa/*.params '
      '(ICOOR_INTERNAL 理想键长/键角)。",')
    w('}')
    return "\n".join(lines)


if __name__ == "__main__":
    print(main())