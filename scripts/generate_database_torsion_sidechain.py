# -*- coding: utf-8 -*-
"""Emit the sidechain-torsion dicts (SIDECHAIN_CHI / SIDECHAIN_IC_DIHEDRAL
/ ROTAMER_BIN / DUNBRACK_ROTAMERS) for
biorazer/database/molecule/bond/dihedral/protein.py from Rosetta fa_standard
params.

Run (from the repo root)::

    /opt/envs/BioRazer/bin/python scripts/generate_database_torsion_sidechain.py > /tmp/sc.py

The emitted text contains the four dicts with the uniform
``{mean, std, lb, up, source}`` record (std/lb/up = np.nan); paste them into
the target module.
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
AAS = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
       "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]


def parse(aa):
    icoor, chi = {}, []
    for line in open(PARAMS + aa + ".params", encoding="utf-8"):
        if line.startswith("ICOOR_INTERNAL"):
            p = line.split()
            icoor[p[1]] = (float(p[2]), float(p[3]), float(p[4]), p[5], p[6], p[7])
        elif line.startswith("CHI "):
            p = line.split()
            chi.append(tuple(p[2:6]))
    return icoor, chi


def build(aa):
    ic, _ = parse(aa)
    N = np.array([0.0, 0.0, 0.0])
    CA = N + 1.458 * np.array([1.0, 0.0, 0.0])
    C = CA + 1.523 * np.array([0.361, 0.932, 0.0])
    coord = {"N": N, "CA": CA, "C": C}
    coord["O"] = _place(N, CA, C, 1.231, 120.8, 180.0)
    while True:
        g = False
        for nm, (dih, th, ln, p1, p2, p3) in ic.items():
            if nm in coord or nm[0] not in HEAVY:
                continue
            if p1 in coord and p2 in coord and p3 in coord:
                coord[nm] = _place(coord[p3], coord[p2], coord[p1], ln, 180.0 - th, dih)
                g = True
        if not g:
            break
    return coord


def main():
    out = []
    a = out.append
    a('import numpy as np')
    a('')
    a('')
    a('AAS = %r' % (AAS,))
    a('')
    a('#: 每种残基的 chi 定义 (原子名四元组), chi1..chiN == Dunbrack 库的 chi1..chiN')
    a('SIDECHAIN_CHI = {')
    for aa in AAS:
        _, chi = parse(aa)
        coord = build(aa)  # real heavy atoms actually placed
        hc = [q for q in chi if set(q) <= set(coord) and "NV" not in q]  # drop virtual NV
        if hc:
            a('    %r: %s,' % (aa, repr(hc)))
        else:
            a('    %r: [],' % aa)
    a('}')
    a('')
    a('')
    a('#: 每种残基侧链的规范 IC-frame 理想二面角 (度), keyed by 生长四元组原子名')
    a('#: (与 molecule.bond.length.protein.AMINO_ACID_SIDECHAIN_BOND 同 key)。')
    a('#: GLY 侧链为空。每条为 {mean, std, lb, up, source}; std/lb/up 为 np.nan')
    a('#: (Rosetta ICOOR 只给理想点值)。')
    a('SIDECHAIN_IC_DIHEDRAL = {')
    for aa in AAS:
        coord = build(aa)
        a('    %r: {' % aa)
        for quad in IC_PATH[aa]:
            if all(n in coord for n in quad):
                a('        %r: {"mean": %.2f, "std": np.nan, "lb": np.nan, '
                  '"up": np.nan, "source": "rosetta_params_408"},'
                  % (quad, _dihedral(*(coord[n] for n in quad))))
        a('    },' )
    a('}')
    a('')
    a('#: 标准 rotamer bin 中心 (理想化 chi 定义), 度。g-/t/g+ 命名与这些中心是')
    a('#: Dunbrack rotamer 库的通用归类; 完整骨架依赖数值表见模块 docstring 说明。')
    a('ROTAMER_BIN = {')
    for name, val in (("g-", -60.0), ("t", 180.0), ("g+", 60.0)):
        a('    %r: {"mean": %.1f, "std": np.nan, "lb": np.nan, "up": np.nan, '
          '"source": "dunbrack_2010"},' % (name, val))
    a('}')
    a('')
    a('#: per 残基的 rotamer 分类框架: 可旋转 chi 轴数 (见 SIDECHAIN_CHI) + 标准 bin 中心。')
    a('#: 完整数值表 (逐 phi/psi) 未内嵌。')
    a('DUNBRACK_ROTAMERS = {')
    for aa in AAS:
        _, chi = parse(aa)
        coord = build(aa)
        n = len([q for q in chi if set(q) <= set(coord) and "NV" not in q])
        a('    %r: {"chi": %d, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)),'
          ' "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},' % (aa, n))
    a('}')
    a('')
    a('SIDECHAIN_DIHE_REFS = {')
    a('    "rosetta_params_408": "Rosetta 408 ... l-caa/*.params (ICOOR_INTERNAL 规范残基理想几何)",')
    a('    "dunbrack_2010": "Shapovalov MV, Dunbrack RL Jr. A smoothed backbone-dependent '
      'rotamer library. Structure 19:844-858, 2011."')
    a('}')
    return "\n".join(out)


if __name__ == "__main__":
    print(main())