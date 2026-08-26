# -*- coding: utf-8 -*-
"""Emit biorazer/database/bond/sidechain.py from Rosetta fa_standard params."""
import numpy as np
from biorazer.structure.objects.internal_coords import _place, _dihedral
from biorazer.database.internal_coord_template._topology import SIDE_CHAIN_IC_PATH

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
    w('# -*- coding: utf-8 -*-')
    w('"""Amino-acid side-chain covalent bond lengths (Å) and bond angles (°).')
    w("")
    w('数据集来源与诚实性说明')
    w('───────────────────────')
    w('本模块记录 20 种标准氨基酸**侧链**(CA 的 CB 及以外重原子, 见')
    w('`internal_coord_template._topology.SIDE_CHAIN_IC_PATH` 的生长路径)的理想键长与键角。')
    w('数值为**理想点值 (ideal point values)**, 翻译自 Rosetta 408 的')
    w('`fa_standard` 残基 .params (`ICOOR_INTERNAL`) 构建的规范残基几何 ——')
    w('Rosetta 的主链/侧链理想几何源自 Engh & Huber (1991) (见')
    w('`bond.backbone` 的交叉核对)。键长即 ICOOR 的 `d`, 键角即 `180° − theta`')
    w('(解码后又在构建出的残基坐标上实测复核)。')
    w("")
    w('注意: 这里**没有**给出 Engh & Huber 的样本标准差 (std) —— Rosetta ICOOR')
    w('只给理想点值, 不给 CSD sigma; 需要 std 时请以 E&H 1991 原始文献为准。')
    w('键长单位 Å, 键角单位 **度 (degree)**。')
    w("")
    w('键/角的 key 采用 `SIDE_CHAIN_IC_PATH` 的生长四元组 (i,j,k,l): 从父原子')
    w('(i,j,k) 生长出 l, 键 `(k,l)`, 键角 `(j,k,l)` (顶点在 k)。')
    w('"""')
    w("")
    w("")
    w('#: 每种氨基酸侧链生长的规范键长 (k,l) 与键角 (j,k,l)')
    w('#: {(res_name): {(i,j,k,l)取原子名的三元/四元: {...}}}')
    w('#: 含第一个 CB 生长四元组 (N,C,CA,CB) (其 CA-CB 键/角度亦见 bond.backbone)')
    w('AMINO_ACID_SIDECHAIN_BOND = {')
    for aa in AAS:
        coord, _ = build_canonical(aa)
        w(f'    "{aa}": {{')
        for quad in SIDE_CHAIN_IC_PATH[aa]:
            if not all(n in coord for n in quad):
                continue
            blen, bang, _ = quad_geom(aa, coord, quad)
            i, j, k, l = quad
            quadkey = "(" + ", ".join(repr(n) for n in quad) + ")"
            w('        %s: {"bond": (%r, %r), "value": %.4f, "mean": %.4f, '
              '"angle_pair": (%r, %r, %r), "angle": %.2f, '
              '"note": "Rosetta fa_standard %s.params ICOOR ideal", '
              '"source": "rosetta_params_408"},'
              % (quadkey, k, l, blen, blen, j, k, l, bang, aa))
        w('    },')
    w('}')
    w("")
    w('BOND_SIDECHAIN_REFS = {')
    w('    "rosetta_params_408": "Rosetta 408 main/database/chemical/residue_type_sets/'
      'fa_standard/residue_types/l-caa/*.params (ICOOR_INTERNAL 理想键长/键角)。",')
    w('}')
    return "\n".join(lines)


if __name__ == "__main__":
    print(main())