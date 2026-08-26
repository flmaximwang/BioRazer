# -*- coding: utf-8 -*-
"""Emit biorazer/database/torsion_angle/sidechain.py: chi definitions +
canonical IC-path sidechain dihedrals (measured from Rosetta ideal build)."""
import numpy as np
from biorazer.structure.objects.internal_coords import _place, _dihedral
from biorazer.database.internal_coord_template._topology import SIDE_CHAIN_IC_PATH

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
                coord[nm] = _place(coord[p3], coord[p2], coord[p1], ln, 180.0 - th, -dih)
                g = True
        if not g:
            break
    return coord


def main():
    out = []
    a = out.append
    a('# -*- coding: utf-8 -*-')
    a('"""Amino-acid side-chain torsion data: chi attenuation-free ideal dihedrals.')
    a("")
    a('本模块放侧链**扭转角**数据 (与 `bond.sidechain` 的键长/键角对应):')
    a("")
    a('* `SIDECHAIN_CHI`            -- 每种残基的官方 chi 扭转角定义 (chi1..chiN 的')
    a('   原子名四元组)，取自 Rosetta 408 fa_standard .params 的 CHI 行，只保留')
    a('   重原子四元组。这些就是 Dunbrack 库所用的 chi 定义。')
    a('* `SIDECHAIN_IC_DIHEDRAL`    -- 每种残基侧链在 `SIDE_CHAIN_IC_PATH` 生长四元组')
    a('   参考系下的**规范理想二面角** (度)。由 Rosetta 408 .params')
    a('   (ICOOR_INTERNAL) 构建的规范残基几何实测 —— 与 `bond.sidechain` 的键长/键角')
    a('   取自同一次规范构建，因此配套使用即可把侧链拼出一个自洽 (可 to_coords 重建)')
    a('   的规范侧面链。注意 IC-frame 二面角随参考骨架/chi 构象轻微耦合 (见')
    a('   `internal_coord_template` 的 docstring)，故它是**一个规范构象**的快照，')
    a('   不是对任何骨架都成立的万能不变常量。')
    a('* `ROTAMER_BIN` / `DUNBRACK_ROTAMERS` -- 侧链 rotamer 的**分类框架**: 标准')
    a('   rotamer bin 中心 (g-/t/g+ = -60/180/+60) 和每种残基的可旋转 chi 轴数。')
    a('   这是 Dunbrack 库所用的 rotamer 命名/归类定义; 完整的**骨架依赖数值表**')
    a('   (逐 phi/psi bin 的均值/方差) 属于外部 Dunbrack 2010 数据集')
    a('   (Shapovalov & Dunbrack 2011, CC BY 4.0), 未内嵌, 需要时可另行 vendor。')
    a("")
    a('所有角度单位 **度 (degree)**。')
    a('"""')
    a("")
    a("")
    a('AAS = %r' % (AAS,))
    a("")
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
    a("")
    a("")
    a('#: 每种残基侧链的规范 IC-frame 理想二面角 (度), keyed by 生长四元组原子名')
    a('#: (与 bond.sidechain.AMINO_ACID_SIDECHAIN_BOND 同 key)。GLY 侧链为空。')
    a('SIDECHAIN_IC_DIHEDRAL = {')
    for aa in AAS:
        coord = build(aa)
        a('    %r: {' % aa)
        for quad in SIDE_CHAIN_IC_PATH[aa]:
            if all(n in coord for n in quad):
                a('        %r: %.2f,' % (quad, _dihedral(*(coord[n] for n in quad))))
        a('    },')
    a('}')
    a("")
    a('#: 标准 rotamer bin 中心 (理想化 chi 定义), 度。g-/t/g+ 命名与这些中心是')
    a('#: Dunbrack rotamer 库的通用归类; 完整骨架依赖数值表见模块 docstring 说明。')
    a('ROTAMER_BIN = {"g-": -60.0, "t": 180.0, "g+": 60.0}')
    a("")
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
    a("")
    a('SIDECHAIN_DIHE_REFS = {')
    a('    "rosetta_params_408": "Rosetta 408 ... l-caa/*.params (ICOOR_INTERNAL 规范残基理想几何)",')
    a('    "dunbrack_2010": "Shapovalov MV, Dunbrack RL Jr. A smoothed backbone-dependent '
      'rotamer library. Structure 19:844-858, 2011."')
    a('}')
    return "\n".join(out)


if __name__ == "__main__":
    print(main())