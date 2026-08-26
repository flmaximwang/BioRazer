# -*- coding: utf-8 -*-
"""Amino-acid side-chain torsion data: chi attenuation-free ideal dihedrals.

本模块放侧链**扭转角**数据 (与 `bond.sidechain` 的键长/键角对应):

* `SIDECHAIN_CHI`            -- 每种残基的官方 chi 扭转角定义 (chi1..chiN 的
   原子名四元组)，取自 Rosetta 408 fa_standard .params 的 CHI 行，只保留
   重原子四元组。这些就是 Dunbrack 库所用的 chi 定义。
* `SIDECHAIN_IC_DIHEDRAL`    -- 每种残基侧链在 `SIDE_CHAIN_IC_PATH` 生长四元组
   参考系下的**规范理想二面角** (度)。由 Rosetta 408 .params
   (ICOOR_INTERNAL) 构建的规范残基几何实测 —— 与 `bond.sidechain` 的键长/键角
   取自同一次规范构建，因此配套使用即可把侧链拼出一个自洽 (可 to_coords 重建)
   的规范侧面链。注意 IC-frame 二面角随参考骨架/chi 构象轻微耦合 (见
   `internal_coord_template` 的 docstring)，故它是**一个规范构象**的快照，
   不是对任何骨架都成立的万能不变常量。
* `ROTAMER_BIN` / `DUNBRACK_ROTAMERS` -- 侧链 rotamer 的**分类框架**: 标准
   rotamer bin 中心 (g-/t/g+ = -60/180/+60) 和每种残基的可旋转 chi 轴数。
   这是 Dunbrack 库所用的 rotamer 命名/归类定义; 完整的**骨架依赖数值表**
   (逐 phi/psi bin 的均值/方差) 属于外部 Dunbrack 2010 数据集
   (Shapovalov & Dunbrack 2011, CC BY 4.0), 未内嵌, 需要时可另行 vendor。

所有角度单位 **度 (degree)**。
"""


AAS = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

#: 每种残基的 chi 定义 (原子名四元组), chi1..chiN == Dunbrack 库的 chi1..chiN
SIDECHAIN_CHI = {
    'ALA': [],
    'ARG': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD'), ('CB', 'CG', 'CD', 'NE'), ('CG', 'CD', 'NE', 'CZ')],
    'ASN': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'OD1')],
    'ASP': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'OD1')],
    'CYS': [('N', 'CA', 'CB', 'SG')],
    'GLN': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD'), ('CB', 'CG', 'CD', 'OE1')],
    'GLU': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD'), ('CB', 'CG', 'CD', 'OE1')],
    'GLY': [],
    'HIS': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'ND1')],
    'ILE': [('N', 'CA', 'CB', 'CG1'), ('CA', 'CB', 'CG1', 'CD1')],
    'LEU': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD1')],
    'LYS': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD'), ('CB', 'CG', 'CD', 'CE'), ('CG', 'CD', 'CE', 'NZ')],
    'MET': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'SD'), ('CB', 'CG', 'SD', 'CE')],
    'PHE': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD1')],
    'PRO': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD')],
    'SER': [('N', 'CA', 'CB', 'OG')],
    'THR': [('N', 'CA', 'CB', 'OG1')],
    'TRP': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD1')],
    'TYR': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD1')],
    'VAL': [('N', 'CA', 'CB', 'CG1')],
}


#: 每种残基侧链的规范 IC-frame 理想二面角 (度), keyed by 生长四元组原子名
#: (与 bond.sidechain.AMINO_ACID_SIDECHAIN_BOND 同 key)。GLY 侧链为空。
SIDECHAIN_IC_DIHEDRAL = {
    'ALA': {
        ('C', 'N', 'CA', 'CB'): -122.80,
    },
    'ARG': {
        ('C', 'N', 'CA', 'CB'): -121.70,
        ('N', 'CA', 'CB', 'CG'): 0.00,
        ('CA', 'CB', 'CG', 'CD'): -0.00,
        ('CB', 'CG', 'CD', 'NE'): -0.00,
        ('CG', 'CD', 'NE', 'CZ'): 0.00,
        ('CD', 'NE', 'CZ', 'NH1'): 0.00,
        ('NH1', 'NE', 'CZ', 'NH2'): -180.00,
    },
    'ASN': {
        ('C', 'N', 'CA', 'CB'): -123.30,
        ('N', 'CA', 'CB', 'CG'): -0.00,
        ('CA', 'CB', 'CG', 'OD1'): 0.00,
        ('OD1', 'CB', 'CG', 'ND2'): -180.00,
    },
    'ASP': {
        ('C', 'N', 'CA', 'CB'): -122.10,
        ('N', 'CA', 'CB', 'CG'): 0.00,
        ('CA', 'CB', 'CG', 'OD1'): -0.00,
        ('OD1', 'CB', 'CG', 'OD2'): 180.00,
    },
    'CYS': {
        ('C', 'N', 'CA', 'CB'): -121.60,
        ('N', 'CA', 'CB', 'SG'): -0.00,
    },
    'GLN': {
        ('C', 'N', 'CA', 'CB'): -122.40,
        ('N', 'CA', 'CB', 'CG'): -0.00,
        ('CA', 'CB', 'CG', 'CD'): -0.00,
        ('CB', 'CG', 'CD', 'OE1'): -0.00,
        ('OE1', 'CG', 'CD', 'NE2'): -180.00,
    },
    'GLU': {
        ('C', 'N', 'CA', 'CB'): -122.20,
        ('N', 'CA', 'CB', 'CG'): -0.00,
        ('CA', 'CB', 'CG', 'CD'): -0.00,
        ('CB', 'CG', 'CD', 'OE1'): -0.00,
        ('OE1', 'CG', 'CD', 'OE2'): -180.00,
    },
    'GLY': {
    },
    'HIS': {
        ('C', 'N', 'CA', 'CB'): -122.20,
        ('N', 'CA', 'CB', 'CG'): -0.00,
        ('CA', 'CB', 'CG', 'ND1'): 0.00,
        ('CB', 'CG', 'ND1', 'CE1'): 180.00,
        ('CG', 'ND1', 'CE1', 'NE2'): 0.04,
        ('ND1', 'CE1', 'NE2', 'CD2'): -0.03,
    },
    'ILE': {
        ('C', 'N', 'CA', 'CB'): -122.10,
        ('N', 'CA', 'CB', 'CG1'): 0.00,
        ('CA', 'CB', 'CG1', 'CD1'): 0.00,
        ('CG1', 'CA', 'CB', 'CG2'): -122.68,
    },
    'LEU': {
        ('C', 'N', 'CA', 'CB'): -121.90,
        ('N', 'CA', 'CB', 'CG'): 0.00,
        ('CA', 'CB', 'CG', 'CD1'): 0.00,
        ('CD1', 'CB', 'CG', 'CD2'): 122.00,
    },
    'LYS': {
        ('C', 'N', 'CA', 'CB'): -122.60,
        ('N', 'CA', 'CB', 'CG'): -0.00,
        ('CA', 'CB', 'CG', 'CD'): 0.00,
        ('CB', 'CG', 'CD', 'CE'): 0.00,
        ('CG', 'CD', 'CE', 'NZ'): -0.00,
    },
    'MET': {
        ('C', 'N', 'CA', 'CB'): -122.80,
        ('N', 'CA', 'CB', 'CG'): 0.00,
        ('CA', 'CB', 'CG', 'SD'): 0.00,
        ('CB', 'CG', 'SD', 'CE'): -0.00,
    },
    'PHE': {
        ('C', 'N', 'CA', 'CB'): -122.00,
        ('N', 'CA', 'CB', 'CG'): -0.00,
        ('CA', 'CB', 'CG', 'CD1'): 0.00,
        ('CB', 'CG', 'CD1', 'CE1'): -179.99,
        ('CG', 'CD1', 'CE1', 'CZ'): -0.03,
        ('CD1', 'CE1', 'CZ', 'CE2'): 0.03,
        ('CE1', 'CZ', 'CE2', 'CD2'): 0.03,
    },
    'PRO': {
        ('C', 'N', 'CA', 'CB'): -119.70,
        ('N', 'CA', 'CB', 'CG'): 30.00,
        ('CA', 'CB', 'CG', 'CD'): -33.90,
    },
    'SER': {
        ('C', 'N', 'CA', 'CB'): -122.00,
        ('N', 'CA', 'CB', 'OG'): -0.00,
    },
    'THR': {
        ('C', 'N', 'CA', 'CB'): -122.40,
        ('N', 'CA', 'CB', 'OG1'): -0.00,
        ('OG1', 'CA', 'CB', 'CG2'): -120.54,
    },
    'TRP': {
        ('C', 'N', 'CA', 'CB'): -121.90,
        ('N', 'CA', 'CB', 'CG'): 0.00,
        ('CA', 'CB', 'CG', 'CD1'): 0.00,
        ('CB', 'CG', 'CD1', 'NE1'): -179.97,
        ('CG', 'CD1', 'NE1', 'CE2'): -0.11,
        ('CD1', 'NE1', 'CE2', 'CZ2'): 179.96,
        ('NE1', 'CE2', 'CZ2', 'CH2'): -179.97,
        ('CE2', 'CZ2', 'CH2', 'CZ3'): 0.13,
        ('CZ2', 'CH2', 'CZ3', 'CE3'): -0.07,
        ('CH2', 'CZ3', 'CE3', 'CD2'): 0.02,
    },
    'TYR': {
        ('C', 'N', 'CA', 'CB'): -122.80,
        ('N', 'CA', 'CB', 'CG'): -0.00,
        ('CA', 'CB', 'CG', 'CD1'): 0.00,
        ('CB', 'CG', 'CD1', 'CE1'): -180.00,
        ('CD1', 'CB', 'CG', 'CD2'): 180.00,
        ('CB', 'CG', 'CD2', 'CE2'): -180.00,
        ('CG', 'CD2', 'CE2', 'CZ'): -0.00,
        ('CD2', 'CE2', 'CZ', 'OH'): 180.00,
    },
    'VAL': {
        ('C', 'N', 'CA', 'CB'): -121.50,
        ('N', 'CA', 'CB', 'CG1'): -0.00,
        ('CG1', 'CA', 'CB', 'CG2'): 122.54,
    },
}

#: 标准 rotamer bin 中心 (理想化 chi 定义), 度。g-/t/g+ 命名与这些中心是
#: Dunbrack rotamer 库的通用归类; 完整骨架依赖数值表见模块 docstring 说明。
ROTAMER_BIN = {"g-": -60.0, "t": 180.0, "g+": 60.0}

#: per 残基的 rotamer 分类框架: 可旋转 chi 轴数 (见 SIDECHAIN_CHI) + 标准 bin 中心。
#: 完整数值表 (逐 phi/psi) 未内嵌。
DUNBRACK_ROTAMERS = {
    'ALA': {"chi": 0, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'ARG': {"chi": 4, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'ASN': {"chi": 2, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'ASP': {"chi": 2, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'CYS': {"chi": 1, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'GLN': {"chi": 3, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'GLU': {"chi": 3, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'GLY': {"chi": 0, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'HIS': {"chi": 2, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'ILE': {"chi": 2, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'LEU': {"chi": 2, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'LYS': {"chi": 4, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'MET': {"chi": 3, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'PHE': {"chi": 2, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'PRO': {"chi": 2, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'SER': {"chi": 1, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'THR': {"chi": 1, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'TRP': {"chi": 2, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'TYR': {"chi": 2, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'VAL': {"chi": 1, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
}

SIDECHAIN_DIHE_REFS = {
    "rosetta_params_408": "Rosetta 408 ... l-caa/*.params (ICOOR_INTERNAL 规范残基理想几何)",
    "dunbrack_2010": "Shapovalov MV, Dunbrack RL Jr. A smoothed backbone-dependent rotamer library. Structure 19:844-858, 2011."
}