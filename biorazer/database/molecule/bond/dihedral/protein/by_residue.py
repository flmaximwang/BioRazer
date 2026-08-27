# -*- coding: utf-8 -*-
"""Per-residue (side-chain) torsion data: chi definitions and rotamer framework.

数据集来源与诚实性说明
────────────────────────
本模块记录蛋白质**侧链** (side-chain) 扭转角数据: 每种残基的 chi 扭转角
定义、规范 IC-frame 理想二面角与 rotamer 分类框架。

* ``SIDECHAIN_CHI`` -- 每种残基的官方 chi 扭转角定义 (chi1..chiN 的
  原子名四元组), 取自 Rosetta 408 fa_standard .params 的 CHI 行, 只保留
  重原子四元组。这些就是 Dunbrack 库所用的 chi 定义。
* ``SIDECHAIN_IC_DIHEDRAL`` -- 每种残基侧链在 :data:`IC_PATH` 生长四元组
  参考系下的**规范理想二面角** (度)。由 Rosetta 408 .params
  (ICOOR_INTERNAL) 构建的规范残基几何实测 —— 与
  :mod:`biorazer.database.molecule.bond.length.protein` 的键长/键角取自
  同一次规范构建, 因此配套使用即可把侧链拼出一个自洽 (可 to_coords
  重建) 的规范侧链。注意 IC-frame 二面角随参考骨架/chi 构象轻微耦合,
  故它是**一个规范构象**的快照, 不是对任何骨架都成立的万能不变常量。
* ``ROTAMER_BIN`` / ``DUNBRACK_ROTAMERS`` -- 侧链 rotamer 的**分类框架**:
  标准 rotamer bin 中心 (g-/t/g+ = -60/180/+60) 和每种残基的可旋转 chi
  轴数。完整的**骨架依赖数值表** (逐 phi/psi bin 的均值/方差) 属于外部
  Dunbrack 2010 数据集 (Shapovalov & Dunbrack 2011, CC BY 4.0), 未内嵌,
  需要时可另行 vendor。
* ``NON_ROTAMERIC_BIN_WIDTH`` -- 非 rotameric 末端 chi (ASN/ASP/GLN/GLU/
  PHE/TRP/HIS/TYR 的 sp2 末端) 的 30 deg 细 bin 宽度 (度)。这类 chi 不是
  g-/t/g+ 的离散 rotamer, 而是宽而对称性差的连续分布; 论文用 kernel 密度
  估计建模, 并额外提供 30 deg 离散 bin (每 bin 的均值/方差/占比) 以兼容
  SCWRL 等旧应用。bin 数见 ``DUNBRACK_ROTAMERS`` 的 ``non_rotameric_bins``。

所有角度单位 **度 (degree)**。数值记录统一为 ``{mean, std, lb, up,
source}``; 查不到 spread 的字段为 ``np.nan`` (Rosetta ICOOR 只给理想点值)。

主链 (按二级结构分类的 phi/psi/omega) 数据见 :mod:`.by_ss`。
"""

import numpy as np

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
#: (i,j,k,l), 即官方二面角定义 (如 chi1=(N,CA,CB,CG))。GLY 侧链为空。
#: 每条为 {mean, std, lb, up, source}; Rosetta ICOOR 只给理想点值,
#: 故 std/lb/up = np.nan。键长/键角见 bond.length/angle.protein 的侧链表
#: (key 分别为 2/3 原子元组)。
SIDECHAIN_IC_DIHEDRAL = {
    'ALA': {
        ('C', 'N', 'CA', 'CB'): {"mean": -122.80, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
    },
    'ARG': {
        ('C', 'N', 'CA', 'CB'): {"mean": -121.70, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"mean": 0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'CD'): {"mean": -0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD', 'NE'): {"mean": -0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CG', 'CD', 'NE', 'CZ'): {"mean": 0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CD', 'NE', 'CZ', 'NH1'): {"mean": 0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('NH1', 'NE', 'CZ', 'NH2'): {"mean": -180.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
    },
    'ASN': {
        ('C', 'N', 'CA', 'CB'): {"mean": -123.30, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"mean": -0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'OD1'): {"mean": 0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('OD1', 'CB', 'CG', 'ND2'): {"mean": -180.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
    },
    'ASP': {
        ('C', 'N', 'CA', 'CB'): {"mean": -122.10, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"mean": 0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'OD1'): {"mean": -0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('OD1', 'CB', 'CG', 'OD2'): {"mean": 180.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
    },
    'CYS': {
        ('C', 'N', 'CA', 'CB'): {"mean": -121.60, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'SG'): {"mean": -0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
    },
    'GLN': {
        ('C', 'N', 'CA', 'CB'): {"mean": -122.40, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"mean": -0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'CD'): {"mean": -0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD', 'OE1'): {"mean": -0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('OE1', 'CG', 'CD', 'NE2'): {"mean": -180.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
    },
    'GLU': {
        ('C', 'N', 'CA', 'CB'): {"mean": -122.20, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"mean": -0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'CD'): {"mean": -0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD', 'OE1'): {"mean": -0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('OE1', 'CG', 'CD', 'OE2'): {"mean": -180.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
    },
    'GLY': {
    },
    'HIS': {
        ('C', 'N', 'CA', 'CB'): {"mean": -122.20, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"mean": -0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'ND1'): {"mean": 0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CB', 'CG', 'ND1', 'CE1'): {"mean": 180.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CG', 'ND1', 'CE1', 'NE2'): {"mean": 0.04, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('ND1', 'CE1', 'NE2', 'CD2'): {"mean": -0.03, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
    },
    'ILE': {
        ('C', 'N', 'CA', 'CB'): {"mean": -122.10, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG1'): {"mean": 0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG1', 'CD1'): {"mean": 0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CG1', 'CA', 'CB', 'CG2'): {"mean": -122.68, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
    },
    'LEU': {
        ('C', 'N', 'CA', 'CB'): {"mean": -121.90, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"mean": 0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'CD1'): {"mean": 0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CD1', 'CB', 'CG', 'CD2'): {"mean": 122.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
    },
    'LYS': {
        ('C', 'N', 'CA', 'CB'): {"mean": -122.60, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"mean": -0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'CD'): {"mean": 0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD', 'CE'): {"mean": 0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CG', 'CD', 'CE', 'NZ'): {"mean": -0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
    },
    'MET': {
        ('C', 'N', 'CA', 'CB'): {"mean": -122.80, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"mean": 0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'SD'): {"mean": 0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CB', 'CG', 'SD', 'CE'): {"mean": -0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
    },
    'PHE': {
        ('C', 'N', 'CA', 'CB'): {"mean": -122.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"mean": -0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'CD1'): {"mean": 0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD1', 'CE1'): {"mean": -179.99, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CG', 'CD1', 'CE1', 'CZ'): {"mean": -0.03, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CD1', 'CE1', 'CZ', 'CE2'): {"mean": 0.03, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CE1', 'CZ', 'CE2', 'CD2'): {"mean": 0.03, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
    },
    'PRO': {
        ('C', 'N', 'CA', 'CB'): {"mean": -119.70, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"mean": 30.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'CD'): {"mean": -33.90, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
    },
    'SER': {
        ('C', 'N', 'CA', 'CB'): {"mean": -122.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'OG'): {"mean": -0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
    },
    'THR': {
        ('C', 'N', 'CA', 'CB'): {"mean": -122.40, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'OG1'): {"mean": -0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('OG1', 'CA', 'CB', 'CG2'): {"mean": -120.54, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
    },
    'TRP': {
        ('C', 'N', 'CA', 'CB'): {"mean": -121.90, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"mean": 0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'CD1'): {"mean": 0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD1', 'NE1'): {"mean": -179.97, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CG', 'CD1', 'NE1', 'CE2'): {"mean": -0.11, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CD1', 'NE1', 'CE2', 'CZ2'): {"mean": 179.96, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('NE1', 'CE2', 'CZ2', 'CH2'): {"mean": -179.97, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CE2', 'CZ2', 'CH2', 'CZ3'): {"mean": 0.13, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CZ2', 'CH2', 'CZ3', 'CE3'): {"mean": -0.07, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CH2', 'CZ3', 'CE3', 'CD2'): {"mean": 0.02, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
    },
    'TYR': {
        ('C', 'N', 'CA', 'CB'): {"mean": -122.80, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"mean": -0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'CD1'): {"mean": 0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD1', 'CE1'): {"mean": -180.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CD1', 'CB', 'CG', 'CD2'): {"mean": 180.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD2', 'CE2'): {"mean": -180.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CG', 'CD2', 'CE2', 'CZ'): {"mean": -0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CD2', 'CE2', 'CZ', 'OH'): {"mean": 180.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
    },
    'VAL': {
        ('C', 'N', 'CA', 'CB'): {"mean": -121.50, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG1'): {"mean": -0.00, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
        ('CG1', 'CA', 'CB', 'CG2'): {"mean": 122.54, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "rosetta_params_408"},
    },
}

#: 标准 rotamer bin 中心 (理想化 chi 定义), 度。g-/t/g+ 命名与这些中心是
#: Dunbrack rotamer 库的通用归类; 完整骨架依赖数值表见模块 docstring 说明。
ROTAMER_BIN = {
    "g-": {"mean": -60.0, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "dunbrack_2010"},
    "t":  {"mean": 180.0, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "dunbrack_2010"},
    "g+": {"mean": 60.0, "std": np.nan, "lb": np.nan, "up": np.nan, "source": "dunbrack_2010"},
}

#: 非 rotameric 末端 chi (sp2 末端) 的 30 deg 细 bin 宽度 (度)。
NON_ROTAMERIC_BIN_WIDTH = 30.0

#: per 残基的 rotamer 分类框架。rotameric_chi = 前导可 rotameric chi 轴数
#: (见 SIDECHAIN_CHI); terminal_non_rotameric = 末端 chi 是否为非 rotameric
#: (sp2 杂化, 用 30 deg 细 bin 连续分布建模); non_rotameric_bins = 末端非
#: rotameric chi 的 30 deg bin 数 (None 表示 rotameric)。完整数值表
#: (逐 phi/psi) 未内嵌。
DUNBRACK_ROTAMERS = {
    'ALA': {"chi": 0, "rotameric_chi": 0, "terminal_non_rotameric": False, "non_rotameric_bins": None, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'ARG': {"chi": 4, "rotameric_chi": 4, "terminal_non_rotameric": False, "non_rotameric_bins": None, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'ASN': {"chi": 2, "rotameric_chi": 1, "terminal_non_rotameric": True, "non_rotameric_bins": 12, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "末端 chi 非 rotameric, 30 deg 细 bin (Table S1 / Simple Mode); 数值表未内嵌 (见模块 docstring)"},
    'ASP': {"chi": 2, "rotameric_chi": 1, "terminal_non_rotameric": True, "non_rotameric_bins": 6, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "末端 chi 非 rotameric, 30 deg 细 bin (Table S1 / Simple Mode); 数值表未内嵌 (见模块 docstring)"},
    'CYS': {"chi": 1, "rotameric_chi": 1, "terminal_non_rotameric": False, "non_rotameric_bins": None, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'GLN': {"chi": 3, "rotameric_chi": 2, "terminal_non_rotameric": True, "non_rotameric_bins": 12, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "末端 chi 非 rotameric, 30 deg 细 bin (Table S1 / Simple Mode); 数值表未内嵌 (见模块 docstring)"},
    'GLU': {"chi": 3, "rotameric_chi": 2, "terminal_non_rotameric": True, "non_rotameric_bins": 6, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "末端 chi 非 rotameric, 30 deg 细 bin (Table S1 / Simple Mode); 数值表未内嵌 (见模块 docstring)"},
    'GLY': {"chi": 0, "rotameric_chi": 0, "terminal_non_rotameric": False, "non_rotameric_bins": None, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'HIS': {"chi": 2, "rotameric_chi": 1, "terminal_non_rotameric": True, "non_rotameric_bins": 12, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "末端 chi 非 rotameric, 30 deg 细 bin (Table S1 / Simple Mode); 数值表未内嵌 (见模块 docstring)"},
    'ILE': {"chi": 2, "rotameric_chi": 2, "terminal_non_rotameric": False, "non_rotameric_bins": None, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'LEU': {"chi": 2, "rotameric_chi": 2, "terminal_non_rotameric": False, "non_rotameric_bins": None, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'LYS': {"chi": 4, "rotameric_chi": 4, "terminal_non_rotameric": False, "non_rotameric_bins": None, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'MET': {"chi": 3, "rotameric_chi": 3, "terminal_non_rotameric": False, "non_rotameric_bins": None, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'PHE': {"chi": 2, "rotameric_chi": 1, "terminal_non_rotameric": True, "non_rotameric_bins": 6, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "末端 chi 非 rotameric, 30 deg 细 bin (Table S1 / Simple Mode); 数值表未内嵌 (见模块 docstring)"},
    'PRO': {"chi": 2, "rotameric_chi": 2, "terminal_non_rotameric": False, "non_rotameric_bins": None, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'SER': {"chi": 1, "rotameric_chi": 1, "terminal_non_rotameric": False, "non_rotameric_bins": None, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'THR': {"chi": 1, "rotameric_chi": 1, "terminal_non_rotameric": False, "non_rotameric_bins": None, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'TRP': {"chi": 2, "rotameric_chi": 1, "terminal_non_rotameric": True, "non_rotameric_bins": 12, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "末端 chi 非 rotameric, 30 deg 细 bin (Table S1 / Simple Mode); 数值表未内嵌 (见模块 docstring)"},
    'TYR': {"chi": 2, "rotameric_chi": 1, "terminal_non_rotameric": True, "non_rotameric_bins": 6, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "末端 chi 非 rotameric, 30 deg 细 bin (Table S1 / Simple Mode); 数值表未内嵌 (见模块 docstring)"},
    'VAL': {"chi": 1, "rotameric_chi": 1, "terminal_non_rotameric": False, "non_rotameric_bins": None, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
}

SIDECHAIN_DIHE_REFS = {
    "rosetta_params_408": "Rosetta 408 ... l-caa/*.params (ICOOR_INTERNAL 规范残基理想几何)",
    "dunbrack_2010": "Shapovalov MV, Dunbrack RL Jr. A smoothed backbone-dependent rotamer library. Structure 19:844-858, 2011."
}
