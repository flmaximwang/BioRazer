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
* ``ROTAMER_BIN`` / ``SIDECHAIN_ROTAMER_LIB`` / ``SIDECHAIN_NON_ROTAMERIC_BINS``
  -- 侧链 rotamer 的**分类框架**: 标准 rotamer bin 中心 (g-/t/g+ =
  -60/180/+60) 和逐残基 rotamer 库 (canonical + 命名 g-/g+/t)。完整的
  **骨架依赖数值表** (逐 phi/psi bin 的均值/方差) 属于外部 Dunbrack 2010
  数据集 (Shapovalov & Dunbrack 2011, CC BY 4.0), 未内嵌, 需要时可另行
  vendor。
* ``NON_ROTAMERIC_BIN_WIDTH`` -- 非 rotameric 末端 chi (ASN/ASP/GLN/GLU/
  PHE/TRP/HIS/TYR 的 sp2 末端) 的 30 deg 细 bin 宽度 (度)。这类 chi 不是
  g-/t/g+ 的离散 rotamer, 而是宽而对称性差的连续分布; 论文用 kernel 密度
  估计建模, 并额外提供 30 deg 离散 bin (每 bin 的均值/方差/占比) 以兼容
  SCWRL 等旧应用。bin 数见 ``SIDECHAIN_NON_ROTAMERIC_BINS`` 的 ``bins``。
* ``SIDECHAIN_ROTAMER_LIB`` -- per-residue rotamer 库。单层键
  ``<RES>_canonical``、``<RES>_g-``、``<RES>_g+``、``<RES>_t`` (2 轴残基为
  ``<RES>_<a>_<b>``) 映射到 {chi_quad: {mean,std,lb,up,source}}。
  ``canonical`` 等于 SIDECHAIN_IC_DIHEDRAL 模板几何 (覆盖 0-deg 状态)。
  命名 rotamer 用 Dunbrack bin 中心。非 rotameric 末端 chi 除 ``canonical``
  外, 另有 ``<RES>_nr<i>`` 条目: 用 ``NON_ROTAMERIC_BIN_WIDTH`` (30 deg)
  将完整周期 (180 或 360 deg) 均匀离散的 ``bins`` 个 bin, 中心在 bin 中点
  (15/45/75/...)。这些是均匀骨架 (source ``dunbrack_2010_uniform_30deg_bin``),
  非 Dunbrack 数据拟合均值。完整骨架依赖数值表未内嵌。
* ``SIDECHAIN_NON_ROTAMERIC_BINS`` -- per-residue 非 rotameric 末端 chi
  (sp3-sp2/芳香) 细 bin 规范: chi 四元组 + 30 deg bin 数 (ASN/GLN/HIS/TRP
  = 12, ASP/GLU/PHE/TYR = 6)。这些 bin 的**中心位置**由
  ``NON_ROTAMERIC_BIN_WIDTH`` 在 ``SIDECHAIN_ROTAMER_LIB`` 的 ``<RES>_nr<i>``
  中生成; 数据拟合的均值/方差 (完整 Dunbrack 数据集许可待定) 尚未内嵌。

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

SIDECHAIN_DIHE_REFS = {
    "rosetta_params_408": "Rosetta 408 ... l-caa/*.params (ICOOR_INTERNAL 规范残基理想几何)",
    "dunbrack_2010": "Shapovalov MV, Dunbrack RL Jr. A smoothed backbone-dependent rotamer library. Structure 19:844-858, 2011.",
    "dunbrack_2010_uniform_30deg_bin": "均匀 30-deg 离散骨架 (bin 中心在 bin 中点, 15/45/75/...), 由 NON_ROTAMERIC_BIN_WIDTH 生成, 兼容 Rosetta backbone-independent bin; 非 Dunbrack 数据拟合均值。"
}

#: 非 rotameric 末端 chi (sp2/芳香末端) 的细 bin 规范: chi 四元组 + 30 deg
#: bin 数 (SCWRL 兼容)。ASN/GLN/HIS/TRP = 12, ASP/GLU/PHE/TYR = 6。这些
#: 末端 chi 不是 g-/g+/t 离散 rotamer, 而是宽而对称性差的连续分布, 用
#: NON_ROTAMERIC_BIN_WIDTH 的细 bin 离散化。bin 的**中心位置**在
#: ``SIDECHAIN_ROTAMER_LIB`` 的 ``<RES>_nr<i>`` 条目中 (均匀 30-deg 骨架);
#: 数据拟合的均值/方差 (完整 Dunbrack 数据集许可待定) 尚未内嵌。
SIDECHAIN_NON_ROTAMERIC_BINS = {
    'ASN': {'chi_quad': ('CA', 'CB', 'CG', 'OD1'), 'bins': 12},
    'ASP': {'chi_quad': ('CA', 'CB', 'CG', 'OD1'), 'bins': 6},
    'GLN': {'chi_quad': ('CB', 'CG', 'CD', 'OE1'),  'bins': 12},
    'GLU': {'chi_quad': ('CB', 'CG', 'CD', 'OE1'),  'bins': 6},
    'HIS': {'chi_quad': ('CA', 'CB', 'CG', 'ND1'),  'bins': 12},
    'PHE': {'chi_quad': ('CA', 'CB', 'CG', 'CD1'),  'bins': 6},
    'TRP': {'chi_quad': ('CA', 'CB', 'CG', 'CD1'),  'bins': 12},
    'TYR': {'chi_quad': ('CA', 'CB', 'CG', 'CD1'),  'bins': 6},
}


def _build_rotamer_lib():
    """Build the per-residue rotamer library programmatically.

    Single-level keys ``<RES>_canonical``, ``<RES>_g-``, ``<RES>_g+``,
    ``<RES>_t`` (and ``<RES>_<a>_<b>`` for 2-axis residues) map to
    {chi_quad: {mean, std, lb, up, source}}.  ``canonical`` equals the
    SIDECHAIN_IC_DIHEDRAL template geometry (covers the 0-deg state);
    named rotamers use Dunbrack bin centers.

    Non-rotameric terminal chi (members of SIDECHAIN_NON_ROTAMERIC_BINS,
    e.g. the sp2 terminal chi of ASN/ASP/GLN/GLU/PHE/HIS/TRP/TYR) appear
    in ``canonical`` AND as ``<RES>_nr<i>`` entries: ``bins`` uniform
    NON_ROTAMERIC_BIN_WIDTH (30 deg) bins covering the full period
    (bins*30 = 180 deg for the 6-bin residues, 360 deg for the 12-bin
    residues).  bin centers sit at the bin midpoint, starting at
    width/2 = 15 deg (15, 45, 75, ...), matching Rosetta's
    backbone-independent 30-deg discretisation.  These centers are a
    uniform scaffold, NOT the data-fitted means of the full Dunbrack
    tables (those are not embedded); hence their source is
    ``dunbrack_2010_uniform_30deg_bin``.

    Named rotamers only cover the first 2 chi (chi1 / chi1 x chi2); chi3+
    stay canonical.  PRO is ring-constrained (its chi cannot rotate into
    g-/g+/t rotamers), so it has no named rotamers.
    """
    labels = ["g-", "g+", "t"]
    centers = {"g-": -60.0, "g+": 60.0, "t": 180.0}
    #: 环约束残基: chi 在环内, 不构成可旋转 rotamer (只有 canonical)。
    RING_CONSTRAINED = {"PRO"}
    #: 均匀 30-deg 非 rotameric bin 的 source 标签 (均匀骨架, 非 Dunbrack 拟合均值)。
    NR_BIN_SOURCE = "dunbrack_2010_uniform_30deg_bin"
    lib = {}
    for res in AAS:
        chis = SIDECHAIN_CHI[res]
        ic = SIDECHAIN_IC_DIHEDRAL[res]
        def rec(mean, source):
            return {"mean": float(mean), "std": np.nan, "lb": np.nan,
                    "up": np.nan, "source": source}
        # canonical: all chi present in IC_DIHEDRAL -> template mean
        canon = {q: rec(ic[q]["mean"], "rosetta_params_408")
                 for q in chis if q in ic}
        lib[f"{res}_canonical"] = canon
        # rotameric chi count for naming; named rotamers cover only the
        # first 2 chi (chi1 x chi2); PRO ring chi are not rotamers.
        if res in RING_CONSTRAINED:
            continue
        rc = len(chis) - (1 if res in SIDECHAIN_NON_ROTAMERIC_BINS else 0)
        rc = min(2, rc)
        if rc == 0:
            pass  # no rotameric chi to rotate (ALA/GLY)
        elif rc == 1:
            lib[f"{res}_g-"] = {chis[0]: rec(centers["g-"], "dunbrack_2010")}
            lib[f"{res}_g+"] = {chis[0]: rec(centers["g+"], "dunbrack_2010")}
            lib[f"{res}_t"]  = {chis[0]: rec(centers["t"],  "dunbrack_2010")}
        else:  # rc == 2
            for a in labels:
                for b in labels:
                    lib[f"{res}_{a}_{b}"] = {
                        chis[0]: rec(centers[a], "dunbrack_2010"),
                        chis[1]: rec(centers[b], "dunbrack_2010"),
                    }
        # non-rotameric terminal chi: bins uniform 30-deg bins over the
        # full period.  The terminal chi is the last entry of chis.
        if res in SIDECHAIN_NON_ROTAMERIC_BINS:
            spec = SIDECHAIN_NON_ROTAMERIC_BINS[res]
            nbin = spec["bins"]
            nq = spec["chi_quad"]
            for i in range(1, nbin + 1):
                center = NON_ROTAMERIC_BIN_WIDTH * (i - 0.5)
                lib[f"{res}_nr{i}"] = {nq: rec(center, NR_BIN_SOURCE)}
    return lib

#: 每种残基的逐残基 rotamer 库: {chi_quad: {mean, std, lb, up, source}}。
#: ``canonical`` 是头等 rotamer (quad->mean == SIDECHAIN_IC_DIHEDRAL 模板值,
#: 覆盖 mean=0 状态)。命名 rotamer (g-/g+/t) 用 Dunbrack bin 中心, 只覆盖
#: 前导可 rotameric chi (chi1 / chi1+chi2)。非 rotameric 末端 chi 除
#: ``canonical`` 外, 还有 ``<RES>_nr<i>`` 条目 (i=1..bins): 用
#: NON_ROTAMERIC_BIN_WIDTH 均匀离散完整周期的 bin, 中心在 bin 中点
#: (15/45/75/...), 为均匀骨架 (非 Dunbrack 拟合均值)。完整骨架依赖数值表
#: 未内嵌。
#:
#: 键形: 单 chi 残基 -> ``<RES>_{canonical,g-,g+,t}``; 双 chi 残基 ->
#: ``<RES>_{g-,g+,t}_{g-,g+,t}`` (仅 chi1 x chi2); 非 rotameric 末端 chi ->
#: ``<RES>_nr1..<RES>_nr{bins}``。值 dict 的键是 chi 扭转角四元组 (原子名),
#: 值各字段含义:
#:   ``mean``   -- bin 中心二面角 (度)
#:   ``std``    -- 标准差 (未内嵌时为 np.nan)
#:   ``lb``/``up`` -- bin 下/上界 (未内嵌时为 np.nan)
#:   ``source`` -- 数据来源 ("rosetta_params_408" / "dunbrack_2010" /
#:                 "dunbrack_2010_uniform_30deg_bin")
#:
#: 示例 (真实内容):
#:   >>> SIDECHAIN_ROTAMER_LIB['VAL_t']          # 单 chi, t (trans)
#:   {('N', 'CA', 'CB', 'CG1'): {'mean': 180.0, 'std': nan,
#:                               'lb': nan, 'up': nan,
#:                               'source': 'dunbrack_2010'}}
#:   >>> SIDECHAIN_ROTAMER_LIB['LEU_t_g-']       # 双 chi: chi1=t, chi2=g-
#:   {('N', 'CA', 'CB', 'CG'): {'mean': 180.0, ...},
#:    ('CA', 'CB', 'CG', 'CD1'): {'mean': -60.0, ...}}
#:   >>> SIDECHAIN_ROTAMER_LIB['SER_canonical']  # canonical == IC 模板几何
#:   {('N', 'CA', 'CB', 'OG'): {'mean': -0.0, 'source': 'rosetta_params_408', ...}}
#:   >>> SIDECHAIN_ROTAMER_LIB['ASP_nr1']        # 非 rotameric chi2 首 bin
#:   {('CA', 'CB', 'CG', 'OD1'): {'mean': 15.0, 'std': nan,
#:                                'lb': nan, 'up': nan,
#:                                'source': 'dunbrack_2010_uniform_30deg_bin'}}
SIDECHAIN_ROTAMER_LIB = _build_rotamer_lib()
