# -*- coding: utf-8 -*-
"""Amino-acid side-chain covalent bond lengths (Å) and bond angles (°).

数据集来源与诚实性说明
───────────────────────
本模块记录 20 种标准氨基酸**侧链**(CA 的 CB 及以外重原子, 见
`internal_coord_template._topology.SIDE_CHAIN_IC_PATH` 的生长路径)的理想键长与键角。
数值为**理想点值 (ideal point values)**, 翻译自 Rosetta 408 的
`fa_standard` 残基 .params (`ICOOR_INTERNAL`) 构建的规范残基几何 ——
Rosetta 的主链/侧链理想几何源自 Engh & Huber (1991) (见
`bond.backbone` 的交叉核对)。键长即 ICOOR 的 `d`, 键角即 `180° − theta`
(解码后又在构建出的残基坐标上实测复核)。

注意: 这里**没有**给出 Engh & Huber 的样本标准差 (std) —— Rosetta ICOOR
只给理想点值, 不给 CSD sigma; 需要 std 时请以 E&H 1991 原始文献为准。
键长单位 Å, 键角单位 **度 (degree)**。

键/角的 key 采用 `SIDE_CHAIN_IC_PATH` 的生长四元组 (i,j,k,l): 从父原子
(i,j,k) 生长出 l, 键 `(k,l)`, 键角 `(j,k,l)` (顶点在 k)。
"""


#: 每种氨基酸侧链生长的规范键长 (k,l) 与键角 (j,k,l)
#: {(res_name): {(i,j,k,l)取原子名的三元/四元: {...}}}
#: 含第一个 CB 生长四元组 (N,C,CA,CB) (其 CA-CB 键/角度亦见 bond.backbone)
AMINO_ACID_SIDECHAIN_BOND = {
    "ALA": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5217, "mean": 1.5217, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.37, "note": "Rosetta fa_standard ALA.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "ARG": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5216, "mean": 1.5216, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.60, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.5204, "mean": 1.5204, "angle_pair": ('CA', 'CB', 'CG'), "angle": 114.87, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'CD'): {"bond": ('CG', 'CD'), "value": 1.4854, "mean": 1.4854, "angle_pair": ('CB', 'CG', 'CD'), "angle": 111.70, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD', 'NE'): {"bond": ('CD', 'NE'), "value": 1.4541, "mean": 1.4541, "angle_pair": ('CG', 'CD', 'NE'), "angle": 111.90, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CD', 'NE', 'CZ'): {"bond": ('NE', 'CZ'), "value": 1.3473, "mean": 1.3473, "angle_pair": ('CD', 'NE', 'CZ'), "angle": 124.60, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD', 'NE', 'CZ', 'NH1'): {"bond": ('CZ', 'NH1'), "value": 1.3146, "mean": 1.3146, "angle_pair": ('NE', 'CZ', 'NH1'), "angle": 120.00, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
        ('NH1', 'NE', 'CZ', 'NH2'): {"bond": ('CZ', 'NH2'), "value": 1.3216, "mean": 1.3216, "angle_pair": ('NE', 'CZ', 'NH2'), "angle": 120.00, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "ASN": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5177, "mean": 1.5177, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.60, "note": "Rosetta fa_standard ASN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.5035, "mean": 1.5035, "angle_pair": ('CA', 'CB', 'CG'), "angle": 112.60, "note": "Rosetta fa_standard ASN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'OD1'): {"bond": ('CG', 'OD1'), "value": 1.2364, "mean": 1.2364, "angle_pair": ('CB', 'CG', 'OD1'), "angle": 120.80, "note": "Rosetta fa_standard ASN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('OD1', 'CB', 'CG', 'ND2'): {"bond": ('CG', 'ND2'), "value": 1.3086, "mean": 1.3086, "angle_pair": ('CB', 'CG', 'ND2'), "angle": 116.50, "note": "Rosetta fa_standard ASN.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "ASP": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5307, "mean": 1.5307, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.52, "note": "Rosetta fa_standard ASP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.5228, "mean": 1.5228, "angle_pair": ('CA', 'CB', 'CG'), "angle": 112.90, "note": "Rosetta fa_standard ASP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'OD1'): {"bond": ('CG', 'OD1'), "value": 1.2082, "mean": 1.2082, "angle_pair": ('CB', 'CG', 'OD1'), "angle": 118.39, "note": "Rosetta fa_standard ASP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('OD1', 'CB', 'CG', 'OD2'): {"bond": ('CG', 'OD2'), "value": 1.2078, "mean": 1.2078, "angle_pair": ('CB', 'CG', 'OD2'), "angle": 118.37, "note": "Rosetta fa_standard ASP.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "CYS": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5289, "mean": 1.5289, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.60, "note": "Rosetta fa_standard CYS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'SG'): {"bond": ('CB', 'SG'), "value": 1.8088, "mean": 1.8088, "angle_pair": ('CA', 'CB', 'SG'), "angle": 114.10, "note": "Rosetta fa_standard CYS.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "GLN": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5311, "mean": 1.5311, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.49, "note": "Rosetta fa_standard GLN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.5191, "mean": 1.5191, "angle_pair": ('CA', 'CB', 'CG'), "angle": 114.19, "note": "Rosetta fa_standard GLN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'CD'): {"bond": ('CG', 'CD'), "value": 1.5169, "mean": 1.5169, "angle_pair": ('CB', 'CG', 'CD'), "angle": 112.43, "note": "Rosetta fa_standard GLN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD', 'OE1'): {"bond": ('CD', 'OE1'), "value": 1.2342, "mean": 1.2342, "angle_pair": ('CG', 'CD', 'OE1'), "angle": 120.95, "note": "Rosetta fa_standard GLN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('OE1', 'CG', 'CD', 'NE2'): {"bond": ('CD', 'NE2'), "value": 1.3281, "mean": 1.3281, "angle_pair": ('CG', 'CD', 'NE2'), "angle": 116.41, "note": "Rosetta fa_standard GLN.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "GLU": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5303, "mean": 1.5303, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.40, "note": "Rosetta fa_standard GLU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.5221, "mean": 1.5221, "angle_pair": ('CA', 'CB', 'CG'), "angle": 114.40, "note": "Rosetta fa_standard GLU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'CD'): {"bond": ('CG', 'CD'), "value": 1.5034, "mean": 1.5034, "angle_pair": ('CB', 'CG', 'CD'), "angle": 112.90, "note": "Rosetta fa_standard GLU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD', 'OE1'): {"bond": ('CD', 'OE1'), "value": 1.2076, "mean": 1.2076, "angle_pair": ('CG', 'CD', 'OE1'), "angle": 118.45, "note": "Rosetta fa_standard GLU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('OE1', 'CG', 'CD', 'OE2'): {"bond": ('CD', 'OE2'), "value": 1.2085, "mean": 1.2085, "angle_pair": ('CG', 'CD', 'OE2'), "angle": 118.36, "note": "Rosetta fa_standard GLU.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "GLY": {
    },
    "HIS": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5321, "mean": 1.5321, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.69, "note": "Rosetta fa_standard HIS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.4972, "mean": 1.4972, "angle_pair": ('CA', 'CB', 'CG'), "angle": 113.69, "note": "Rosetta fa_standard HIS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'ND1'): {"bond": ('CG', 'ND1'), "value": 1.3792, "mean": 1.3792, "angle_pair": ('CB', 'CG', 'ND1'), "angle": 122.63, "note": "Rosetta fa_standard HIS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'ND1', 'CE1'): {"bond": ('ND1', 'CE1'), "value": 1.3219, "mean": 1.3219, "angle_pair": ('CG', 'ND1', 'CE1'), "angle": 109.29, "note": "Rosetta fa_standard HIS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'ND1', 'CE1', 'NE2'): {"bond": ('CE1', 'NE2'), "value": 1.3204, "mean": 1.3204, "angle_pair": ('ND1', 'CE1', 'NE2'), "angle": 108.38, "note": "Rosetta fa_standard HIS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('ND1', 'CE1', 'NE2', 'CD2'): {"bond": ('NE2', 'CD2'), "value": 1.3732, "mean": 1.3732, "angle_pair": ('CE1', 'NE2', 'CD2'), "angle": 109.01, "note": "Rosetta fa_standard HIS.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "ILE": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5396, "mean": 1.5396, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.90, "note": "Rosetta fa_standard ILE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG1'): {"bond": ('CB', 'CG1'), "value": 1.5309, "mean": 1.5309, "angle_pair": ('CA', 'CB', 'CG1'), "angle": 110.41, "note": "Rosetta fa_standard ILE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG1', 'CD1'): {"bond": ('CG1', 'CD1'), "value": 1.5117, "mean": 1.5117, "angle_pair": ('CB', 'CG1', 'CD1'), "angle": 113.83, "note": "Rosetta fa_standard ILE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG1', 'CA', 'CB', 'CG2'): {"bond": ('CB', 'CG2'), "value": 1.5209, "mean": 1.5209, "angle_pair": ('CA', 'CB', 'CG2'), "angle": 110.47, "note": "Rosetta fa_standard ILE.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "LEU": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5339, "mean": 1.5339, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.21, "note": "Rosetta fa_standard LEU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.5340, "mean": 1.5340, "angle_pair": ('CA', 'CB', 'CG'), "angle": 115.70, "note": "Rosetta fa_standard LEU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'CD1'): {"bond": ('CG', 'CD1'), "value": 1.5227, "mean": 1.5227, "angle_pair": ('CB', 'CG', 'CD1'), "angle": 109.50, "note": "Rosetta fa_standard LEU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD1', 'CB', 'CG', 'CD2'): {"bond": ('CG', 'CD2'), "value": 1.5214, "mean": 1.5214, "angle_pair": ('CB', 'CG', 'CD2'), "angle": 109.50, "note": "Rosetta fa_standard LEU.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "LYS": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5295, "mean": 1.5295, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.50, "note": "Rosetta fa_standard LYS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.5229, "mean": 1.5229, "angle_pair": ('CA', 'CB', 'CG'), "angle": 114.40, "note": "Rosetta fa_standard LYS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'CD'): {"bond": ('CG', 'CD'), "value": 1.5213, "mean": 1.5213, "angle_pair": ('CB', 'CG', 'CD'), "angle": 111.30, "note": "Rosetta fa_standard LYS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD', 'CE'): {"bond": ('CD', 'CE'), "value": 1.5216, "mean": 1.5216, "angle_pair": ('CG', 'CD', 'CE'), "angle": 111.37, "note": "Rosetta fa_standard LYS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CD', 'CE', 'NZ'): {"bond": ('CE', 'NZ'), "value": 1.4881, "mean": 1.4881, "angle_pair": ('CD', 'CE', 'NZ'), "angle": 111.96, "note": "Rosetta fa_standard LYS.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "MET": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5274, "mean": 1.5274, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.21, "note": "Rosetta fa_standard MET.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.5222, "mean": 1.5222, "angle_pair": ('CA', 'CB', 'CG'), "angle": 114.44, "note": "Rosetta fa_standard MET.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'SD'): {"bond": ('CG', 'SD'), "value": 1.8038, "mean": 1.8038, "angle_pair": ('CB', 'CG', 'SD'), "angle": 112.67, "note": "Rosetta fa_standard MET.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'SD', 'CE'): {"bond": ('SD', 'CE'), "value": 1.7904, "mean": 1.7904, "angle_pair": ('CG', 'SD', 'CE'), "angle": 100.89, "note": "Rosetta fa_standard MET.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "PHE": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5298, "mean": 1.5298, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.44, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.5022, "mean": 1.5022, "angle_pair": ('CA', 'CB', 'CG'), "angle": 113.79, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'CD1'): {"bond": ('CG', 'CD1'), "value": 1.3870, "mean": 1.3870, "angle_pair": ('CB', 'CG', 'CD1'), "angle": 120.69, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD1', 'CE1'): {"bond": ('CD1', 'CE1'), "value": 1.3822, "mean": 1.3822, "angle_pair": ('CG', 'CD1', 'CE1'), "angle": 120.73, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CD1', 'CE1', 'CZ'): {"bond": ('CE1', 'CZ'), "value": 1.3786, "mean": 1.3786, "angle_pair": ('CD1', 'CE1', 'CZ'), "angle": 120.03, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD1', 'CE1', 'CZ', 'CE2'): {"bond": ('CZ', 'CE2'), "value": 1.3805, "mean": 1.3805, "angle_pair": ('CE1', 'CZ', 'CE2'), "angle": 119.87, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CE1', 'CZ', 'CE2', 'CD2'): {"bond": ('CE2', 'CD2'), "value": 1.3813, "mean": 1.3813, "angle_pair": ('CZ', 'CE2', 'CD2'), "angle": 119.97, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "PRO": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5320, "mean": 1.5320, "angle_pair": ('N', 'CA', 'CB'), "angle": 103.00, "note": "Rosetta fa_standard PRO.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.4906, "mean": 1.4906, "angle_pair": ('CA', 'CB', 'CG'), "angle": 104.20, "note": "Rosetta fa_standard PRO.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'CD'): {"bond": ('CG', 'CD'), "value": 1.5055, "mean": 1.5055, "angle_pair": ('CB', 'CG', 'CD'), "angle": 104.40, "note": "Rosetta fa_standard PRO.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "SER": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5163, "mean": 1.5163, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.14, "note": "Rosetta fa_standard SER.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'OG'): {"bond": ('CB', 'OG'), "value": 1.4012, "mean": 1.4012, "angle_pair": ('CA', 'CB', 'OG'), "angle": 111.00, "note": "Rosetta fa_standard SER.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "THR": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5399, "mean": 1.5399, "angle_pair": ('N', 'CA', 'CB'), "angle": 111.10, "note": "Rosetta fa_standard THR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'OG1'): {"bond": ('CB', 'OG1'), "value": 1.4335, "mean": 1.4335, "angle_pair": ('CA', 'CB', 'OG1'), "angle": 109.58, "note": "Rosetta fa_standard THR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('OG1', 'CA', 'CB', 'CG2'): {"bond": ('CB', 'CG2'), "value": 1.5210, "mean": 1.5210, "angle_pair": ('CA', 'CB', 'CG2'), "angle": 110.53, "note": "Rosetta fa_standard THR.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "TRP": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5298, "mean": 1.5298, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.40, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.4987, "mean": 1.4987, "angle_pair": ('CA', 'CB', 'CG'), "angle": 113.53, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'CD1'): {"bond": ('CG', 'CD1'), "value": 1.3627, "mean": 1.3627, "angle_pair": ('CB', 'CG', 'CD1'), "angle": 126.74, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD1', 'NE1'): {"bond": ('CD1', 'NE1'), "value": 1.3729, "mean": 1.3729, "angle_pair": ('CG', 'CD1', 'NE1'), "angle": 110.16, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CD1', 'NE1', 'CE2'): {"bond": ('NE1', 'CE2'), "value": 1.3721, "mean": 1.3721, "angle_pair": ('CD1', 'NE1', 'CE2'), "angle": 108.90, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD1', 'NE1', 'CE2', 'CZ2'): {"bond": ('CE2', 'CZ2'), "value": 1.3859, "mean": 1.3859, "angle_pair": ('NE1', 'CE2', 'CZ2'), "angle": 130.14, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('NE1', 'CE2', 'CZ2', 'CH2'): {"bond": ('CZ2', 'CH2'), "value": 1.3950, "mean": 1.3950, "angle_pair": ('CE2', 'CZ2', 'CH2'), "angle": 117.50, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CE2', 'CZ2', 'CH2', 'CZ3'): {"bond": ('CH2', 'CZ3'), "value": 1.3721, "mean": 1.3721, "angle_pair": ('CZ2', 'CH2', 'CZ3'), "angle": 121.50, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CZ2', 'CH2', 'CZ3', 'CE3'): {"bond": ('CZ3', 'CE3'), "value": 1.3898, "mean": 1.3898, "angle_pair": ('CH2', 'CZ3', 'CE3'), "angle": 121.00, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CH2', 'CZ3', 'CE3', 'CD2'): {"bond": ('CE3', 'CD2'), "value": 1.4004, "mean": 1.4004, "angle_pair": ('CZ3', 'CE3', 'CD2'), "angle": 118.73, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "TYR": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5304, "mean": 1.5304, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.47, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.5127, "mean": 1.5127, "angle_pair": ('CA', 'CB', 'CG'), "angle": 113.80, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'CD1'): {"bond": ('CG', 'CD1'), "value": 1.3872, "mean": 1.3872, "angle_pair": ('CB', 'CG', 'CD1'), "angle": 120.00, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD1', 'CE1'): {"bond": ('CD1', 'CE1'), "value": 1.3816, "mean": 1.3816, "angle_pair": ('CG', 'CD1', 'CE1'), "angle": 120.00, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD1', 'CB', 'CG', 'CD2'): {"bond": ('CG', 'CD2'), "value": 1.3869, "mean": 1.3869, "angle_pair": ('CB', 'CG', 'CD2'), "angle": 120.00, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD2', 'CE2'): {"bond": ('CD2', 'CE2'), "value": 1.3814, "mean": 1.3814, "angle_pair": ('CG', 'CD2', 'CE2'), "angle": 120.00, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CD2', 'CE2', 'CZ'): {"bond": ('CE2', 'CZ'), "value": 1.3800, "mean": 1.3800, "angle_pair": ('CD2', 'CE2', 'CZ'), "angle": 120.00, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD2', 'CE2', 'CZ', 'OH'): {"bond": ('CZ', 'OH'), "value": 1.3760, "mean": 1.3760, "angle_pair": ('CE2', 'CZ', 'OH'), "angle": 120.00, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "VAL": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5402, "mean": 1.5402, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.30, "note": "Rosetta fa_standard VAL.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG1'): {"bond": ('CB', 'CG1'), "value": 1.5214, "mean": 1.5214, "angle_pair": ('CA', 'CB', 'CG1'), "angle": 110.49, "note": "Rosetta fa_standard VAL.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG1', 'CA', 'CB', 'CG2'): {"bond": ('CB', 'CG2'), "value": 1.5211, "mean": 1.5211, "angle_pair": ('CA', 'CB', 'CG2'), "angle": 109.90, "note": "Rosetta fa_standard VAL.params ICOOR ideal", "source": "rosetta_params_408"},
    },
}

BOND_SIDECHAIN_REFS = {
    "rosetta_params_408": "Rosetta 408 main/database/chemical/residue_type_sets/fa_standard/residue_types/l-caa/*.params (ICOOR_INTERNAL 理想键长/键角)。",
}