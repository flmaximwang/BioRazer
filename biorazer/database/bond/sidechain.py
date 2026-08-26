# -*- coding: utf-8 -*-
"""Amino-acid side-chain covalent bond lengths (Å) and bond angles (°).

数据集来源与诚实性说明
───────────────────────
本模块记录 20 种标准氨基酸**侧链**(CA 的 CB 及以外重原子, 见
`amino_acid_internal_coords.SIDE_CHAIN_IC_PATH` 的生长路径)的理想键长与键角。
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
        ('N', 'C', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5217, "mean": 1.5217, "angle_pair": ('C', 'CA', 'CB'), "angle": 110.35, "note": "Rosetta fa_standard ALA.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "ARG": {
        ('N', 'C', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5216, "mean": 1.5216, "angle_pair": ('C', 'CA', 'CB'), "angle": 109.37, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'N', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.5204, "mean": 1.5204, "angle_pair": ('N', 'CB', 'CG'), "angle": 81.02, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CA', 'CG', 'CD'): {"bond": ('CG', 'CD'), "value": 1.4854, "mean": 1.4854, "angle_pair": ('CA', 'CG', 'CD'), "angle": 79.12, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CB', 'CD', 'NE'): {"bond": ('CD', 'NE'), "value": 1.4541, "mean": 1.4541, "angle_pair": ('CB', 'CD', 'NE'), "angle": 77.30, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD', 'CG', 'NE', 'CZ'): {"bond": ('NE', 'CZ'), "value": 1.3473, "mean": 1.3473, "angle_pair": ('CG', 'NE', 'CZ'), "angle": 90.14, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
        ('NE', 'CD', 'CZ', 'NH1'): {"bond": ('CZ', 'NH1'), "value": 1.3146, "mean": 1.3146, "angle_pair": ('CD', 'CZ', 'NH1'), "angle": 91.15, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
        ('NE', 'NH1', 'CZ', 'NH2'): {"bond": ('CZ', 'NH2'), "value": 1.3216, "mean": 1.3216, "angle_pair": ('NH1', 'CZ', 'NH2'), "angle": 120.00, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "ASN": {
        ('N', 'C', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5177, "mean": 1.5177, "angle_pair": ('C', 'CA', 'CB'), "angle": 110.62, "note": "Rosetta fa_standard ASN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'N', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.5035, "mean": 1.5035, "angle_pair": ('N', 'CB', 'CG'), "angle": 78.70, "note": "Rosetta fa_standard ASN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CA', 'CG', 'OD1'): {"bond": ('CG', 'OD1'), "value": 1.2364, "mean": 1.2364, "angle_pair": ('CA', 'CG', 'OD1'), "angle": 86.92, "note": "Rosetta fa_standard ASN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'OD1', 'CG', 'ND2'): {"bond": ('CG', 'ND2'), "value": 1.3086, "mean": 1.3086, "angle_pair": ('OD1', 'CG', 'ND2'), "angle": 122.70, "note": "Rosetta fa_standard ASN.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "ASP": {
        ('N', 'C', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5307, "mean": 1.5307, "angle_pair": ('C', 'CA', 'CB'), "angle": 109.72, "note": "Rosetta fa_standard ASP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'N', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.5228, "mean": 1.5228, "angle_pair": ('N', 'CB', 'CG'), "angle": 79.13, "note": "Rosetta fa_standard ASP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CA', 'CG', 'OD1'): {"bond": ('CG', 'OD1'), "value": 1.2082, "mean": 1.2082, "angle_pair": ('CA', 'CG', 'OD1'), "angle": 84.75, "note": "Rosetta fa_standard ASP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'OD1', 'CG', 'OD2'): {"bond": ('CG', 'OD2'), "value": 1.2078, "mean": 1.2078, "angle_pair": ('OD1', 'CG', 'OD2'), "angle": 123.24, "note": "Rosetta fa_standard ASP.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "CYS": {
        ('N', 'C', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5289, "mean": 1.5289, "angle_pair": ('C', 'CA', 'CB'), "angle": 109.29, "note": "Rosetta fa_standard CYS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'N', 'CB', 'SG'): {"bond": ('CB', 'SG'), "value": 1.8088, "mean": 1.8088, "angle_pair": ('N', 'CB', 'SG'), "angle": 80.34, "note": "Rosetta fa_standard CYS.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "GLN": {
        ('N', 'C', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5311, "mean": 1.5311, "angle_pair": ('C', 'CA', 'CB'), "angle": 109.97, "note": "Rosetta fa_standard GLN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'N', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.5191, "mean": 1.5191, "angle_pair": ('N', 'CB', 'CG'), "angle": 80.41, "note": "Rosetta fa_standard GLN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CA', 'CG', 'CD'): {"bond": ('CG', 'CD'), "value": 1.5169, "mean": 1.5169, "angle_pair": ('CA', 'CG', 'CD'), "angle": 79.38, "note": "Rosetta fa_standard GLN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CB', 'CD', 'OE1'): {"bond": ('CD', 'OE1'), "value": 1.2342, "mean": 1.2342, "angle_pair": ('CB', 'CD', 'OE1'), "angle": 87.14, "note": "Rosetta fa_standard GLN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'OE1', 'CD', 'NE2'): {"bond": ('CD', 'NE2'), "value": 1.3281, "mean": 1.3281, "angle_pair": ('OE1', 'CD', 'NE2'), "angle": 122.64, "note": "Rosetta fa_standard GLN.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "GLU": {
        ('N', 'C', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5303, "mean": 1.5303, "angle_pair": ('C', 'CA', 'CB'), "angle": 109.87, "note": "Rosetta fa_standard GLU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'N', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.5221, "mean": 1.5221, "angle_pair": ('N', 'CB', 'CG'), "angle": 80.56, "note": "Rosetta fa_standard GLU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CA', 'CG', 'CD'): {"bond": ('CG', 'CD'), "value": 1.5034, "mean": 1.5034, "angle_pair": ('CA', 'CG', 'CD'), "angle": 80.00, "note": "Rosetta fa_standard GLU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CB', 'CD', 'OE1'): {"bond": ('CD', 'OE1'), "value": 1.2076, "mean": 1.2076, "angle_pair": ('CB', 'CD', 'OE1'), "angle": 84.66, "note": "Rosetta fa_standard GLU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'OE1', 'CD', 'OE2'): {"bond": ('CD', 'OE2'), "value": 1.2085, "mean": 1.2085, "angle_pair": ('OE1', 'CD', 'OE2'), "angle": 123.19, "note": "Rosetta fa_standard GLU.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "GLY": {
    },
    "HIS": {
        ('N', 'C', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5321, "mean": 1.5321, "angle_pair": ('C', 'CA', 'CB'), "angle": 109.71, "note": "Rosetta fa_standard HIS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'N', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.4972, "mean": 1.4972, "angle_pair": ('N', 'CB', 'CG'), "angle": 80.02, "note": "Rosetta fa_standard HIS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CA', 'CG', 'ND1'): {"bond": ('CG', 'ND1'), "value": 1.3792, "mean": 1.3792, "angle_pair": ('CA', 'CG', 'ND1'), "angle": 89.04, "note": "Rosetta fa_standard HIS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CB', 'ND1', 'CE1'): {"bond": ('ND1', 'CE1'), "value": 1.3219, "mean": 1.3219, "angle_pair": ('CB', 'ND1', 'CE1'), "angle": 139.26, "note": "Rosetta fa_standard HIS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('ND1', 'CG', 'CE1', 'NE2'): {"bond": ('CE1', 'NE2'), "value": 1.3204, "mean": 1.3204, "angle_pair": ('CG', 'CE1', 'NE2'), "angle": 72.16, "note": "Rosetta fa_standard HIS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CE1', 'ND1', 'NE2', 'CD2'): {"bond": ('NE2', 'CD2'), "value": 1.3732, "mean": 1.3732, "angle_pair": ('ND1', 'NE2', 'CD2'), "angle": 73.17, "note": "Rosetta fa_standard HIS.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "ILE": {
        ('N', 'C', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5396, "mean": 1.5396, "angle_pair": ('C', 'CA', 'CB'), "angle": 109.52, "note": "Rosetta fa_standard ILE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'N', 'CB', 'CG1'): {"bond": ('CB', 'CG1'), "value": 1.5309, "mean": 1.5309, "angle_pair": ('N', 'CB', 'CG1'), "angle": 76.93, "note": "Rosetta fa_standard ILE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CA', 'CG1', 'CD1'): {"bond": ('CG1', 'CD1'), "value": 1.5117, "mean": 1.5117, "angle_pair": ('CA', 'CG1', 'CD1'), "angle": 78.92, "note": "Rosetta fa_standard ILE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CG1', 'CB', 'CG2'): {"bond": ('CB', 'CG2'), "value": 1.5209, "mean": 1.5209, "angle_pair": ('CG1', 'CB', 'CG2'), "angle": 110.62, "note": "Rosetta fa_standard ILE.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "LEU": {
        ('N', 'C', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5339, "mean": 1.5339, "angle_pair": ('C', 'CA', 'CB'), "angle": 109.73, "note": "Rosetta fa_standard LEU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'N', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.5340, "mean": 1.5340, "angle_pair": ('N', 'CB', 'CG'), "angle": 81.82, "note": "Rosetta fa_standard LEU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CA', 'CG', 'CD1'): {"bond": ('CG', 'CD1'), "value": 1.5227, "mean": 1.5227, "angle_pair": ('CA', 'CG', 'CD1'), "angle": 77.35, "note": "Rosetta fa_standard LEU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CD1', 'CG', 'CD2'): {"bond": ('CG', 'CD2'), "value": 1.5214, "mean": 1.5214, "angle_pair": ('CD1', 'CG', 'CD2'), "angle": 111.06, "note": "Rosetta fa_standard LEU.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "LYS": {
        ('N', 'C', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5295, "mean": 1.5295, "angle_pair": ('C', 'CA', 'CB'), "angle": 110.13, "note": "Rosetta fa_standard LYS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'N', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.5229, "mean": 1.5229, "angle_pair": ('N', 'CB', 'CG'), "angle": 80.60, "note": "Rosetta fa_standard LYS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CA', 'CG', 'CD'): {"bond": ('CG', 'CD'), "value": 1.5213, "mean": 1.5213, "angle_pair": ('CA', 'CG', 'CD'), "angle": 78.42, "note": "Rosetta fa_standard LYS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CB', 'CD', 'CE'): {"bond": ('CD', 'CE'), "value": 1.5216, "mean": 1.5216, "angle_pair": ('CB', 'CD', 'CE'), "angle": 77.00, "note": "Rosetta fa_standard LYS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD', 'CG', 'CE', 'NZ'): {"bond": ('CE', 'NZ'), "value": 1.4881, "mean": 1.4881, "angle_pair": ('CG', 'CE', 'NZ'), "angle": 77.65, "note": "Rosetta fa_standard LYS.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "MET": {
        ('N', 'C', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5274, "mean": 1.5274, "angle_pair": ('C', 'CA', 'CB'), "angle": 110.44, "note": "Rosetta fa_standard MET.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'N', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.5222, "mean": 1.5222, "angle_pair": ('N', 'CB', 'CG'), "angle": 80.47, "note": "Rosetta fa_standard MET.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CA', 'CG', 'SD'): {"bond": ('CG', 'SD'), "value": 1.8038, "mean": 1.8038, "angle_pair": ('CA', 'CG', 'SD'), "angle": 79.83, "note": "Rosetta fa_standard MET.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CB', 'SD', 'CE'): {"bond": ('SD', 'CE'), "value": 1.7904, "mean": 1.7904, "angle_pair": ('CB', 'SD', 'CE'), "angle": 70.46, "note": "Rosetta fa_standard MET.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "PHE": {
        ('N', 'C', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5298, "mean": 1.5298, "angle_pair": ('C', 'CA', 'CB'), "angle": 109.69, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'N', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.5022, "mean": 1.5022, "angle_pair": ('N', 'CB', 'CG'), "angle": 79.97, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CA', 'CG', 'CD1'): {"bond": ('CG', 'CD1'), "value": 1.3870, "mean": 1.3870, "angle_pair": ('CA', 'CG', 'CD1'), "angle": 87.24, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CB', 'CD1', 'CE1'): {"bond": ('CD1', 'CE1'), "value": 1.3822, "mean": 1.3822, "angle_pair": ('CB', 'CD1', 'CE1'), "angle": 151.69, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD1', 'CG', 'CE1', 'CZ'): {"bond": ('CE1', 'CZ'), "value": 1.3786, "mean": 1.3786, "angle_pair": ('CG', 'CE1', 'CZ'), "angle": 90.34, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CE1', 'CD1', 'CZ', 'CE2'): {"bond": ('CZ', 'CE2'), "value": 1.3805, "mean": 1.3805, "angle_pair": ('CD1', 'CZ', 'CE2'), "angle": 89.84, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CZ', 'CE1', 'CE2', 'CD2'): {"bond": ('CE2', 'CD2'), "value": 1.3813, "mean": 1.3813, "angle_pair": ('CE1', 'CE2', 'CD2'), "angle": 89.93, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "PRO": {
        ('N', 'C', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5320, "mean": 1.5320, "angle_pair": ('C', 'CA', 'CB'), "angle": 111.65, "note": "Rosetta fa_standard PRO.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'N', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.4906, "mean": 1.4906, "angle_pair": ('N', 'CB', 'CG'), "angle": 71.66, "note": "Rosetta fa_standard PRO.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CA', 'CG', 'CD'): {"bond": ('CG', 'CD'), "value": 1.5055, "mean": 1.5055, "angle_pair": ('CA', 'CG', 'CD'), "angle": 72.18, "note": "Rosetta fa_standard PRO.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "SER": {
        ('N', 'C', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5163, "mean": 1.5163, "angle_pair": ('C', 'CA', 'CB'), "angle": 109.85, "note": "Rosetta fa_standard SER.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'N', 'CB', 'OG'): {"bond": ('CB', 'OG'), "value": 1.4012, "mean": 1.4012, "angle_pair": ('N', 'CB', 'OG'), "angle": 76.85, "note": "Rosetta fa_standard SER.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "THR": {
        ('N', 'C', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5399, "mean": 1.5399, "angle_pair": ('C', 'CA', 'CB'), "angle": 109.64, "note": "Rosetta fa_standard THR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'N', 'CB', 'OG1'): {"bond": ('CB', 'OG1'), "value": 1.4335, "mean": 1.4335, "angle_pair": ('N', 'CB', 'OG1'), "angle": 76.20, "note": "Rosetta fa_standard THR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'OG1', 'CB', 'CG2'): {"bond": ('CB', 'CG2'), "value": 1.5210, "mean": 1.5210, "angle_pair": ('OG1', 'CB', 'CG2'), "angle": 109.32, "note": "Rosetta fa_standard THR.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "TRP": {
        ('N', 'C', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5298, "mean": 1.5298, "angle_pair": ('C', 'CA', 'CB'), "angle": 109.63, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'N', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.4987, "mean": 1.4987, "angle_pair": ('N', 'CB', 'CG'), "angle": 79.69, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CA', 'CG', 'CD1'): {"bond": ('CG', 'CD1'), "value": 1.3627, "mean": 1.3627, "angle_pair": ('CA', 'CG', 'CD1'), "angle": 93.12, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CB', 'CD1', 'NE1'): {"bond": ('CD1', 'NE1'), "value": 1.3729, "mean": 1.3729, "angle_pair": ('CB', 'CD1', 'NE1'), "angle": 138.15, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD1', 'CG', 'NE1', 'CE2'): {"bond": ('NE1', 'CE2'), "value": 1.3721, "mean": 1.3721, "angle_pair": ('CG', 'NE1', 'CE2'), "angle": 74.13, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('NE1', 'CD1', 'CE2', 'CZ2'): {"bond": ('CE2', 'CZ2'), "value": 1.3859, "mean": 1.3859, "angle_pair": ('CD1', 'CE2', 'CZ2'), "angle": 165.70, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CE2', 'NE1', 'CZ2', 'CH2'): {"bond": ('CZ2', 'CH2'), "value": 1.3950, "mean": 1.3950, "angle_pair": ('NE1', 'CZ2', 'CH2'), "angle": 142.30, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CZ2', 'CE2', 'CH2', 'CZ3'): {"bond": ('CH2', 'CZ3'), "value": 1.3721, "mean": 1.3721, "angle_pair": ('CE2', 'CH2', 'CZ3'), "angle": 90.36, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CH2', 'CZ2', 'CZ3', 'CE3'): {"bond": ('CZ3', 'CE3'), "value": 1.3898, "mean": 1.3898, "angle_pair": ('CZ2', 'CZ3', 'CE3'), "angle": 91.48, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CZ3', 'CH2', 'CE3', 'CD2'): {"bond": ('CE3', 'CD2'), "value": 1.4004, "mean": 1.4004, "angle_pair": ('CH2', 'CE3', 'CD2'), "angle": 89.44, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "TYR": {
        ('N', 'C', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5304, "mean": 1.5304, "angle_pair": ('C', 'CA', 'CB'), "angle": 110.30, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'N', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.5127, "mean": 1.5127, "angle_pair": ('N', 'CB', 'CG'), "angle": 80.00, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CA', 'CG', 'CD1'): {"bond": ('CG', 'CD1'), "value": 1.3872, "mean": 1.3872, "angle_pair": ('CA', 'CG', 'CD1'), "angle": 86.68, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CB', 'CD1', 'CE1'): {"bond": ('CD1', 'CE1'), "value": 1.3816, "mean": 1.3816, "angle_pair": ('CB', 'CD1', 'CE1'), "angle": 151.43, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CD1', 'CG', 'CD2'): {"bond": ('CG', 'CD2'), "value": 1.3869, "mean": 1.3869, "angle_pair": ('CD1', 'CG', 'CD2'), "angle": 120.00, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CB', 'CD2', 'CE2'): {"bond": ('CD2', 'CE2'), "value": 1.3814, "mean": 1.3814, "angle_pair": ('CB', 'CD2', 'CE2'), "angle": 151.43, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD2', 'CG', 'CE2', 'CZ'): {"bond": ('CE2', 'CZ'), "value": 1.3800, "mean": 1.3800, "angle_pair": ('CG', 'CE2', 'CZ'), "angle": 89.93, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CE2', 'CD2', 'CZ', 'OH'): {"bond": ('CZ', 'OH'), "value": 1.3760, "mean": 1.3760, "angle_pair": ('CD2', 'CZ', 'OH'), "angle": 150.02, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "VAL": {
        ('N', 'C', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5402, "mean": 1.5402, "angle_pair": ('C', 'CA', 'CB'), "angle": 109.37, "note": "Rosetta fa_standard VAL.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'N', 'CB', 'CG1'): {"bond": ('CB', 'CG1'), "value": 1.5214, "mean": 1.5214, "angle_pair": ('N', 'CB', 'CG1'), "angle": 76.73, "note": "Rosetta fa_standard VAL.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CG1', 'CB', 'CG2'): {"bond": ('CB', 'CG2'), "value": 1.5211, "mean": 1.5211, "angle_pair": ('CG1', 'CB', 'CG2'), "angle": 110.77, "note": "Rosetta fa_standard VAL.params ICOOR ideal", "source": "rosetta_params_408"},
    },
}

BOND_SIDECHAIN_REFS = {
    "rosetta_params_408": "Rosetta 408 main/database/chemical/residue_type_sets/fa_standard/residue_types/l-caa/*.params (ICOOR_INTERNAL 理想键长/键角)。",
}
