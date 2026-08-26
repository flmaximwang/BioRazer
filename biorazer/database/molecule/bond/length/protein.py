# -*- coding: utf-8 -*-
"""Protein-specific covalent bond lengths (Å): per-residue refinements.

数据集来源与诚实性说明
────────────────────────
本模块记录蛋白质**逐残基细分**的共价键长理想值:

* ``AMINO_ACID_BOND_LENGTH_BY_RESIDUE`` -- Engh & Huber (1991) 对同一化学
  键按残基区分的目标值 (key 为单残基型 Gly/Pro/Ala 或同类型残基组
  VIT=Val,Ile,Thr), 覆盖 :mod:`..length.generic` 主表的通用值。
* ``AMINO_ACID_SIDECHAIN_BOND`` -- 20 种标准氨基酸**侧链** (CA 的 CB 及
  以外重原子, 见 :data:`biorazer.database.molecule.icoor.protein.topology.IC_PATH`
  的生长路径) 的理想键长。数值为**理想点值 (ideal point values)**,
  翻译自 Rosetta 408 的 ``fa_standard`` 残基 .params (``ICOOR_INTERNAL``)
  构建的规范残基几何 —— Rosetta 的主链/侧链理想几何源自 Engh & Huber
  (1991) (见 :mod:`..length.generic` 的交叉核对)。键长即 ICOOR 的 ``d``。

注意: 侧链表**没有** Engh & Huber 的样本标准差 —— Rosetta ICOOR 只给
理想点值, 不给 CSD sigma; 需要 std 时请以 E&H 1991 原始文献为准。
故侧链条目 std/lb/up 为 ``np.nan``。

键长单位 **埃 (Å)**。侧链键/角的 key 采用 :data:`IC_PATH` 的生长四元组
(i,j,k,l): 从父原子 (i,j,k) 生长出 l, 键 ``(k,l)``, 键角 ``(j,k,l)``
(顶点在 k, 见 :mod:`..angle.protein`)。
"""

import numpy as np

from .generic import BOND_REFS

#: 键长按残基细分 (Engh & Huber 1991)。key: 单残基型 (Gly/Pro/Ala) 或
#: 同类型残基组 (VIT=Val,Ile,Thr)。一个键写多条同值时以残基覆盖通用主表。
AMINO_ACID_BOND_LENGTH_BY_RESIDUE = {
    "Gly": {
        ("CA", "C"): {"mean": 1.516, "std": 0.018, "lb": 1.462, "up": 1.570,
                       "note": "Gly CH2G+C 1.516±0.018。", "source": "engh_huber_1991"},
        ("N", "CA"): {"mean": 1.451, "std": 0.016, "lb": 1.403, "up": 1.499,
                       "note": "Gly NH1-CH2G+ 1.451±0.016。", "source": "engh_huber_1991"},
    },
    "Pro": {
        ("C", "N"): {"mean": 1.341, "std": 0.016, "lb": 1.293, "up": 1.389,
                       "note": "Pro C-N 1.341±0.016。", "source": "engh_huber_1991"},
        ("N", "CA"): {"mean": 1.466, "std": 0.015, "lb": 1.421, "up": 1.511,
                       "note": "Pro N-CH1E 1.466±0.015。", "source": "engh_huber_1991"},
    },
    "Ala": {
        ("CA", "CB"): {"mean": 1.521, "std": 0.033, "lb": 1.422, "up": 1.620,
                        "note": "Ala CH1E-CH3E 1.521±0.033。", "source": "engh_huber_1991"},
    },
    "VIT": {  # Val, Ile, Thr: Cbeta 为 CH1E (分支)
        ("CA", "CB"): {"mean": 1.540, "std": 0.027, "lb": 1.459, "up": 1.621,
                        "note": "Val/Ile/Thr CH1E-CH1E 1.540±0.027。", "source": "engh_huber_1991"},
    },
}

#: 每种氨基酸侧链生长的规范键长 (k,l)
#: {(res_name): {(i,j,k,l)取原子名的四元组: {...}}}
#: 含第一个 CB 生长四元组 (C,N,CA,CB) (其 CA-CB 键/角度亦见 generic 主表)。
#: std/lb/up 为 np.nan (Rosetta ICOOR 只给理想点值, 无 CSD sigma)。
AMINO_ACID_SIDECHAIN_BOND = {
    "ALA": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5217, "mean": 1.5217, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.37, "note": "Rosetta fa_standard ALA.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "ARG": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5216, "mean": 1.5216, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.60, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.5204, "mean": 1.5204, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CA', 'CB', 'CG'), "angle": 114.87, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'CD'): {"bond": ('CG', 'CD'), "value": 1.4854, "mean": 1.4854, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CB', 'CG', 'CD'), "angle": 111.70, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD', 'NE'): {"bond": ('CD', 'NE'), "value": 1.4541, "mean": 1.4541, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CG', 'CD', 'NE'), "angle": 111.90, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CD', 'NE', 'CZ'): {"bond": ('NE', 'CZ'), "value": 1.3473, "mean": 1.3473, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CD', 'NE', 'CZ'), "angle": 124.60, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD', 'NE', 'CZ', 'NH1'): {"bond": ('CZ', 'NH1'), "value": 1.3146, "mean": 1.3146, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('NE', 'CZ', 'NH1'), "angle": 120.00, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
        ('NH1', 'NE', 'CZ', 'NH2'): {"bond": ('CZ', 'NH2'), "value": 1.3216, "mean": 1.3216, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('NE', 'CZ', 'NH2'), "angle": 120.00, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "ASN": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5177, "mean": 1.5177, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.60, "note": "Rosetta fa_standard ASN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.5035, "mean": 1.5035, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CA', 'CB', 'CG'), "angle": 112.60, "note": "Rosetta fa_standard ASN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'OD1'): {"bond": ('CG', 'OD1'), "value": 1.2364, "mean": 1.2364, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CB', 'CG', 'OD1'), "angle": 120.80, "note": "Rosetta fa_standard ASN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('OD1', 'CB', 'CG', 'ND2'): {"bond": ('CG', 'ND2'), "value": 1.3086, "mean": 1.3086, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CB', 'CG', 'ND2'), "angle": 116.50, "note": "Rosetta fa_standard ASN.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "ASP": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5307, "mean": 1.5307, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.52, "note": "Rosetta fa_standard ASP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.5228, "mean": 1.5228, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CA', 'CB', 'CG'), "angle": 112.90, "note": "Rosetta fa_standard ASP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'OD1'): {"bond": ('CG', 'OD1'), "value": 1.2082, "mean": 1.2082, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CB', 'CG', 'OD1'), "angle": 118.39, "note": "Rosetta fa_standard ASP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('OD1', 'CB', 'CG', 'OD2'): {"bond": ('CG', 'OD2'), "value": 1.2078, "mean": 1.2078, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CB', 'CG', 'OD2'), "angle": 118.37, "note": "Rosetta fa_standard ASP.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "CYS": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5289, "mean": 1.5289, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.60, "note": "Rosetta fa_standard CYS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'SG'): {"bond": ('CB', 'SG'), "value": 1.8088, "mean": 1.8088, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CA', 'CB', 'SG'), "angle": 114.10, "note": "Rosetta fa_standard CYS.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "GLN": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5311, "mean": 1.5311, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.49, "note": "Rosetta fa_standard GLN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.5191, "mean": 1.5191, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CA', 'CB', 'CG'), "angle": 114.19, "note": "Rosetta fa_standard GLN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'CD'): {"bond": ('CG', 'CD'), "value": 1.5169, "mean": 1.5169, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CB', 'CG', 'CD'), "angle": 112.43, "note": "Rosetta fa_standard GLN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD', 'OE1'): {"bond": ('CD', 'OE1'), "value": 1.2342, "mean": 1.2342, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CG', 'CD', 'OE1'), "angle": 120.95, "note": "Rosetta fa_standard GLN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('OE1', 'CG', 'CD', 'NE2'): {"bond": ('CD', 'NE2'), "value": 1.3281, "mean": 1.3281, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CG', 'CD', 'NE2'), "angle": 116.41, "note": "Rosetta fa_standard GLN.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "GLU": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5303, "mean": 1.5303, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.40, "note": "Rosetta fa_standard GLU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.5221, "mean": 1.5221, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CA', 'CB', 'CG'), "angle": 114.40, "note": "Rosetta fa_standard GLU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'CD'): {"bond": ('CG', 'CD'), "value": 1.5034, "mean": 1.5034, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CB', 'CG', 'CD'), "angle": 112.90, "note": "Rosetta fa_standard GLU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD', 'OE1'): {"bond": ('CD', 'OE1'), "value": 1.2076, "mean": 1.2076, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CG', 'CD', 'OE1'), "angle": 118.45, "note": "Rosetta fa_standard GLU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('OE1', 'CG', 'CD', 'OE2'): {"bond": ('CD', 'OE2'), "value": 1.2085, "mean": 1.2085, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CG', 'CD', 'OE2'), "angle": 118.36, "note": "Rosetta fa_standard GLU.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "GLY": {
    },
    "HIS": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5321, "mean": 1.5321, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.69, "note": "Rosetta fa_standard HIS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.4972, "mean": 1.4972, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CA', 'CB', 'CG'), "angle": 113.69, "note": "Rosetta fa_standard HIS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'ND1'): {"bond": ('CG', 'ND1'), "value": 1.3792, "mean": 1.3792, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CB', 'CG', 'ND1'), "angle": 122.63, "note": "Rosetta fa_standard HIS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'ND1', 'CE1'): {"bond": ('ND1', 'CE1'), "value": 1.3219, "mean": 1.3219, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CG', 'ND1', 'CE1'), "angle": 109.29, "note": "Rosetta fa_standard HIS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'ND1', 'CE1', 'NE2'): {"bond": ('CE1', 'NE2'), "value": 1.3204, "mean": 1.3204, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('ND1', 'CE1', 'NE2'), "angle": 108.38, "note": "Rosetta fa_standard HIS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('ND1', 'CE1', 'NE2', 'CD2'): {"bond": ('NE2', 'CD2'), "value": 1.3732, "mean": 1.3732, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CE1', 'NE2', 'CD2'), "angle": 109.01, "note": "Rosetta fa_standard HIS.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "ILE": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5396, "mean": 1.5396, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.90, "note": "Rosetta fa_standard ILE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG1'): {"bond": ('CB', 'CG1'), "value": 1.5309, "mean": 1.5309, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CA', 'CB', 'CG1'), "angle": 110.41, "note": "Rosetta fa_standard ILE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG1', 'CD1'): {"bond": ('CG1', 'CD1'), "value": 1.5117, "mean": 1.5117, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CB', 'CG1', 'CD1'), "angle": 113.83, "note": "Rosetta fa_standard ILE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG1', 'CA', 'CB', 'CG2'): {"bond": ('CB', 'CG2'), "value": 1.5209, "mean": 1.5209, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CA', 'CB', 'CG2'), "angle": 110.47, "note": "Rosetta fa_standard ILE.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "LEU": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5339, "mean": 1.5339, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.21, "note": "Rosetta fa_standard LEU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.5340, "mean": 1.5340, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CA', 'CB', 'CG'), "angle": 115.70, "note": "Rosetta fa_standard LEU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'CD1'): {"bond": ('CG', 'CD1'), "value": 1.5227, "mean": 1.5227, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CB', 'CG', 'CD1'), "angle": 109.50, "note": "Rosetta fa_standard LEU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD1', 'CB', 'CG', 'CD2'): {"bond": ('CG', 'CD2'), "value": 1.5214, "mean": 1.5214, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CB', 'CG', 'CD2'), "angle": 109.50, "note": "Rosetta fa_standard LEU.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "LYS": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5295, "mean": 1.5295, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.50, "note": "Rosetta fa_standard LYS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.5229, "mean": 1.5229, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CA', 'CB', 'CG'), "angle": 114.40, "note": "Rosetta fa_standard LYS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'CD'): {"bond": ('CG', 'CD'), "value": 1.5213, "mean": 1.5213, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CB', 'CG', 'CD'), "angle": 111.30, "note": "Rosetta fa_standard LYS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD', 'CE'): {"bond": ('CD', 'CE'), "value": 1.5216, "mean": 1.5216, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CG', 'CD', 'CE'), "angle": 111.37, "note": "Rosetta fa_standard LYS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CD', 'CE', 'NZ'): {"bond": ('CE', 'NZ'), "value": 1.4881, "mean": 1.4881, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CD', 'CE', 'NZ'), "angle": 111.96, "note": "Rosetta fa_standard LYS.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "MET": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5274, "mean": 1.5274, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.21, "note": "Rosetta fa_standard MET.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.5222, "mean": 1.5222, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CA', 'CB', 'CG'), "angle": 114.44, "note": "Rosetta fa_standard MET.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'SD'): {"bond": ('CG', 'SD'), "value": 1.8038, "mean": 1.8038, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CB', 'CG', 'SD'), "angle": 112.67, "note": "Rosetta fa_standard MET.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'SD', 'CE'): {"bond": ('SD', 'CE'), "value": 1.7904, "mean": 1.7904, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CG', 'SD', 'CE'), "angle": 100.89, "note": "Rosetta fa_standard MET.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "PHE": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5298, "mean": 1.5298, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.44, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.5022, "mean": 1.5022, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CA', 'CB', 'CG'), "angle": 113.79, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'CD1'): {"bond": ('CG', 'CD1'), "value": 1.3870, "mean": 1.3870, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CB', 'CG', 'CD1'), "angle": 120.69, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD1', 'CE1'): {"bond": ('CD1', 'CE1'), "value": 1.3822, "mean": 1.3822, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CG', 'CD1', 'CE1'), "angle": 120.73, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CD1', 'CE1', 'CZ'): {"bond": ('CE1', 'CZ'), "value": 1.3786, "mean": 1.3786, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CD1', 'CE1', 'CZ'), "angle": 120.03, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD1', 'CE1', 'CZ', 'CE2'): {"bond": ('CZ', 'CE2'), "value": 1.3805, "mean": 1.3805, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CE1', 'CZ', 'CE2'), "angle": 119.87, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CE1', 'CZ', 'CE2', 'CD2'): {"bond": ('CE2', 'CD2'), "value": 1.3813, "mean": 1.3813, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CZ', 'CE2', 'CD2'), "angle": 119.97, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "PRO": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5320, "mean": 1.5320, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('N', 'CA', 'CB'), "angle": 103.00, "note": "Rosetta fa_standard PRO.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.4906, "mean": 1.4906, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CA', 'CB', 'CG'), "angle": 104.20, "note": "Rosetta fa_standard PRO.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'CD'): {"bond": ('CG', 'CD'), "value": 1.5055, "mean": 1.5055, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CB', 'CG', 'CD'), "angle": 104.40, "note": "Rosetta fa_standard PRO.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "SER": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5163, "mean": 1.5163, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.14, "note": "Rosetta fa_standard SER.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'OG'): {"bond": ('CB', 'OG'), "value": 1.4012, "mean": 1.4012, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CA', 'CB', 'OG'), "angle": 111.00, "note": "Rosetta fa_standard SER.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "THR": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5399, "mean": 1.5399, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('N', 'CA', 'CB'), "angle": 111.10, "note": "Rosetta fa_standard THR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'OG1'): {"bond": ('CB', 'OG1'), "value": 1.4335, "mean": 1.4335, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CA', 'CB', 'OG1'), "angle": 109.58, "note": "Rosetta fa_standard THR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('OG1', 'CA', 'CB', 'CG2'): {"bond": ('CB', 'CG2'), "value": 1.5210, "mean": 1.5210, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CA', 'CB', 'CG2'), "angle": 110.53, "note": "Rosetta fa_standard THR.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "TRP": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5298, "mean": 1.5298, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.40, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.4987, "mean": 1.4987, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CA', 'CB', 'CG'), "angle": 113.53, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'CD1'): {"bond": ('CG', 'CD1'), "value": 1.3627, "mean": 1.3627, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CB', 'CG', 'CD1'), "angle": 126.74, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD1', 'NE1'): {"bond": ('CD1', 'NE1'), "value": 1.3729, "mean": 1.3729, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CG', 'CD1', 'NE1'), "angle": 110.16, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CD1', 'NE1', 'CE2'): {"bond": ('NE1', 'CE2'), "value": 1.3721, "mean": 1.3721, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CD1', 'NE1', 'CE2'), "angle": 108.90, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD1', 'NE1', 'CE2', 'CZ2'): {"bond": ('CE2', 'CZ2'), "value": 1.3859, "mean": 1.3859, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('NE1', 'CE2', 'CZ2'), "angle": 130.14, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('NE1', 'CE2', 'CZ2', 'CH2'): {"bond": ('CZ2', 'CH2'), "value": 1.3950, "mean": 1.3950, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CE2', 'CZ2', 'CH2'), "angle": 117.50, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CE2', 'CZ2', 'CH2', 'CZ3'): {"bond": ('CH2', 'CZ3'), "value": 1.3721, "mean": 1.3721, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CZ2', 'CH2', 'CZ3'), "angle": 121.50, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CZ2', 'CH2', 'CZ3', 'CE3'): {"bond": ('CZ3', 'CE3'), "value": 1.3898, "mean": 1.3898, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CH2', 'CZ3', 'CE3'), "angle": 121.00, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CH2', 'CZ3', 'CE3', 'CD2'): {"bond": ('CE3', 'CD2'), "value": 1.4004, "mean": 1.4004, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CZ3', 'CE3', 'CD2'), "angle": 118.73, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "TYR": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5304, "mean": 1.5304, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.47, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG'): {"bond": ('CB', 'CG'), "value": 1.5127, "mean": 1.5127, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CA', 'CB', 'CG'), "angle": 113.80, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG', 'CD1'): {"bond": ('CG', 'CD1'), "value": 1.3872, "mean": 1.3872, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CB', 'CG', 'CD1'), "angle": 120.00, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD1', 'CE1'): {"bond": ('CD1', 'CE1'), "value": 1.3816, "mean": 1.3816, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CG', 'CD1', 'CE1'), "angle": 120.00, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD1', 'CB', 'CG', 'CD2'): {"bond": ('CG', 'CD2'), "value": 1.3869, "mean": 1.3869, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CB', 'CG', 'CD2'), "angle": 120.00, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD2', 'CE2'): {"bond": ('CD2', 'CE2'), "value": 1.3814, "mean": 1.3814, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CG', 'CD2', 'CE2'), "angle": 120.00, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CD2', 'CE2', 'CZ'): {"bond": ('CE2', 'CZ'), "value": 1.3800, "mean": 1.3800, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CD2', 'CE2', 'CZ'), "angle": 120.00, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD2', 'CE2', 'CZ', 'OH'): {"bond": ('CZ', 'OH'), "value": 1.3760, "mean": 1.3760, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CE2', 'CZ', 'OH'), "angle": 120.00, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "VAL": {
        ('C', 'N', 'CA', 'CB'): {"bond": ('CA', 'CB'), "value": 1.5402, "mean": 1.5402, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('N', 'CA', 'CB'), "angle": 110.30, "note": "Rosetta fa_standard VAL.params ICOOR ideal", "source": "rosetta_params_408"},
        ('N', 'CA', 'CB', 'CG1'): {"bond": ('CB', 'CG1'), "value": 1.5214, "mean": 1.5214, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CA', 'CB', 'CG1'), "angle": 110.49, "note": "Rosetta fa_standard VAL.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG1', 'CA', 'CB', 'CG2'): {"bond": ('CB', 'CG2'), "value": 1.5211, "mean": 1.5211, "std": np.nan, "lb": np.nan, "up": np.nan, "angle_pair": ('CA', 'CB', 'CG2'), "angle": 109.90, "note": "Rosetta fa_standard VAL.params ICOOR ideal", "source": "rosetta_params_408"},
    },
}

BOND_SIDECHAIN_REFS = {
    "rosetta_params_408": "Rosetta 408 main/database/chemical/residue_type_sets/fa_standard/residue_types/l-caa/*.params (ICOOR_INTERNAL 理想键长/键角)。",
}
