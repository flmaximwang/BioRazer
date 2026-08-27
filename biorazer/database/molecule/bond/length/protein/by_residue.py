# -*- coding: utf-8 -*-
"""Protein bond-length reference data: generic, per-residue, and side-chain.

数据集来源与诚实性说明
────────────────────────
本模块 (``length/protein/by_residue``) 汇聚蛋白质共价键长的全部参考数据,
按「通用 → 逐残基 → 侧链」组织在同一模块:

* ``AMINO_ACID_BOND_LENGTH`` -- 通用主表 (Engh & Huber 1991 对同一化学键
  最「一般」的目标值, key 为两端原子名二元组, 如 ``("C", "N")``)。
  每个 key 至少被 2 个残基引用 (键长只按残基细分到 Gly/Pro/Ala 及 VIT=
  Val,Ile,Thr 组, 其余残基直接取该通用值)。
* ``AMINO_ACID_BOND_LENGTH_BY_RESIDUE`` -- 同一化学键按残基细分的目标值
  (key 为单残基型 Gly/Pro/Ala 或同类型残基组 VIT), 覆盖通用主表。
* ``AMINO_ACID_SIDECHAIN_BOND`` -- 20 种标准氨基酸**侧链** (CA 的 CB 及
  以外重原子) 的理想键长。数值为**理想点值 (ideal point values)**, 翻译自
  Rosetta 408 的 ``fa_standard`` 残基 .params (``ICOOR_INTERNAL``) 构建的
  规范残基几何。侧链表**没有** Engh & Huber 的样本标准差, 故 std/lb/up 为
  ``np.nan``。

数值与标准差**全部来自 Engh & Huber (1991) X 射线蛋白质结构精修参数**
(统计自 CSD 小分子晶体结构), 主表 (含 sigma) 取自 PROCHECK 手册
Appendix A (Laskowski et al., 1993), 并与 Rosetta 408 氨基酸 .params
交叉核对一致 (见 :data:`BOND_REFS` 各条)。

键长单位 **埃 (Å)**。侧链键的 key 为 **2 原子元组 ``(k, l)``** (被生长
原子 ``l`` 与其成键父原子 ``k``); 对应键角见
:data:`~biorazer.database.molecule.bond.angle.protein.AMINO_ACID_BOND_ANGLE`
(key 为 3 原子元组), 二面角见
:data:`~biorazer.database.molecule.bond.dihedral.protein.SIDECHAIN_IC_DIHEDRAL`
(key 为生长四元组)。
"""

import numpy as np

#: 文献库 (provenance)。每条记录的 ``source`` 指向这里的键。
BOND_REFS = {
    "engh_huber_1991": "Engh RA, Huber R. Accurate bond and angle parameters for X-ray protein structure refinement. Acta Cryst. A47:392-400, 1991. (键长/键角主参数, 统计自 CSD 小分子晶体结构)",
    "procheck_appendix_a": "Laskowski RA, MacArthur MW, Moss DS, Thornton JM. PROCHECK: a program to check the stereochemical quality of protein structures. J. Appl. Cryst. 26:283-291, 1993. (手册 Appendix A.2 原封转载 Engh & Huber 1991 的 X-PLOR 标号键长/键角表)",
    "rosetta_params_408": "Rosetta 408 main/database/chemical/residue_type_sets/fa_standard/residue_types/l-caa/*.params (ICOOR_INTERNAL 主链内坐标), 与本模块主表交叉核对一致。",
}

# ---------------------------------------------------------------------------
# 键长主表 (通用, 单位: Å)
# ---------------------------------------------------------------------------
# 每个键以两端原子名二元组为 key。mean/std 为 Engh & Huber (1991) 目标值与
# CSD 样本标准差; lb/up = mean ± 3*std。每个 key 至少被 2 个残基引用。
AMINO_ACID_BOND_LENGTH = {
    # 肽键 C-N (非 Pro)。Pro 的 C-N 为 1.341±0.016, 见 protein._BY_RESIDUE。
    ("C", "N"): {
        "mean": 1.329, "std": 0.014, "lb": 1.287, "up": 1.371,
        "note": "肽键 C(=O)-N。Engh-Huber 主表 C-NH1 (except Pro) 1.329±0.014; Pro 1.341±0.016。",
        "source": ("engh_huber_1991", "procheck_appendix_a"),
    },
    # 羰基 C=O
    ("C", "O"): {
        "mean": 1.231, "std": 0.020, "lb": 1.171, "up": 1.291,
        "note": "肽羰基 C=O。C-O 1.231±0.020 (所有残基同一值)。",
        "source": ("engh_huber_1991", "procheck_appendix_a"),
    },
    # CA-C (非 Gly)
    ("CA", "C"): {
        "mean": 1.525, "std": 0.021, "lb": 1.462, "up": 1.588,
        "note": "Calpha-C。Engh-Huber CH1E-C (except Gly) 1.525±0.021; Gly 1.516±0.018。",
        "source": ("engh_huber_1991", "procheck_appendix_a"),
    },
    # N-CA (非 Gly/Pro)
    ("N", "CA"): {
        "mean": 1.458, "std": 0.019, "lb": 1.401, "up": 1.515,
        "note": "N-Calpha。Engh-Huber NH1-CH1E (except Gly,Pro) 1.458±0.019; Gly 1.451±0.016; Pro 1.466±0.015 (见 protein._BY_RESIDUE)。",
        "source": ("engh_huber_1991", "procheck_appendix_a"),
    },
    # CA-CB (其余残基, 非 Ala / 非 Val,Ile,Thr)
    ("CA", "CB"): {
        "mean": 1.530, "std": 0.020, "lb": 1.470, "up": 1.590,
        "note": "Calpha-Cbeta。Engh-Huber CH1E-CH2E (the rest) 1.530±0.020; Ala 1.521±0.033; Val/Ile/Thr 1.540±0.027 (见 protein._BY_RESIDUE)。",
        "source": ("engh_huber_1991", "procheck_appendix_a"),
    },
}

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

AMINO_ACID_SIDECHAIN_BOND = {
    "ALA": {
        ('CA', 'CB'): {"mean": 1.5217, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ALA.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "ARG": {
        ('CA', 'CB'): {"mean": 1.5216, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG'): {"mean": 1.5204, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CD'): {"mean": 1.4854, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD', 'NE'): {"mean": 1.4541, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
        ('NE', 'CZ'): {"mean": 1.3473, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CZ', 'NH1'): {"mean": 1.3146, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CZ', 'NH2'): {"mean": 1.3216, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "ASN": {
        ('CA', 'CB'): {"mean": 1.5177, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ASN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG'): {"mean": 1.5035, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ASN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'OD1'): {"mean": 1.2364, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ASN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'ND2'): {"mean": 1.3086, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ASN.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "ASP": {
        ('CA', 'CB'): {"mean": 1.5307, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ASP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG'): {"mean": 1.5228, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ASP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'OD1'): {"mean": 1.2082, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ASP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'OD2'): {"mean": 1.2078, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ASP.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "CYS": {
        ('CA', 'CB'): {"mean": 1.5289, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard CYS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'SG'): {"mean": 1.8088, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard CYS.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "GLN": {
        ('CA', 'CB'): {"mean": 1.5311, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard GLN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG'): {"mean": 1.5191, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard GLN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CD'): {"mean": 1.5169, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard GLN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD', 'OE1'): {"mean": 1.2342, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard GLN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD', 'NE2'): {"mean": 1.3281, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard GLN.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "GLU": {
        ('CA', 'CB'): {"mean": 1.5303, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard GLU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG'): {"mean": 1.5221, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard GLU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CD'): {"mean": 1.5034, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard GLU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD', 'OE1'): {"mean": 1.2076, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard GLU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD', 'OE2'): {"mean": 1.2085, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard GLU.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "GLY": {
    },
    "HIS": {
        ('CA', 'CB'): {"mean": 1.5321, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard HIS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG'): {"mean": 1.4972, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard HIS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'ND1'): {"mean": 1.3792, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard HIS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('ND1', 'CE1'): {"mean": 1.3219, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard HIS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CE1', 'NE2'): {"mean": 1.3204, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard HIS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('NE2', 'CD2'): {"mean": 1.3732, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard HIS.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "ILE": {
        ('CA', 'CB'): {"mean": 1.5396, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ILE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG1'): {"mean": 1.5309, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ILE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG1', 'CD1'): {"mean": 1.5117, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ILE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG2'): {"mean": 1.5209, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ILE.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "LEU": {
        ('CA', 'CB'): {"mean": 1.5339, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard LEU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG'): {"mean": 1.5340, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard LEU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CD1'): {"mean": 1.5227, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard LEU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CD2'): {"mean": 1.5214, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard LEU.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "LYS": {
        ('CA', 'CB'): {"mean": 1.5295, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard LYS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG'): {"mean": 1.5229, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard LYS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CD'): {"mean": 1.5213, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard LYS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD', 'CE'): {"mean": 1.5216, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard LYS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CE', 'NZ'): {"mean": 1.4881, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard LYS.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "MET": {
        ('CA', 'CB'): {"mean": 1.5274, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard MET.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG'): {"mean": 1.5222, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard MET.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'SD'): {"mean": 1.8038, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard MET.params ICOOR ideal", "source": "rosetta_params_408"},
        ('SD', 'CE'): {"mean": 1.7904, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard MET.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "PHE": {
        ('CA', 'CB'): {"mean": 1.5298, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG'): {"mean": 1.5022, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CD1'): {"mean": 1.3870, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD1', 'CE1'): {"mean": 1.3822, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CE1', 'CZ'): {"mean": 1.3786, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CZ', 'CE2'): {"mean": 1.3805, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CE2', 'CD2'): {"mean": 1.3813, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "PRO": {
        ('CA', 'CB'): {"mean": 1.5320, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard PRO.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG'): {"mean": 1.4906, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard PRO.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CD'): {"mean": 1.5055, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard PRO.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "SER": {
        ('CA', 'CB'): {"mean": 1.5163, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard SER.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'OG'): {"mean": 1.4012, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard SER.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "THR": {
        ('CA', 'CB'): {"mean": 1.5399, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard THR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'OG1'): {"mean": 1.4335, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard THR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG2'): {"mean": 1.5210, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard THR.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "TRP": {
        ('CA', 'CB'): {"mean": 1.5298, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG'): {"mean": 1.4987, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CD1'): {"mean": 1.3627, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD1', 'NE1'): {"mean": 1.3729, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('NE1', 'CE2'): {"mean": 1.3721, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CE2', 'CZ2'): {"mean": 1.3859, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CZ2', 'CH2'): {"mean": 1.3950, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CH2', 'CZ3'): {"mean": 1.3721, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CZ3', 'CE3'): {"mean": 1.3898, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CE3', 'CD2'): {"mean": 1.4004, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "TYR": {
        ('CA', 'CB'): {"mean": 1.5304, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG'): {"mean": 1.5127, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CD1'): {"mean": 1.3872, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD1', 'CE1'): {"mean": 1.3816, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CD2'): {"mean": 1.3869, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD2', 'CE2'): {"mean": 1.3814, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CE2', 'CZ'): {"mean": 1.3800, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CZ', 'OH'): {"mean": 1.3760, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "VAL": {
        ('CA', 'CB'): {"mean": 1.5402, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard VAL.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG1'): {"mean": 1.5214, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard VAL.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG2'): {"mean": 1.5211, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard VAL.params ICOOR ideal", "source": "rosetta_params_408"},
    },
}

BOND_SIDECHAIN_REFS = {
    "rosetta_params_408": "Rosetta 408 main/database/chemical/residue_type_sets/fa_standard/residue_types/l-caa/*.params (ICOOR_INTERNAL 理想键长/键角)。",
}
