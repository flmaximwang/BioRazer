# -*- coding: utf-8 -*-
"""Protein-specific covalent bond angles (°): per-residue, complete.

数据集来源与诚实性说明
────────────────────────
本模块记录蛋白质**逐残基**的共价键角理想值,key 为 20 种标准氨基酸残基名,
值为该残基的**完整键角表** ``{(j, k, l): {mean, std, lb, up, source}}``:

* **骨架** 键角对每个残基都写全。通用值直接引用
  :data:`~biorazer.database.molecule.bond.angle.generic.AMINO_ACID_BACKBONE_BOND_ANGLE`
  (``AMINO_ACID_BACKBONE_BOND_ANGLE[key]``); 对同一化学键按残基细分
  (Engh & Huber 1991) 时以字面量覆盖 (Gly/Pro/Ala 及同类型残基组
  VIT=Val,Ile,Thr 的相关角)。
* **侧链** (CA 的 CB 及以外重原子, 见
  :data:`biorazer.database.molecule.icoor.protein.topology.IC_PATH` 的
  生长路径) 理想键角。数值为**理想点值 (ideal point values)**, 翻译自
  Rosetta 408 的 ``fa_standard`` 残基 .params (``ICOOR_INTERNAL``) 构建的
  规范残基几何 —— Rosetta 的主链/侧链理想几何源自 Engh & Huber (1991)
  (见 :mod:`..angle.generic` 的交叉核对)。键角即 ``180° − theta`` (解码后
  又在构建出的残基坐标上实测复核)。

注意: 侧链表**没有** Engh & Huber 的样本标准差 —— Rosetta ICOOR 只给
理想点值, 不给 CSD sigma; 需要 std 时请以 E&H 1991 原始文献为准。
故侧链条目 std/lb/up 为 ``np.nan``。

键角单位 **度 (degree)**。键角的 key 为 **3 原子元组 ``(j, k, l)``** (顶点
在 ``k``); 对应键长见
:data:`~biorazer.database.molecule.bond.length.protein.AMINO_ACID_SIDECHAIN_BOND`
(key 为 2 原子元组), 二面角见
:data:`~biorazer.database.molecule.bond.dihedral.protein.SIDECHAIN_IC_DIHEDRAL`
(key 为生长四元组)。

通用骨架角的引用 (避免逐残基重复硬编码): 每个残基的骨架部分通过
``AMINO_ACID_BACKBONE_BOND_ANGLE[(j, k, l)]`` 取 :mod:`..angle.generic`
的通用值, 仅在 Gly/Pro/Ala/VIT 有细化时以字面量覆盖。
"""

import numpy as np

from .generic import (
    AMINO_ACID_BACKBONE_BOND_ANGLE,  # noqa: F401  (generic backbone refs)
    BOND_REFS,                       # noqa: F401  (shared provenance)
)

#: 逐残基完整键角表 (Engh & Huber 1991 骨架 + Rosetta 408 侧链理想点值)。
#: key 为 20 种标准氨基酸残基名。骨架通用值引用 generic 主表, 细分值字面量
#: 覆盖; 侧链条目 std/lb/up 为 np.nan。
AMINO_ACID_BOND_ANGLE = {
    "ALA": {
        # ---- 骨架 (通用值引用 generic; 本残基无细化) ----
        ("N", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("N", "CA", "C")],
        ("CA", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "N")],
        ("CA", "C", "O"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "O")],
        ("C", "N", "CA"): AMINO_ACID_BACKBONE_BOND_ANGLE[("C", "N", "CA")],
        ("O", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("O", "C", "N")],
        ("CB", "CA", "C"): {"mean": 110.5, "std": 1.5, "lb": 106.0, "up": 115.0,
                            "note": "Ala CH3E-CH1E-C 110.5±1.5。", "source": "engh_huber_1991"},
        # ---- 侧链 (Rosetta 408 ICOOR 理想点值) ----
        ('N', 'CA', 'CB'): {"mean": 110.37, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ALA.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "ARG": {
        # ---- 骨架 (通用值引用 generic; 本残基无细化) ----
        ("N", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("N", "CA", "C")],
        ("CA", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "N")],
        ("CA", "C", "O"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "O")],
        ("C", "N", "CA"): AMINO_ACID_BACKBONE_BOND_ANGLE[("C", "N", "CA")],
        ("O", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("O", "C", "N")],
        ("CB", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CB", "CA", "C")],
        # ---- 侧链 (Rosetta 408 ICOOR 理想点值) ----
        ('N', 'CA', 'CB'): {"mean": 110.60, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG'): {"mean": 114.87, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD'): {"mean": 111.70, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CD', 'NE'): {"mean": 111.90, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD', 'NE', 'CZ'): {"mean": 124.60, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
        ('NE', 'CZ', 'NH1'): {"mean": 120.00, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
        ('NE', 'CZ', 'NH2'): {"mean": 120.00, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ARG.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "ASN": {
        # ---- 骨架 (通用值引用 generic; 本残基无细化) ----
        ("N", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("N", "CA", "C")],
        ("CA", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "N")],
        ("CA", "C", "O"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "O")],
        ("C", "N", "CA"): AMINO_ACID_BACKBONE_BOND_ANGLE[("C", "N", "CA")],
        ("O", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("O", "C", "N")],
        ("CB", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CB", "CA", "C")],
        # ---- 侧链 (Rosetta 408 ICOOR 理想点值) ----
        ('N', 'CA', 'CB'): {"mean": 110.60, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ASN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG'): {"mean": 112.60, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ASN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'OD1'): {"mean": 120.80, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ASN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'ND2'): {"mean": 116.50, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ASN.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "ASP": {
        # ---- 骨架 (通用值引用 generic; 本残基无细化) ----
        ("N", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("N", "CA", "C")],
        ("CA", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "N")],
        ("CA", "C", "O"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "O")],
        ("C", "N", "CA"): AMINO_ACID_BACKBONE_BOND_ANGLE[("C", "N", "CA")],
        ("O", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("O", "C", "N")],
        ("CB", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CB", "CA", "C")],
        # ---- 侧链 (Rosetta 408 ICOOR 理想点值) ----
        ('N', 'CA', 'CB'): {"mean": 110.52, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ASP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG'): {"mean": 112.90, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ASP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'OD1'): {"mean": 118.39, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ASP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'OD2'): {"mean": 118.37, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ASP.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "CYS": {
        # ---- 骨架 (通用值引用 generic; 本残基无细化) ----
        ("N", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("N", "CA", "C")],
        ("CA", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "N")],
        ("CA", "C", "O"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "O")],
        ("C", "N", "CA"): AMINO_ACID_BACKBONE_BOND_ANGLE[("C", "N", "CA")],
        ("O", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("O", "C", "N")],
        ("CB", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CB", "CA", "C")],
        # ---- 侧链 (Rosetta 408 ICOOR 理想点值) ----
        ('N', 'CA', 'CB'): {"mean": 110.60, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard CYS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'SG'): {"mean": 114.10, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard CYS.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "GLN": {
        # ---- 骨架 (通用值引用 generic; 本残基无细化) ----
        ("N", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("N", "CA", "C")],
        ("CA", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "N")],
        ("CA", "C", "O"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "O")],
        ("C", "N", "CA"): AMINO_ACID_BACKBONE_BOND_ANGLE[("C", "N", "CA")],
        ("O", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("O", "C", "N")],
        ("CB", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CB", "CA", "C")],
        # ---- 侧链 (Rosetta 408 ICOOR 理想点值) ----
        ('N', 'CA', 'CB'): {"mean": 110.49, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard GLN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG'): {"mean": 114.19, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard GLN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD'): {"mean": 112.43, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard GLN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CD', 'OE1'): {"mean": 120.95, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard GLN.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CD', 'NE2'): {"mean": 116.41, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard GLN.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "GLU": {
        # ---- 骨架 (通用值引用 generic; 本残基无细化) ----
        ("N", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("N", "CA", "C")],
        ("CA", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "N")],
        ("CA", "C", "O"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "O")],
        ("C", "N", "CA"): AMINO_ACID_BACKBONE_BOND_ANGLE[("C", "N", "CA")],
        ("O", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("O", "C", "N")],
        ("CB", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CB", "CA", "C")],
        # ---- 侧链 (Rosetta 408 ICOOR 理想点值) ----
        ('N', 'CA', 'CB'): {"mean": 110.40, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard GLU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG'): {"mean": 114.40, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard GLU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD'): {"mean": 112.90, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard GLU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CD', 'OE1'): {"mean": 118.45, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard GLU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CD', 'OE2'): {"mean": 118.36, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard GLU.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "GLY": {
        # ---- 骨架 (Gly 细化: 112.5/116.4/120.8/120.6) ----
        ("N", "CA", "C"): {"mean": 112.5, "std": 2.9, "lb": 103.8, "up": 121.2,
                            "note": "Gly NH1-CH2G+C 112.5±2.9。", "source": "engh_huber_1991"},
        ("CA", "C", "N"): {"mean": 116.4, "std": 2.1, "lb": 110.1, "up": 122.7,
                            "note": "Gly CH2G+C-NH1 116.4±2.1。", "source": "engh_huber_1991"},
        ("CA", "C", "O"): {"mean": 120.8, "std": 2.1, "lb": 114.5, "up": 127.1,
                            "note": "Gly CH2G+C-O 120.8±2.1。", "source": "engh_huber_1991"},
        ("C", "N", "CA"): {"mean": 120.6, "std": 1.7, "lb": 115.5, "up": 125.7,
                            "note": "Gly C-NH1-CH2G+ 120.6±1.7。", "source": "engh_huber_1991"},
        ("O", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("O", "C", "N")],
        # ---- 侧链 (Gly 无 CB) ----
    },
    "HIS": {
        # ---- 骨架 (通用值引用 generic; 本残基无细化) ----
        ("N", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("N", "CA", "C")],
        ("CA", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "N")],
        ("CA", "C", "O"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "O")],
        ("C", "N", "CA"): AMINO_ACID_BACKBONE_BOND_ANGLE[("C", "N", "CA")],
        ("O", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("O", "C", "N")],
        ("CB", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CB", "CA", "C")],
        # ---- 侧链 (Rosetta 408 ICOOR 理想点值) ----
        ('N', 'CA', 'CB'): {"mean": 110.69, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard HIS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG'): {"mean": 113.69, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard HIS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'ND1'): {"mean": 122.63, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard HIS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'ND1', 'CE1'): {"mean": 109.29, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard HIS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('ND1', 'CE1', 'NE2'): {"mean": 108.38, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard HIS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CE1', 'NE2', 'CD2'): {"mean": 109.01, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard HIS.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "ILE": {
        # ---- 骨架 (通用值引用 generic; 本残基无细化) ----
        ("N", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("N", "CA", "C")],
        ("CA", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "N")],
        ("CA", "C", "O"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "O")],
        ("C", "N", "CA"): AMINO_ACID_BACKBONE_BOND_ANGLE[("C", "N", "CA")],
        ("O", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("O", "C", "N")],
        ("CB", "CA", "C"): {"mean": 109.1, "std": 2.2, "lb": 102.5, "up": 115.7,
                            "note": "Val/Ile/Thr CH1E-CH1E-C 109.1±2.2。", "source": "engh_huber_1991"},
        # ---- 侧链 (Rosetta 408 ICOOR 理想点值) ----
        ('N', 'CA', 'CB'): {"mean": 110.90, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ILE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG1'): {"mean": 110.41, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ILE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG1', 'CD1'): {"mean": 113.83, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ILE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG2'): {"mean": 110.47, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard ILE.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "LEU": {
        # ---- 骨架 (通用值引用 generic; 本残基无细化) ----
        ("N", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("N", "CA", "C")],
        ("CA", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "N")],
        ("CA", "C", "O"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "O")],
        ("C", "N", "CA"): AMINO_ACID_BACKBONE_BOND_ANGLE[("C", "N", "CA")],
        ("O", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("O", "C", "N")],
        ("CB", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CB", "CA", "C")],
        # ---- 侧链 (Rosetta 408 ICOOR 理想点值) ----
        ('N', 'CA', 'CB'): {"mean": 110.21, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard LEU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG'): {"mean": 115.70, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard LEU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD1'): {"mean": 109.50, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard LEU.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD2'): {"mean": 109.50, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard LEU.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "LYS": {
        # ---- 骨架 (通用值引用 generic; 本残基无细化) ----
        ("N", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("N", "CA", "C")],
        ("CA", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "N")],
        ("CA", "C", "O"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "O")],
        ("C", "N", "CA"): AMINO_ACID_BACKBONE_BOND_ANGLE[("C", "N", "CA")],
        ("O", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("O", "C", "N")],
        ("CB", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CB", "CA", "C")],
        # ---- 侧链 (Rosetta 408 ICOOR 理想点值) ----
        ('N', 'CA', 'CB'): {"mean": 110.50, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard LYS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG'): {"mean": 114.40, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard LYS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD'): {"mean": 111.30, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard LYS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CD', 'CE'): {"mean": 111.37, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard LYS.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD', 'CE', 'NZ'): {"mean": 111.96, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard LYS.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "MET": {
        # ---- 骨架 (通用值引用 generic; 本残基无细化) ----
        ("N", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("N", "CA", "C")],
        ("CA", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "N")],
        ("CA", "C", "O"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "O")],
        ("C", "N", "CA"): AMINO_ACID_BACKBONE_BOND_ANGLE[("C", "N", "CA")],
        ("O", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("O", "C", "N")],
        ("CB", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CB", "CA", "C")],
        # ---- 侧链 (Rosetta 408 ICOOR 理想点值) ----
        ('N', 'CA', 'CB'): {"mean": 110.21, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard MET.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG'): {"mean": 114.44, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard MET.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'SD'): {"mean": 112.67, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard MET.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'SD', 'CE'): {"mean": 100.89, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard MET.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "PHE": {
        # ---- 骨架 (通用值引用 generic; 本残基无细化) ----
        ("N", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("N", "CA", "C")],
        ("CA", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "N")],
        ("CA", "C", "O"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "O")],
        ("C", "N", "CA"): AMINO_ACID_BACKBONE_BOND_ANGLE[("C", "N", "CA")],
        ("O", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("O", "C", "N")],
        ("CB", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CB", "CA", "C")],
        # ---- 侧链 (Rosetta 408 ICOOR 理想点值) ----
        ('N', 'CA', 'CB'): {"mean": 110.44, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG'): {"mean": 113.79, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD1'): {"mean": 120.69, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CD1', 'CE1'): {"mean": 120.73, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD1', 'CE1', 'CZ'): {"mean": 120.03, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CE1', 'CZ', 'CE2'): {"mean": 119.87, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CZ', 'CE2', 'CD2'): {"mean": 119.97, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard PHE.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "PRO": {
        # ---- 骨架 (Pro 细化: 111.8/116.9/122.6/122.0) ----
        ("N", "CA", "C"): {"mean": 111.8, "std": 2.5, "lb": 104.3, "up": 119.3,
                            "note": "Pro N-CH1E-C 111.8±2.5。", "source": "engh_huber_1991"},
        ("CA", "C", "N"): {"mean": 116.9, "std": 1.5, "lb": 112.4, "up": 121.4,
                            "note": "Pro CH1E-C-N 116.9±1.5。", "source": "engh_huber_1991"},
        ("CA", "C", "O"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "O")],
        ("C", "N", "CA"): {"mean": 122.6, "std": 5.0, "lb": 107.6, "up": 137.6,
                            "note": "Pro C-N-CH1E 122.6±5.0 (pyrrolidine 环构象浮动, sigma 大)。", "source": "engh_huber_1991"},
        ("O", "C", "N"): {"mean": 122.0, "std": 1.4, "lb": 117.8, "up": 126.2,
                            "note": "Pro O-C-N 122.0±1.4。", "source": "engh_huber_1991"},
        ("CB", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CB", "CA", "C")],
        # ---- 侧链 (Rosetta 408 ICOOR 理想点值) ----
        ('N', 'CA', 'CB'): {"mean": 103.00, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard PRO.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG'): {"mean": 104.20, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard PRO.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD'): {"mean": 104.40, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard PRO.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "SER": {
        # ---- 骨架 (通用值引用 generic; 本残基无细化) ----
        ("N", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("N", "CA", "C")],
        ("CA", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "N")],
        ("CA", "C", "O"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "O")],
        ("C", "N", "CA"): AMINO_ACID_BACKBONE_BOND_ANGLE[("C", "N", "CA")],
        ("O", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("O", "C", "N")],
        ("CB", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CB", "CA", "C")],
        # ---- 侧链 (Rosetta 408 ICOOR 理想点值) ----
        ('N', 'CA', 'CB'): {"mean": 110.14, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard SER.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'OG'): {"mean": 111.00, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard SER.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "THR": {
        # ---- 骨架 (通用值引用 generic; 本残基无细化) ----
        ("N", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("N", "CA", "C")],
        ("CA", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "N")],
        ("CA", "C", "O"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "O")],
        ("C", "N", "CA"): AMINO_ACID_BACKBONE_BOND_ANGLE[("C", "N", "CA")],
        ("O", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("O", "C", "N")],
        ("CB", "CA", "C"): {"mean": 109.1, "std": 2.2, "lb": 102.5, "up": 115.7,
                            "note": "Val/Ile/Thr CH1E-CH1E-C 109.1±2.2。", "source": "engh_huber_1991"},
        # ---- 侧链 (Rosetta 408 ICOOR 理想点值) ----
        ('N', 'CA', 'CB'): {"mean": 111.10, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard THR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'OG1'): {"mean": 109.58, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard THR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG2'): {"mean": 110.53, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard THR.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "TRP": {
        # ---- 骨架 (通用值引用 generic; 本残基无细化) ----
        ("N", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("N", "CA", "C")],
        ("CA", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "N")],
        ("CA", "C", "O"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "O")],
        ("C", "N", "CA"): AMINO_ACID_BACKBONE_BOND_ANGLE[("C", "N", "CA")],
        ("O", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("O", "C", "N")],
        ("CB", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CB", "CA", "C")],
        # ---- 侧链 (Rosetta 408 ICOOR 理想点值) ----
        ('N', 'CA', 'CB'): {"mean": 110.40, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG'): {"mean": 113.53, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD1'): {"mean": 126.74, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CD1', 'NE1'): {"mean": 110.16, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD1', 'NE1', 'CE2'): {"mean": 108.90, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('NE1', 'CE2', 'CZ2'): {"mean": 130.14, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CE2', 'CZ2', 'CH2'): {"mean": 117.50, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CZ2', 'CH2', 'CZ3'): {"mean": 121.50, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CH2', 'CZ3', 'CE3'): {"mean": 121.00, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CZ3', 'CE3', 'CD2'): {"mean": 118.73, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TRP.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "TYR": {
        # ---- 骨架 (通用值引用 generic; 本残基无细化) ----
        ("N", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("N", "CA", "C")],
        ("CA", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "N")],
        ("CA", "C", "O"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "O")],
        ("C", "N", "CA"): AMINO_ACID_BACKBONE_BOND_ANGLE[("C", "N", "CA")],
        ("O", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("O", "C", "N")],
        ("CB", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CB", "CA", "C")],
        # ---- 侧链 (Rosetta 408 ICOOR 理想点值) ----
        ('N', 'CA', 'CB'): {"mean": 110.47, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG'): {"mean": 113.80, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD1'): {"mean": 120.00, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CD1', 'CE1'): {"mean": 120.00, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CB', 'CG', 'CD2'): {"mean": 120.00, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CG', 'CD2', 'CE2'): {"mean": 120.00, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CD2', 'CE2', 'CZ'): {"mean": 120.00, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CE2', 'CZ', 'OH'): {"mean": 120.00, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard TYR.params ICOOR ideal", "source": "rosetta_params_408"},
    },
    "VAL": {
        # ---- 骨架 (通用值引用 generic; 本残基无细化) ----
        ("N", "CA", "C"): AMINO_ACID_BACKBONE_BOND_ANGLE[("N", "CA", "C")],
        ("CA", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "N")],
        ("CA", "C", "O"): AMINO_ACID_BACKBONE_BOND_ANGLE[("CA", "C", "O")],
        ("C", "N", "CA"): AMINO_ACID_BACKBONE_BOND_ANGLE[("C", "N", "CA")],
        ("O", "C", "N"): AMINO_ACID_BACKBONE_BOND_ANGLE[("O", "C", "N")],
        ("CB", "CA", "C"): {"mean": 109.1, "std": 2.2, "lb": 102.5, "up": 115.7,
                            "note": "Val/Ile/Thr CH1E-CH1E-C 109.1±2.2。", "source": "engh_huber_1991"},
        # ---- 侧链 (Rosetta 408 ICOOR 理想点值) ----
        ('N', 'CA', 'CB'): {"mean": 110.30, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard VAL.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG1'): {"mean": 110.49, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard VAL.params ICOOR ideal", "source": "rosetta_params_408"},
        ('CA', 'CB', 'CG2'): {"mean": 109.90, "std": np.nan, "lb": np.nan, "up": np.nan, "note": "Rosetta fa_standard VAL.params ICOOR ideal", "source": "rosetta_params_408"},
    },
}
