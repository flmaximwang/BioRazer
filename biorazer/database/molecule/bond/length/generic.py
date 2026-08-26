# -*- coding: utf-8 -*-
"""Molecule-agnostic covalent bond lengths (Å).

数据集来源与诚实性说明
────────────────────────
本模块记录**通用共价键长 (Å) 理想值**,用于结构构建/校验/生成时的
几何目标。数值与标准差**全部来自 Engh & Huber (1991) X 射线蛋白质结构
精修参数**,该套参数统计自剑桥晶体数据库 (CSD) 的小分子晶体结构:

* 主表 (含 sigma) 取自 PROCHECK 软件手册 Appendix A (Laskowski et al.,
  1993; 该表原封转载 Engh & Huber 1991 的 X-PLOR 标号体系数值)。
* 与 Rosetta 408 的氨基酸 .params (residue_type_sets/fa_standard/
  residue_types/l-caa/*.params) 内坐标 (ICOOR_INTERNAL) 交叉核对完成,
  **逐项一致** (见 :mod:`biorazer.database.molecule.bond.length.protein`
  的核对说明)。

这里的每条 key 是**两端原子名二元组** (如 ``("C", "N")``),值为最「一般」
的 Engh-Huber 条目 (即「除 Gly/Pro」或「除 Gly」那个)。蛋白质残基细分
值 (Gly/Pro/Ala/VIT...) 见 :mod:`.protein` 的 ``*_BY_RESIDUE`` 表。

单位: **埃 (Å)**。每条记录含 ``mean`` (目标值)、``std`` (CSD 样本
标准差, 可直接用作约束权重/容差)、``lb``/``up`` (由 ``mean ± 3*std``
得到的合理上下界, 供生成器/校验器作硬阈值)、``note``、``source``
(指向 :data:`BOND_REFS`)。查不到的字段为 ``np.nan``。

与 Rosetta 的核对
────────────────────────
Rosetta 408 主链 .params 用 ICOOR_INTERNAL 定义内坐标: ``d`` 列即键长。
对 ALA 解析得 N-CA 1.4580 / CA-C 1.5233 / C=O 1.2310 / 肽键 C-N 1.3287 /
CA-CB 1.5217, 与 Engh-Huber 主表各差 ≤0.006 Å —— 即 Rosetta 主链理想
几何直接源自 Engh & Huber 1991, 数值无实质差异。
"""

import numpy as np

#: 文献库 (provenance)。每条记录的 ``source`` 指向这里的键。
BOND_REFS = {
    "engh_huber_1991": "Engh RA, Huber R. Accurate bond and angle parameters for X-ray protein structure refinement. Acta Cryst. A47:392-400, 1991. (键长/键角主参数, 统计自 CSD 小分子晶体结构)",
    "procheck_appendix_a": "Laskowski RA, MacArthur MW, Moss DS, Thornton JM. PROCHECK: a program to check the stereochemical quality of protein structures. J. Appl. Cryst. 26:283-291, 1993. (手册 Appendix A.2 原封转载 Engh & Huber 1991 的 X-PLOR 标号键长/键角表)",
    "rosetta_params_408": "Rosetta 408 main/database/chemical/residue_type_sets/fa_standard/residue_types/l-caa/*.params (ICOOR_INTERNAL 主链内坐标), 与本模块主表交叉核对一致。",
}

# ---------------------------------------------------------------------------
# 键长主表 (单位: Å)
# ---------------------------------------------------------------------------
# 每个键以两端原子名二元组为 key。mean/std 为 Engh & Huber (1991) 目标值与
# CSD 样本标准差; lb/up = mean ± 3*std。
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
