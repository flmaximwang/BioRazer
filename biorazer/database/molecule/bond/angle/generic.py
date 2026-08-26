# -*- coding: utf-8 -*-
"""Molecule-agnostic covalent bond angles (°).

数据集来源与诚实性说明
────────────────────────
本模块记录**通用共价键角 (°) 理想值**,用于结构构建/校验/生成时的
几何目标。数值与标准差**全部来自 Engh & Huber (1991) X 射线蛋白质结构
精修参数**,该套参数统计自剑桥晶体数据库 (CSD) 的小分子晶体结构:

* 主表 (含 sigma) 取自 PROCHECK 软件手册 Appendix A (Laskowski et al.,
  1993; 该表原封转载 Engh & Huber 1991 的 X-PLOR 标号体系数值)。
* 与 Rosetta 408 的氨基酸 .params (residue_type_sets/fa_standard/
  residue_types/l-caa/*.params) 内坐标 (ICOOR_INTERNAL) 交叉核对完成,
  **逐项一致** (见 :mod:`biorazer.database.molecule.bond.angle.protein`
  的核对说明)。

这里的每条 key 是**三个原子名三元组** (如 ``("N", "CA", "C")``, 角顶点
在中间原子), 值为最「一般」的 Engh-Huber 条目 (即「除 Gly/Pro」或
「除 Gly」那个)。蛋白质残基细分值 (Gly/Pro/Ala/VIT...) 见
:mod:`.protein` 的 ``*_BY_RESIDUE`` 表。

单位: **度 (degree)**。每条记录含 ``mean`` (目标值)、``std`` (CSD 样本
标准差)、``lb``/``up`` (由 ``mean ± 3*std`` 得到的合理上下界)、
``note``、``source`` (指向 :data:`BOND_REFS`)。查不到的字段为 ``np.nan``。

与 Rosetta 的核对
────────────────────────
Rosetta 408 主链 .params 用 ICOOR_INTERNAL 定义内坐标: 连键的 ``theta``
列满足 真实键角 = 180° − theta。对 ALA 解析得 N-CA-C 111.2 / CA-C-N 116.2 /
CA-C-O 120.8 / N-CA-CB 110.4 —— 与 Engh-Huber (非 Gly) 行完全一致。
"""

import numpy as np

from ..length.generic import BOND_REFS  # noqa: F401  (shared provenance)

# ---------------------------------------------------------------------------
# 键角主表 (单位: 度)
# ---------------------------------------------------------------------------
# 每个键角以三个原子名三元组为 key (角顶点在中间原子)。数值同 Engh & Huber。
AMINO_ACID_BOND_ANGLE = {
    # N-CA-C: 通用 (非 Gly/Pro) —— 骨架核心角度 τ
    ("N", "CA", "C"): {
        "mean": 111.2, "std": 2.8, "lb": 102.8, "up": 119.6,
        "note": "N-Calpha-C (τ)。Engh-Huber NH1-CH1E-C 111.2±2.8 (except Gly,Pro); Gly 112.5±2.9; Pro 111.8±2.5 (见 protein._BY_RESIDUE)。",
        "source": ("engh_huber_1991", "procheck_appendix_a"),
    },
    # CA-C-N: 通用 (非 Gly/Pro)
    ("CA", "C", "N"): {
        "mean": 116.2, "std": 2.0, "lb": 110.2, "up": 122.2,
        "note": "Calpha-C-N。Engh-Huber CH1E-C-NH1 116.2±2.0 (except Gly,Pro); Gly 116.4±2.1; Pro 116.9±1.5。",
        "source": ("engh_huber_1991", "procheck_appendix_a"),
    },
    # CA-C-O: 通用 (非 Gly)
    ("CA", "C", "O"): {
        "mean": 120.8, "std": 1.7, "lb": 115.7, "up": 125.9,
        "note": "Calpha-C-O。Engh-Huber CH1E-C-O 120.8±1.7 (except Gly); Gly 120.8±2.1。",
        "source": ("engh_huber_1991", "procheck_appendix_a"),
    },
    # C-N-CA: 通用 (非 Gly/Pro)
    ("C", "N", "CA"): {
        "mean": 121.7, "std": 1.8, "lb": 116.3, "up": 127.1,
        "note": "C-N-Calpha。Engh-Huber C-NH1-CH1E 121.7±1.8 (except Gly,Pro); Gly 120.6±1.7; Pro 122.6±5.0。",
        "source": ("engh_huber_1991", "procheck_appendix_a"),
    },
    # O-C-N: 通用 (非 Pro)
    ("O", "C", "N"): {
        "mean": 123.0, "std": 1.6, "lb": 118.2, "up": 127.8,
        "note": "O-C-N。Engh-Huber O-C-NH1 123.0±1.6 (except Pro); Pro 122.0±1.4。",
        "source": ("engh_huber_1991", "procheck_appendix_a"),
    },
    # CB-CA-C: 通用 (其余残基)
    ("CB", "CA", "C"): {
        "mean": 110.1, "std": 1.9, "lb": 104.4, "up": 115.8,
        "note": "Cbeta-Calpha-C。Engh-Huber CH2E-CH1E-C (the rest) 110.1±1.9; Ala 110.5±1.5; Val/Ile/Thr 109.1±2.2 (见 protein._BY_RESIDUE)。",
        "source": ("engh_huber_1991", "procheck_appendix_a"),
    },
    # N-CA-CB: 通用 (其余残基, 非 Ala/Pro/VIT)
    ("N", "CA", "CB"): {
        "mean": 110.5, "std": 1.7, "lb": 105.4, "up": 115.6,
        "note": "N-Calpha-Cbeta。Engh-Huber NH1-CH1E-CH2E (the rest) 110.5±1.7; Ala 110.4±1.5; Val/Ile/Thr 111.5±1.7; Pro N-CH1E-CH2E 103.0±1.1 (见 protein._BY_RESIDUE)。",
        "source": ("engh_huber_1991", "procheck_appendix_a"),
    },
}
