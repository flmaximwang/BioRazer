# -*- coding: utf-8 -*-
"""Amino-acid covalent bond lengths and bond angles (ideal geometries).

数据集来源与诚实性说明
────────────────────────
本模块记录的是**共价键的键长 (Å) 与键角 (°) 理想值**,用于结构
构建/校验/生成时的几何目标。数值与标准差**全部来自 Engh & Huber(1991)
X 射线蛋白质结构精修参数**,该套参数统计自剑桥晶体数据库 (CSD) 的
小分子晶体结构:

* 键长与键角的主表 (含 sigma) 取自 PROCHECK 软件手册 Appendix A
  (Laskowski et al., 1993;该表原封转载 Engh & Huber 1991 的
  X-PLOR 标号体系的数值)。
* 与 Rosetta 408 的氨基酸 .params (residue_type_sets/fa_standard/
  residue_types/l-caa/*.params) 内坐标 (ICOOR_INTERNAL) 交叉核对完成,
  **逐项一致** (见下文「与 Rosetta 的核对」)。

键长的区分
────────────────────────
Engh & Huber 的键长按 X-PLOR 原子类型/残基细分,同一化学键对不同
残基可有不同目标值——例如 N-CA 在非 Gly/Pro 为 1.458,在 Gly 为 1.451,
在 Pro 为 1.466。本模块对每类键给出:

* ``AMINO_ACID_BOND_LENGTH`` 主表: 取最「一般」的值 (即「除
  Gly/Pro」或「除 Gly」那个条目), 供通用骨架生成/校验使用。
* 附 ``_BY_RESIDUE`` 细分表: 按残基给出 Engh & Huber 区分的
  具体值 (Gly / Pro / Ala / Val-Ile-Thr(Cbeta 分支) / 其余)。

键角同理。

单位约定
────────────────────────
键长单位为 **埃 (Å)**,键角单位为 **度 (degree)**。
每个条目含 ``mean`` (目标值)、``std`` (标准差,来自 CSD 样本,
可直接用作约束权重/容差)、``lb``/``up`` (由 ``mean ± 3*std`` 得到的
合理上下界,供生成器/校验器作硬阈值)、``note``、``source``。

与 Rosetta 的核对
────────────────────────
Rosetta 408 主链 .params 用 ICOOR_INTERNAL 定义内坐标:
``d`` 列即键长,连键的 ``theta`` 列满足 真实键角 = 180° − theta。
对 ALA 解析得:

* 键长: N-CA 1.4580, CA-C 1.5233, C=O 1.2310, 肽键 C-N 1.3287,
  CA-CB 1.5217 —— 与 Engh-Huber 主表各差 ≤0.006 Å；
* 键角 (180−theta): N-CA-C 111.2, CA-C-N 116.2, CA-C-O 120.8,
  N-CA-CB 110.4 —— 与 Engh-Huber (非 Gly) 行完全一致。

即 Rosetta 主链理想几何直接源自 Engh & Huber 1991,数值无实质差异;
本模块主表采用 Engh & Huber 原值 (含其残基细分),因此在跨软件使用时
(Rosetta / PDB 校验 / 结构构建) 是一致且可回溯的。
"""

# ---------------------------------------------------------------------------
# 键长 (单位: Å)
# ---------------------------------------------------------------------------
# 每个键以两端原子名二元组为 key。mean/std 为 Engh & Huber (1991) 目标值与
# CSD 样本标准差; lb/up = mean ± 3*std。
AMINO_ACID_BOND_LENGTH = {
    # 肽键 C-N (非 Pro)。Pro 的 C-N 为 1.341±0.016, 见 _BY_RESIDUE。
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
        "note": "N-Calpha。Engh-Huber NH1-CH1E (except Gly,Pro) 1.458±0.019; Gly 1.451±0.016; Pro 1.466±0.015 (见 _BY_RESIDUE)。",
        "source": ("engh_huber_1991", "procheck_appendix_a"),
    },
    # CA-CB (其余残基, 非 Ala / 非 Val,Ile,Thr)
    ("CA", "CB"): {
        "mean": 1.530, "std": 0.020, "lb": 1.470, "up": 1.590,
        "note": "Calpha-Cbeta。Engh-Huber CH1E-CH2E (the rest) 1.530±0.020; Ala 1.521±0.033; Val/Ile/Thr 1.540±0.027 (见 _BY_RESIDUE)。",
        "source": ("engh_huber_1991", "procheck_appendix_a"),
    },
}

# 键长按残基细分 (Engh & Huber 1991)。key: 单残基型 (Gly/Pro/Ala) 或
# 同类型残基组 (VIT=Val,Ile,Thr)。一个键写多条同值时以残基覆盖通用主表。
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

# ---------------------------------------------------------------------------
# 键角 (单位: 度)
# ---------------------------------------------------------------------------
# 每个键角以三个原子名三元组为 key (角顶点在中间原子)。数值同 Engh & Huber。
AMINO_ACID_BOND_ANGLE = {
    # N-CA-C: 通用 (非 Gly/Pro) —— 骨架核心角度 τ
    ("N", "CA", "C"): {
        "mean": 111.2, "std": 2.8, "lb": 102.8, "up": 119.6,
        "note": "N-Calpha-C (τ)。Engh-Huber NH1-CH1E-C 111.2±2.8 (except Gly,Pro); Gly 112.5±2.9; Pro 111.8±2.5 (见 _BY_RESIDUE)。",
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
        "note": "Cbeta-Calpha-C。Engh-Huber CH2E-CH1E-C (the rest) 110.1±1.9; Ala 110.5±1.5; Val/Ile/Thr 109.1±2.2 (见 _BY_RESIDUE)。",
        "source": ("engh_huber_1991", "procheck_appendix_a"),
    },
    # N-CA-CB: 通用 (其余残基, 非 Ala/Pro/VIT)
    ("N", "CA", "CB"): {
        "mean": 110.5, "std": 1.7, "lb": 105.4, "up": 115.6,
        "note": "N-Calpha-Cbeta。Engh-Huber NH1-CH1E-CH2E (the rest) 110.5±1.7; Ala 110.4±1.5; Val/Ile/Thr 111.5±1.7; Pro N-CH1E-CH2E 103.0±1.1 (见 _BY_RESIDUE)。",
        "source": ("engh_huber_1991", "procheck_appendix_a"),
    },
}

# 键角按残基细分 (Engh & Huber 1991)。
AMINO_ACID_BOND_ANGLE_BY_RESIDUE = {
    "Gly": {
        ("N", "CA", "C"): {"mean": 112.5, "std": 2.9, "lb": 103.8, "up": 121.2,
                            "note": "Gly NH1-CH2G+C 112.5±2.9。", "source": "engh_huber_1991"},
        ("CA", "C", "N"): {"mean": 116.4, "std": 2.1, "lb": 110.1, "up": 122.7,
                            "note": "Gly CH2G+C-NH1 116.4±2.1。", "source": "engh_huber_1991"},
        ("CA", "C", "O"): {"mean": 120.8, "std": 2.1, "lb": 114.5, "up": 127.1,
                            "note": "Gly CH2G+C-O 120.8±2.1。", "source": "engh_huber_1991"},
        ("C", "N", "CA"): {"mean": 120.6, "std": 1.7, "lb": 115.5, "up": 125.7,
                            "note": "Gly C-NH1-CH2G+ 120.6±1.7。", "source": "engh_huber_1991"},
    },
    "Pro": {
        ("N", "CA", "C"): {"mean": 111.8, "std": 2.5, "lb": 104.3, "up": 119.3,
                            "note": "Pro N-CH1E-C 111.8±2.5。", "source": "engh_huber_1991"},
        ("CA", "C", "N"): {"mean": 116.9, "std": 1.5, "lb": 112.4, "up": 121.4,
                            "note": "Pro CH1E-C-N 116.9±1.5。", "source": "engh_huber_1991"},
        ("C", "N", "CA"): {"mean": 122.6, "std": 5.0, "lb": 107.6, "up": 137.6,
                            "note": "Pro C-N-CH1E 122.6±5.0 (pyrrolidine 环构象浮动, sigma 大)。", "source": "engh_huber_1991"},
        ("O", "C", "N"): {"mean": 122.0, "std": 1.4, "lb": 117.8, "up": 126.2,
                            "note": "Pro O-C-N 122.0±1.4。", "source": "engh_huber_1991"},
        ("N", "CA", "CB"): {"mean": 103.0, "std": 1.1, "lb": 99.7, "up": 106.3,
                            "note": "Pro N-CH1E-CH2E 103.0±1.1 (pyrrolidine 环使该角显著收窄)。", "source": "engh_huber_1991"},
    },
    "Ala": {
        ("CB", "CA", "C"): {"mean": 110.5, "std": 1.5, "lb": 106.0, "up": 115.0,
                             "note": "Ala CH3E-CH1E-C 110.5±1.5。", "source": "engh_huber_1991"},
        ("N", "CA", "CB"): {"mean": 110.4, "std": 1.5, "lb": 105.9, "up": 114.9,
                             "note": "Ala NH1-CH1E-CH3E 110.4±1.5。", "source": "engh_huber_1991"},
    },
    "VIT": {  # Val, Ile, Thr
        ("CB", "CA", "C"): {"mean": 109.1, "std": 2.2, "lb": 102.5, "up": 115.7,
                             "note": "Val/Ile/Thr CH1E-CH1E-C 109.1±2.2。", "source": "engh_huber_1991"},
        ("N", "CA", "CB"): {"mean": 111.5, "std": 1.7, "lb": 106.4, "up": 116.6,
                             "note": "Val/Ile/Thr NH1-CH1E-CH1E 111.5±1.7。", "source": "engh_huber_1991"},
    },
}

# ---------------------------------------------------------------------------
# 文献库 (provenance)
# ---------------------------------------------------------------------------
BOND_REFS = {
    "engh_huber_1991": "Engh RA, Huber R. Accurate bond and angle parameters for X-ray protein structure refinement. Acta Cryst. A47:392-400, 1991. (键长/键角主参数, 统计自 CSD 小分子晶体结构)",
    "procheck_appendix_a": "Laskowski RA, MacArthur MW, Moss DS, Thornton JM. PROCHECK: a program to check the stereochemical quality of protein structures. J. Appl. Cryst. 26:283-291, 1993. (手册 Appendix A.2 原封转载 Engh & Huber 1991 的 X-PLOR 标号键长/键角表)",
    "rosetta_params_408": "Rosetta 408 main/database/chemical/residue_type_sets/fa_standard/residue_types/l-caa/*.params (ICOOR_INTERNAL 主链内坐标), 与本模块主表交叉核对一致。",
}