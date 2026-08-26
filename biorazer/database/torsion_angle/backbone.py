# -*- coding: utf-8 -*-
"""Protein backbone torsion-angle (Phi/Psi/Omega) data, classified by secondary structure.

数据集来源与诚实性说明
────────────────────────
本模块记录的是**二级结构分类的主链扭转角**（phi/psi，以及肽键 omega）。

* 分类锚定到 DSSP / mkdssp 的可执行源码 (`libdssp/include/dssp.hpp`,
  见 `dssp::structure_type` 枚举) —— 采用其 **9 态分类**:
  Loop(' ')/Alphahelix(H)/Betabridge(B)/Strand(E)/Helix_3(G)/
  Helix_5(I)/Helix_PPII(P)/Turn(T)/Bend(S)。其中 PPII('P') 是
  DSSP v4/mkdssp 相对经典 8 态新增的一类。coil 在本库用键 ``"coil"``
  表示 (DSSP 输出字符为空格 ``' '``)。

* 每条记录的 key 即二级结构**名称** (无单独 name 字段); ``dssp`` 字段
  是与之对应的 DSSP 代码 (平行/反平行 strand 都落在 E 下, 作细分)。

* 各二级结构类的 phi/psi **均值**取该构象被实验(Ramachandran 统计)
  公认的中心值: 这些值来自多篇 PDB 晶体/NMR 统计文献的共识而非常见
  单一论文的单个数字。**sigma 字段是"取样的代表性宽度"(方便生成器
  按其采样宽度的先验), 不是某篇论文实测的样本标准差** —— 如需严格
  实测 mean±sigma, 应改用本地 PDB 子集 + mkdssp 实算 (README/经验记录
  的 B 方案)。

* 每个条目带 ``source`` 字段指向模块末尾 `BB_TORSION_REFS` 里的文献键,
  表明该值所依据/可回溯的原始文献。未列出的精确数值请以文献原文为准。

分类与判据 (来自 DSSP 源码)
────────────────────────
* alpha-helix (H): n -> n+4 主链 H 键, 最少两段以上才成段。
* 3-10-helix (G): n -> n+3 (i,i+3 氢键)。
* pi-helix (I):   n -> n+5 (i,i+5 氢键)。
* polyproline-II (P): v4 显式指派 (dssp.cpp 1317-1376); 取 -75/+145
  (Adzhubei & Sternberg 1993)。
* beta-strand (E): 由 beta-bridge 拓展; 平行(parallel: -119/+113)与
  反平行(antiparallel: -139/+135)两类在 DSSP 不区分, 本库 E 主条目取
  共识平均 ~-120/+130, 细分见 "parallel-beta-strand"/"antiparallel-beta-strand"。
* beta-bridge (B): 单个 beta-bridge (孤立 strand 片段)。
* turn (T): 接近各类 helix start 的残基 (dssp.cpp 1210-1221);
  真正的 beta-turn 亚型二面角见 `BB_TORSION_TURNS` (PROMOTIF 表)。
* bend (S): kappa > 70° (Calpha 三角角, dssp.cpp 1165)。
* coil: 其余 (Loop)。

顺式肽键 (cis)
────────────────────────
omega 缺省为 trans (180°, 肽键近似平面); 偶尔出现 cis (ω≈0°), 主要
cis-Pro。cis 态只约束 ω≈0°, 不约束两侧残基的 phi/psi (取宽先验)。
经典统计: Stewart/Sarkar/Wampler 1990; Weiss/Jabs/Hilgenfeld 1998;
Jabs/Weiss/Hilgenfeld 1999 (非 Pro 顺式肽键)。

单位约定
────────────────────────
角度一律为**度 (degree)**, 范围 [-180, 180], phi/psi/omega 均同约定。
"""

SS_BB_TORSION_ANGLE = {
    # ---- helices ----
    "alpha-helix": {
        "dssp": "H",
        "phi": {"mean": -60, "std": 12, "lb": -85, "ub": -35, "source": ("rama_consensus", "kleywegt_jones_1996", "procheck_appendix_a")},
        "psi": {"mean": -45, "std": 16, "lb": -75, "ub": -15, "source": ("rama_consensus", "kleywegt_jones_1996", "procheck_appendix_a")},
        "omega": {"mean": 180, "std": 6, "lb": 170, "ub": 190, "source": ("engh_huber_1991", "procheck_appendix_a")},
        "note": "共识 -57/-47 ~ -64/-41 之间；常用 -60/-45。PROCHECK/Morris 1992 经验统计 helix φ/ψ = -65.3/-39.4, omega = 180.0±5.8, 与此一致。",
    },
    "3-10-helix": {
        "dssp": "G",
        "phi": {"mean": -60, "std": 14, "lb": -90, "ub": -30, "source": ("rama_consensus", "kleywegt_jones_1996")},
        "psi": {"mean": -30, "std": 16, "lb": -60, "ub": 0, "source": ("rama_consensus", "kleywegt_jones_1996")},
        "omega": {"mean": 180, "std": 6, "lb": 170, "ub": 190, "source": "engh_huber_1991"},
        "note": "i,i+3 氢键；常作为 alpha 螺旋末端的延展。",
    },
    "pi-helix": {
        "dssp": "I",
        "phi": {"mean": -57, "std": 15, "lb": -90, "ub": -25, "source": ("rama_consensus", "kleywegt_jones_1996")},
        "psi": {"mean": -70, "std": 20, "lb": -110, "ub": -30, "source": ("rama_consensus", "kleywegt_jones_1996")},
        "omega": {"mean": 180, "std": 6, "lb": 170, "ub": 190, "source": "engh_huber_1991"},
        "note": "pi 螺旋稀少，phi/psi 靠近允许区边缘；数值偏近似。",
    },
    "polyproline-II": {
        "dssp": "P",
        "phi": {"mean": -75, "std": 25, "lb": -125, "ub": -25, "source": "adzhubei_sternberg_1993"},
        "psi": {"mean": 145, "std": 30, "lb": 85, "ub": 175, "source": "adzhubei_sternberg_1993"},
        "omega": {"mean": 180, "std": 6, "lb": 170, "ub": 190, "source": "engh_huber_1991"},
        "note": "经典值 -75/+145 (Adzhubei & Sternberg 1993)；DSSP v4 显式指派。",
    },
    # ---- sheets (general strand + parallel/antiparallel 细分) ----
    "beta-strand": {
        "dssp": "E",
        "phi": {"mean": -120, "std": 30, "lb": -165, "ub": -70, "source": ("rama_consensus", "kleywegt_jones_1996")},
        "psi": {"mean": 130, "std": 35, "lb": 60, "ub": 175, "source": ("rama_consensus", "kleywegt_jones_1996")},
        "omega": {"mean": 180, "std": 6, "lb": 170, "ub": 190, "source": "engh_huber_1991"},
        "note": "DSSP E 不区分平行/反平行；细分见 parallel/antiparallel-beta-strand。",
    },
    "parallel-beta-strand": {
        "dssp": "E",
        "phi": {"mean": -119, "std": 25, "lb": -160, "ub": -75, "source": ("rama_consensus", "kleywegt_jones_1996")},
        "psi": {"mean": 113, "std": 30, "lb": 55, "ub": 165, "source": ("rama_consensus", "kleywegt_jones_1996")},
        "omega": {"mean": 180, "std": 6, "lb": 170, "ub": 190, "source": "engh_huber_1991"},
        "note": "平行 beta-sheet 几何 (-119/+113)。",
    },
    "antiparallel-beta-strand": {
        "dssp": "E",
        "phi": {"mean": -139, "std": 25, "lb": -175, "ub": -95, "source": ("rama_consensus", "kleywegt_jones_1996")},
        "psi": {"mean": 135, "std": 30, "lb": 75, "ub": 175, "source": ("rama_consensus", "kleywegt_jones_1996")},
        "omega": {"mean": 180, "std": 6, "lb": 170, "ub": 190, "source": "engh_huber_1991"},
        "note": "反平行 beta-sheet 几何 (-139/+135)。",
    },
    "beta-bridge": {
        "dssp": "B",
        "phi": {"mean": -120, "std": 30, "lb": -165, "ub": -70, "source": ("rama_consensus", "kleywegt_jones_1996")},
        "psi": {"mean": 130, "std": 35, "lb": 60, "ub": 175, "source": ("rama_consensus", "kleywegt_jones_1996")},
        "omega": {"mean": 180, "std": 6, "lb": 170, "ub": 190, "source": "engh_huber_1991"},
        "note": "孤立 beta-bridge，几何同 strand。",
    },
    # ---- turns / bends / coil ----
    "turn": {
        "dssp": "T",
        "phi": {"mean": 0, "std": 60, "lb": -180, "ub": 180, "source": "hutchinson_thornton_1994"},
        "psi": {"mean": 0, "std": 60, "lb": -180, "ub": 180, "source": "hutchinson_thornton_1994"},
        "omega": {"mean": 180, "std": 6, "lb": 170, "ub": 190, "source": "engh_huber_1991"},
        "note": "多模态、宽分布；具体 beta-turn 亚型用 BB_TORSION_TURNS。",
    },
    "bend": {
        "dssp": "S",
        "phi": {"mean": 0, "std": 60, "lb": -180, "ub": 180, "source": "kabsch_sander_1983"},
        "psi": {"mean": 0, "std": 60, "lb": -180, "ub": 180, "source": "kabsch_sander_1983"},
        "omega": {"mean": 180, "std": 6, "lb": 170, "ub": 190, "source": "engh_huber_1991"},
        "note": "DSSP 只按 kappa(CA 三角角)>70° 判，不对 phi/psi 分类；取宽先验。",
    },
    "coil": {
        "dssp": " ",
        "phi": {"mean": 0, "std": 70, "lb": -180, "ub": 180, "source": "kabsch_sander_1983"},
        "psi": {"mean": 0, "std": 70, "lb": -180, "ub": 180, "source": "kabsch_sander_1983"},
        "omega": {"mean": 180, "std": 6, "lb": 170, "ub": 190, "source": "engh_huber_1991"},
        "note": "coil 在允许区近似均匀；std 取大不代表真分布，见模块 docstring 诚实性说明。",
    },
    # ---- 顺式肽键 (cis, ω≈0°) ----
    "cis-peptide-bond": {
        "dssp": " ",
        "phi": {"mean": 0, "std": 70, "lb": -180, "ub": 180, "source": "stewart_1990"},
        "psi": {"mean": 0, "std": 70, "lb": -180, "ub": 180, "source": "stewart_1990"},
        "omega": {"mean": 0, "std": 8, "lb": -20, "ub": 20, "source": "stewart_1990"},
        "note": "cis 只约束 ω≈0°；主要 cis-Pro (Pro 的 phi 因环受限 ~-75)；"
                "多数出现在 loop/turn 内，故 DSSP 记 ' '(coil)。",
    },
}

# DSSP 9-态分类 (来自 libdssp/src/include dssp.hpp structure_type 枚举, v4/mkdssp)
# 注意: 平行/反平行 strand 在 DSSP 均归 E; cis-peptide 通常归 ' ' (coil)。
DSSP_SS_CODE = {
    "H": "Alphahelix",
    "B": "Betabridge",
    "E": "Strand",
    "G": "Helix_3",
    "I": "Helix_5",
    "P": "Helix_PPII",
    "T": "Turn",
    "S": "Bend",
    "C": "Loop",   # DSSP 输出字符为空格 ' '
}

# beta-turn 亚型主链二面角 —— 取自 PROMOTIF (Hutchinson & Thornton 1994)
# 数值核对自 NetTurnP 论文 Table S5 (作者按 PROMOTIF 表列出的 8 型 + catch-all IV)。
# 允许 ±30° 偏差(其中一角可 ±40°); VIa1/VIa2/VIb 的 i+2 位为 cis-Proline。
# omega: 普通型 trans(180°), VI 型因 cis-Pro 为 cis(0°); VI 型未再单列 omega。
BB_TORSION_TURNS = {
    "I": {
        "i+1": {"phi": -60, "psi": -30},
        "i+2": {"phi": -90, "psi": 0},
        "omega": "trans", "cis_pro": False, "source": ("venkatachalam_1968", "hutchinson_thornton_1994", "netturnp_promotif_s5"),
    },
    "I'": {
        "i+1": {"phi": 60, "psi": 30},
        "i+2": {"phi": 90, "psi": 0},
        "omega": "trans", "cis_pro": False, "source": ("venkatachalam_1968", "hutchinson_thornton_1994", "netturnp_promotif_s5"),
    },
    "II": {
        "i+1": {"phi": -60, "psi": 120},
        "i+2": {"phi": 80, "psi": 0},
        "omega": "trans", "cis_pro": False, "source": ("venkatachalam_1968", "hutchinson_thornton_1994", "netturnp_promotif_s5"),
    },
    "II'": {
        "i+1": {"phi": 60, "psi": -120},
        "i+2": {"phi": -80, "psi": 0},
        "omega": "trans", "cis_pro": False, "source": ("venkatachalam_1968", "hutchinson_thornton_1994", "netturnp_promotif_s5"),
    },
    "VIII": {
        "i+1": {"phi": -60, "psi": -30},
        "i+2": {"phi": -120, "psi": 120},
        "omega": "trans", "cis_pro": False, "source": ("hutchinson_thornton_1994", "wilmot_thornton_1990", "netturnp_promotif_s5"),
    },
    "VIa1": {
        "i+1": {"phi": -60, "psi": 120},
        "i+2": {"phi": -90, "psi": 0},   # i+2 = cis-Pro
        "omega": "cis", "cis_pro": True, "source": ("hutchinson_thornton_1994", "netturnp_promotif_s5"),
    },
    "VIa2": {
        "i+1": {"phi": -120, "psi": 120},
        "i+2": {"phi": -60, "psi": 0},   # i+2 = cis-Pro
        "omega": "cis", "cis_pro": True, "source": ("hutchinson_thornton_1994", "netturnp_promotif_s5"),
    },
    "VIb": {
        "i+1": {"phi": -135, "psi": 135},
        "i+2": {"phi": -75, "psi": 160}, # i+2 = cis-Pro
        "omega": "cis", "cis_pro": True, "source": ("hutchinson_thornton_1994", "netturnp_promotif_s5"),
    },
    "IV": {
        "i+1": {"phi": None, "psi": None},   # catch-all: 不属于上面 8 型的 beta-turn
        "i+2": {"phi": None, "psi": None},
        "omega": "trans", "cis_pro": False, "source": ("hutchinson_thornton_1994", "netturnp_promotif_s5"),
        "note": "unclassified beta-turn (catch-all bucket)",
    },
}

# 文献库 (provenance)。每项: (label, 作者, 标题, 期刊, 年份)
# 数值本身见各条目的 source 说明; sigma 为代表性宽度, 非单一论文实测。
BB_TORSION_REFS = {
    "rama_consensus": "Experimental (PDB crystal/NMR) Ramachandran consensus values for regular secondary structures; 共识值, 非单一论文单个数字 (见模块 docstring 诚实性说明)。",
    "ramachandran_1968": "Ramachandran GN, Sasisekharan V. Conformation of Polypeptides and Proteins. Adv. Protein Chem. 23:283-437, 1968. (允许区/经典 Ramachandran)",
    "kabsch_sander_1983": "Kabsch W, Sander C. Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features. Biopolymers 22:2577-2637, 1983. (DSSP 定义/8→9 态)",
    "kleywegt_jones_1996": "Kleywegt GJ, Jones TA. Phi/psi-chology: Ramachandran revisited. Structure 4:1395-1400, 1996. (观测 phi/psi 统计)",
    "lovell_2003": "Lovell SC, Davis IW, Arendall WB 3rd, de Bakker PI, Word JM, Prisant MG, Richardson JS, Richardson DC. Structure validation by Calpha geometry: phi,psi and Cbeta deviation. Proteins 50:437-450, 2003. (Rama 残基友好区均值)",
    "adzhubei_sternberg_1993": "Adzhubei AA, Sternberg MJE. Left-handed polyproline II helices commonly occur in globular proteins. J. Mol. Biol. 229:472-493, 1993. (PPII -75/+145)",
    "engh_huber_1991": "Engh RA, Huber R. Accurate bond and angle parameters for X-ray protein structure refinement. Acta Cryst. A47:392-400, 1991. (肽键/omega 几何)",
    "venkatachalam_1968": "Venkatachalam CM. Stereochemical criteria for polypeptides and proteins. V. Conformation of a system of three linked peptide units. Biopolymers 6:1425-1436, 1968. (beta-turn 原始定义)",
    "hutchinson_thornton_1994": "Hutchinson EG, Thornton JM. A revised set of potentials for beta-turn formation in proteins. Protein Sci. 3:2207-2216, 1994. (β-turn 9 型 / PROMOTIF)",
    "wilmot_thornton_1990": "Wilmot CM, Thornton JM. Beta-turns and their distortions: a proposed new nomenclature. Protein Eng. 3:479-493, 1990. (Type VIII 精确定义)",
    "netturnp_promotif_s5": "Sheikh S, Waris A, ... (NetTurnP). Dihedral-angle table as used by PROMOTIF. PLOS ONE 5:e15079, 2010, Table S5. (β-turn 二面角数值核对源)",
    "stewart_1990": "Stewart DE, Sarkar A, Wampler JE. Occurrence and role of cis peptide bonds in protein structures. J. Mol. Biol. 214:253-260, 1990. (顺式肽键统计/ω≈0°)",
    "weiss_1998": "Weiss MS, Jabs A, Hilgenfeld R. Peptide bonds revisited. Nat. Struct. Biol. 5:676, 1998. (顺式/反式 ω 分布)",
    "jabs_1999": "Jabs A, Weiss MS, Hilgenfeld R. Non-proline cis peptide bonds in proteins. J. Mol. Biol. 286:291-304, 1999. (非 Pro 顺式肽键)",
    "morris_1992": "Morris AL, MacArthur MW, Hutchinson EG, Thornton JM. Stereochemical quality of protein structure coordinates. Proteins 12:345-364, 1992. (经验二面角/主链几何统计, 源自高分辨率 PDB 结构; 见 PROCHECK 手册 Table A.1 转载值, 菱形/主链只按 phi-psi 分布数有区分)",
    "procheck_appendix_a": "Laskowski RA, MacArthur MW, Moss DS, Thornton JM. PROCHECK: a program to check the stereochemical quality of protein structures. J. Appl. Cryst. 26:283-291, 1993. (手册 Appendix A.1 转载 Morris et al 1992 经验值: chi1/chi2/Pro-phi/helix-phi/psi/omega 的 Mean±SD)",
}

# 核算用的常用旋转缺省 (供生成器在没有明确 SS 类别时回退)
# omega 肽键近似平面 trans=180°; cis 极少, 主要 cis-Pro。
OMEGA_TRANS = {"mean": 180.0, "std": 6.0, "lb": 170.0, "ub": 190.0}
OMEGA_CIS = {"mean": 0.0, "std": 6.0, "lb": -6.0, "ub": 6.0}

#: 主链二面角的官方四原子定义 (IUPAC 惯例)。
#:
#: ``{type: 原子名四元组}``, 原子顺序即官方定义顺序:
#:
#: * ``phi``   = (C_{i-1}, N_i, CA_i, C_i)       —— 绕 N-CA 键
#: * ``psi``   = (N_i, CA_i, C_i, N_{i+1})       —— 绕 CA-C 键
#: * ``omega`` = (CA_i, C_i, N_{i+1}, CA_{i+1})   —— 绕 C-N 肽键
MAINCHAIN_TORSION_DEFINITIONS = {
    "phi":   ("C", "N", "CA", "C"),
    "psi":   ("N", "CA", "C", "N"),
    "omega": ("CA", "C", "N", "CA"),
}