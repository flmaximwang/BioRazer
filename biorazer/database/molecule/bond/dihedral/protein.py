# -*- coding: utf-8 -*-
"""Protein torsion (dihedral) data: backbone phi/psi/omega by secondary
structure, beta-turn dihedrals, chi definitions and rotamer framework.

数据集来源与诚实性说明
────────────────────────
本模块合并了原 ``torsion_angle.backbone`` (主链) 与
``torsion_angle.sidechain`` (侧链) 的扭转角数据:

主链部分
~~~~~~~~
* 分类锚定到 DSSP / mkdssp 的可执行源码 (`libdssp/include/dssp.hpp`,
  见 `dssp::structure_type` 枚举) —— 采用其 **9 态分类**:
  Loop(' ')/Alphahelix(H)/Betabridge(B)/Strand(E)/Helix_3(G)/
  Helix_5(I)/Helix_PPII(P)/Turn(T)/Bend(S)。其中 PPII('P') 是
  DSSP v4/mkdssp 相对经典 8 态新增的一类。coil 在本库用键 ``"coil"``
  表示 (DSSP 输出字符为空格 ``' '``)。
* 每条记录的 key 即二级结构**名称**; ``dssp`` 字段是与之对应的 DSSP
  代码 (平行/反平行 strand 都落在 E 下, 作细分)。
* 各二级结构类的 phi/psi **均值**取该构象被实验 (Ramachandran 统计)
  公认的中心值: 来自多篇 PDB 晶体/NMR 统计文献的共识而非常见单一论文的
  单个数字。**sigma 字段是\"取样的代表性宽度\"(方便生成器按其采样宽度
  的先验), 不是某篇论文实测的样本标准差** —— 如需严格实测 mean±sigma,
  应改用本地 PDB 子集 + mkdssp 实算 (B 方案)。
* 每个条目带 ``source`` 字段指向 :data:`BB_TORSION_REFS` 里的文献键。

侧链部分
~~~~~~~~
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
* ``ROTAMER_BIN`` / ``DUNBRACK_ROTAMERS`` -- 侧链 rotamer 的**分类框架**:
  标准 rotamer bin 中心 (g-/t/g+ = -60/180/+60) 和每种残基的可旋转 chi
  轴数。完整的**骨架依赖数值表** (逐 phi/psi bin 的均值/方差) 属于外部
  Dunbrack 2010 数据集 (Shapovalov & Dunbrack 2011, CC BY 4.0), 未内嵌,
  需要时可另行 vendor。

所有角度单位 **度 (degree)**。数值记录统一为 ``{mean, std, lb, up,
source}``; 查不到 spread 的字段为 ``np.nan``。
"""

import numpy as np

# ---------------------------------------------------------------------------
# 主链: 二级结构分类的 phi/psi/omega
# ---------------------------------------------------------------------------
# 每条记录: {dssp, phi: {...}, psi: {...}, omega: {...}, note}
SS_BB_TORSION_ANGLE = {
    # ---- helices ----
    "alpha-helix": {
        "dssp": "H",
        "phi": {"mean": -60, "std": 12, "lb": -85, "up": -35, "source": ("rama_consensus", "kleywegt_jones_1996", "procheck_appendix_a")},
        "psi": {"mean": -45, "std": 16, "lb": -75, "up": -15, "source": ("rama_consensus", "kleywegt_jones_1996", "procheck_appendix_a")},
        "omega": {"mean": 180, "std": 6, "lb": 170, "up": 190, "source": ("engh_huber_1991", "procheck_appendix_a")},
        "note": "共识 -57/-47 ~ -64/-41 之间；常用 -60/-45。PROCHECK/Morris 1992 经验统计 helix φ/ψ = -65.3/-39.4, omega = 180.0±5.8, 与此一致。",
    },
    "3-10-helix": {
        "dssp": "G",
        "phi": {"mean": -60, "std": 14, "lb": -90, "up": -30, "source": ("rama_consensus", "kleywegt_jones_1996")},
        "psi": {"mean": -30, "std": 16, "lb": -60, "up": 0, "source": ("rama_consensus", "kleywegt_jones_1996")},
        "omega": {"mean": 180, "std": 6, "lb": 170, "up": 190, "source": "engh_huber_1991"},
        "note": "i,i+3 氢键；常作为 alpha 螺旋末端的延展。",
    },
    "pi-helix": {
        "dssp": "I",
        "phi": {"mean": -57, "std": 15, "lb": -90, "up": -25, "source": ("rama_consensus", "kleywegt_jones_1996")},
        "psi": {"mean": -70, "std": 20, "lb": -110, "up": -30, "source": ("rama_consensus", "kleywegt_jones_1996")},
        "omega": {"mean": 180, "std": 6, "lb": 170, "up": 190, "source": "engh_huber_1991"},
        "note": "pi 螺旋稀少，phi/psi 靠近允许区边缘；数值偏近似。",
    },
    "polyproline-II": {
        "dssp": "P",
        "phi": {"mean": -75, "std": 25, "lb": -125, "up": -25, "source": "adzhubei_sternberg_1993"},
        "psi": {"mean": 145, "std": 30, "lb": 85, "up": 175, "source": "adzhubei_sternberg_1993"},
        "omega": {"mean": 180, "std": 6, "lb": 170, "up": 190, "source": "engh_huber_1991"},
        "note": "经典值 -75/+145 (Adzhubei & Sternberg 1993)；DSSP v4 显式指派。",
    },
    # ---- sheets (general strand + parallel/antiparallel 细分) ----
    "beta-strand": {
        "dssp": "E",
        "phi": {"mean": -120, "std": 30, "lb": -165, "up": -70, "source": ("rama_consensus", "kleywegt_jones_1996")},
        "psi": {"mean": 130, "std": 35, "lb": 60, "up": 175, "source": ("rama_consensus", "kleywegt_jones_1996")},
        "omega": {"mean": 180, "std": 6, "lb": 170, "up": 190, "source": "engh_huber_1991"},
        "note": "DSSP E 不区分平行/反平行；细分见 parallel/antiparallel-beta-strand。",
    },
    "parallel-beta-strand": {
        "dssp": "E",
        "phi": {"mean": -119, "std": 25, "lb": -160, "up": -75, "source": ("rama_consensus", "kleywegt_jones_1996")},
        "psi": {"mean": 113, "std": 30, "lb": 55, "up": 165, "source": ("rama_consensus", "kleywegt_jones_1996")},
        "omega": {"mean": 180, "std": 6, "lb": 170, "up": 190, "source": "engh_huber_1991"},
        "note": "平行 beta-sheet 几何 (-119/+113)。",
    },
    "antiparallel-beta-strand": {
        "dssp": "E",
        "phi": {"mean": -139, "std": 25, "lb": -175, "up": -95, "source": ("rama_consensus", "kleywegt_jones_1996")},
        "psi": {"mean": 135, "std": 30, "lb": 75, "up": 175, "source": ("rama_consensus", "kleywegt_jones_1996")},
        "omega": {"mean": 180, "std": 6, "lb": 170, "up": 190, "source": "engh_huber_1991"},
        "note": "反平行 beta-sheet 几何 (-139/+135)。",
    },
    "beta-bridge": {
        "dssp": "B",
        "phi": {"mean": -120, "std": 30, "lb": -165, "up": -70, "source": ("rama_consensus", "kleywegt_jones_1996")},
        "psi": {"mean": 130, "std": 35, "lb": 60, "up": 175, "source": ("rama_consensus", "kleywegt_jones_1996")},
        "omega": {"mean": 180, "std": 6, "lb": 170, "up": 190, "source": "engh_huber_1991"},
        "note": "孤立 beta-bridge，几何同 strand。",
    },
    # ---- turns / bends / coil ----
    "turn": {
        "dssp": "T",
        "phi": {"mean": 0, "std": 60, "lb": -180, "up": 180, "source": "hutchinson_thornton_1994"},
        "psi": {"mean": 0, "std": 60, "lb": -180, "up": 180, "source": "hutchinson_thornton_1994"},
        "omega": {"mean": 180, "std": 6, "lb": 170, "up": 190, "source": "engh_huber_1991"},
        "note": "多模态、宽分布；具体 beta-turn 亚型用 BB_TORSION_TURNS。",
    },
    "bend": {
        "dssp": "S",
        "phi": {"mean": 0, "std": 60, "lb": -180, "up": 180, "source": "kabsch_sander_1983"},
        "psi": {"mean": 0, "std": 60, "lb": -180, "up": 180, "source": "kabsch_sander_1983"},
        "omega": {"mean": 180, "std": 6, "lb": 170, "up": 190, "source": "engh_huber_1991"},
        "note": "DSSP 只按 kappa(CA 三角角)>70° 判，不对 phi/psi 分类；取宽先验。",
    },
    "coil": {
        "dssp": " ",
        "phi": {"mean": 0, "std": 70, "lb": -180, "up": 180, "source": "kabsch_sander_1983"},
        "psi": {"mean": 0, "std": 70, "lb": -180, "up": 180, "source": "kabsch_sander_1983"},
        "omega": {"mean": 180, "std": 6, "lb": 170, "up": 190, "source": "engh_huber_1991"},
        "note": "coil 在允许区近似均匀；std 取大不代表真分布，见模块 docstring 诚实性说明。",
    },
    # ---- 顺式肽键 (cis, ω≈0°) ----
    "cis-peptide-bond": {
        "dssp": " ",
        "phi": {"mean": 0, "std": 70, "lb": -180, "up": 180, "source": "stewart_1990"},
        "psi": {"mean": 0, "std": 70, "lb": -180, "up": 180, "source": "stewart_1990"},
        "omega": {"mean": 0, "std": 8, "lb": -20, "up": 20, "source": "stewart_1990"},
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
# phi/psi 为 {mean, std, lb, up, source} 记录; 无 spread 时 std/lb/up = np.nan。
BB_TORSION_TURNS = {
    "I": {
        "i+1": {"mean_phi": -60, "mean_psi": -30},
        "i+2": {"mean_phi": -90, "mean_psi": 0},
        "omega": "trans", "cis_pro": False, "source": ("venkatachalam_1968", "hutchinson_thornton_1994", "netturnp_promotif_s5"),
    },
    "I'": {
        "i+1": {"mean_phi": 60, "mean_psi": 30},
        "i+2": {"mean_phi": 90, "mean_psi": 0},
        "omega": "trans", "cis_pro": False, "source": ("venkatachalam_1968", "hutchinson_thornton_1994", "netturnp_promotif_s5"),
    },
    "II": {
        "i+1": {"mean_phi": -60, "mean_psi": 120},
        "i+2": {"mean_phi": 80, "mean_psi": 0},
        "omega": "trans", "cis_pro": False, "source": ("venkatachalam_1968", "hutchinson_thornton_1994", "netturnp_promotif_s5"),
    },
    "II'": {
        "i+1": {"mean_phi": 60, "mean_psi": -120},
        "i+2": {"mean_phi": -80, "mean_psi": 0},
        "omega": "trans", "cis_pro": False, "source": ("venkatachalam_1968", "hutchinson_thornton_1994", "netturnp_promotif_s5"),
    },
    "VIII": {
        "i+1": {"mean_phi": -60, "mean_psi": -30},
        "i+2": {"mean_phi": -120, "mean_psi": 120},
        "omega": "trans", "cis_pro": False, "source": ("hutchinson_thornton_1994", "wilmot_thornton_1990", "netturnp_promotif_s5"),
    },
    "VIa1": {
        "i+1": {"mean_phi": -60, "mean_psi": 120},
        "i+2": {"mean_phi": -90, "mean_psi": 0},   # i+2 = cis-Pro
        "omega": "cis", "cis_pro": True, "source": ("hutchinson_thornton_1994", "netturnp_promotif_s5"),
    },
    "VIa2": {
        "i+1": {"mean_phi": -120, "mean_psi": 120},
        "i+2": {"mean_phi": -60, "mean_psi": 0},   # i+2 = cis-Pro
        "omega": "cis", "cis_pro": True, "source": ("hutchinson_thornton_1994", "netturnp_promotif_s5"),
    },
    "VIb": {
        "i+1": {"mean_phi": -135, "mean_psi": 135},
        "i+2": {"mean_phi": -75, "mean_psi": 160}, # i+2 = cis-Pro
        "omega": "cis", "cis_pro": True, "source": ("hutchinson_thornton_1994", "netturnp_promotif_s5"),
    },
    "IV": {
        "i+1": {"mean_phi": None, "mean_psi": None},   # catch-all: 不属于上面 8 型的 beta-turn
        "i+2": {"mean_phi": None, "mean_psi": None},
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
OMEGA_TRANS = {"mean": 180.0, "std": 6.0, "lb": 170.0, "up": 190.0, "source": "engh_huber_1991"}
OMEGA_CIS = {"mean": 0.0, "std": 6.0, "lb": -6.0, "up": 6.0, "source": "stewart_1990"}

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

# ---------------------------------------------------------------------------
# 侧链: chi 定义 + 规范 IC-frame 二面角 + rotamer 框架
# ---------------------------------------------------------------------------

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
#: (与 :data:`~biorazer.database.molecule.bond.length.protein.AMINO_ACID_SIDECHAIN_BOND`
#: 同 key)。GLY 侧链为空。每条为 {mean, std, lb, up, source}; Rosetta ICOOR
#: 只给理想点值, 故 std/lb/up = np.nan。
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

#: per 残基的 rotamer 分类框架: 可旋转 chi 轴数 (见 SIDECHAIN_CHI) + 标准 bin 中心。
#: 完整数值表 (逐 phi/psi) 未内嵌。
DUNBRACK_ROTAMERS = {
    'ALA': {"chi": 0, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'ARG': {"chi": 4, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'ASN': {"chi": 2, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'ASP': {"chi": 2, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'CYS': {"chi": 1, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'GLN': {"chi": 3, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'GLU': {"chi": 3, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'GLY': {"chi": 0, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'HIS': {"chi": 2, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'ILE': {"chi": 2, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'LEU': {"chi": 2, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'LYS': {"chi": 4, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'MET': {"chi": 3, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'PHE': {"chi": 2, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'PRO': {"chi": 2, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'SER': {"chi": 1, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'THR': {"chi": 1, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'TRP': {"chi": 2, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'TYR': {"chi": 2, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
    'VAL': {"chi": 1, "bins": (("g-", -60.0), ("t", 180.0), ("g+", 60.0)), "note": "rotamer 分类框架; 数值表未内嵌 (见模块 docstring)"},
}

SIDECHAIN_DIHE_REFS = {
    "rosetta_params_408": "Rosetta 408 ... l-caa/*.params (ICOOR_INTERNAL 规范残基理想几何)",
    "dunbrack_2010": "Shapovalov MV, Dunbrack RL Jr. A smoothed backbone-dependent rotamer library. Structure 19:844-858, 2011."
}
