# -*- coding: utf-8 -*-
"""原子选择器: 规则表 / 选择表的语法、匹配, 以及它的 GUI 编辑器。

**两个 csv** (规则表与选择表是同一套语法的两种用法):

1. **规则表** (:meth:`AtomArraySelector.from_csv` 的 ``mode="rule"``): 每行一条**选择式**,
   5 列 ``ins_code``/``chain``/``resi``/``name``/``altloc``, 每格是一条**模式**
   (通配符/正则/范围/列表), 一行展开成一条**原子序列**。
2. **选择表** (``mode="selection"``): **严格每行一个原子**, 字段列是具体值, 由规则表按目标
   展开而来。下游只认这个文件。

具体值 (字面值) 展开后就是它自己一个原子, 所以**选择表可以当规则表打开** —— 拿它反查
"这些原子在结构里都在吗", 或把某一行改成范围批量扩。

| 模式 | 框里怎么填 | 例 |
|---|---|---|
| 任意 | 不用填 | 该字段不参与筛选 |
| 字面值 | 精确值 | ``A`` / ``CA`` / ``45``; 含 ``,`` ``-`` ``*`` ``?`` ``[]`` 时自动转义成正则 |
| 通配符 | fnmatch | ``A*`` / ``C*`` / ``?`` |
| 列表 | 逗号分隔 | ``1,5-9`` / ``CA,CB`` |
| 范围 | 起止, 含两端 | ``1-10`` / ``A-D`` / ``H1-H20`` |
| 正则 | ``re:`` 前缀 | ``re:[AB]`` (全匹配; 整格一条正则, 不按逗号拆) |

``ins_code`` 与 ``altloc`` 两个特例: 字面值 + 空 = 只命中**没有**插入码 / **没有** altloc 的原子
(``,A,45,CA`` 不会连 ``45A`` 一起命中); 任意 (空框, 写成 ``*``) = 管它有没有。字面值**区分大小写**。

``altloc`` 的值来自 biotite 的 ``altloc_id`` 标注, 只有按 ``altloc="all"`` 读进来的数组才有
(:class:`biorazer.structure.io.protein.StructureFile_AtomArray` 就是那么读的)。
:class:`~biorazer.structure.objects.InternalCoord` **不带** altloc (它的记录一个原子名一格, 见
:class:`~biorazer.structure.bridge.AtomArray_InternalCoord`), 所以 IC 目标与手工拼的
``AtomArray`` 一样属于"没有这个 category": 带 altloc 约束的规则会报"无法判定"而不是静默当命中。
具体见 :class:`_TargetView`。

用法::

    sel = AtomArraySelector.from_csv("rules.csv", mode="rule")   # 或自己造 rules
    mask = sel.mask(atom_array)              # 与 atom_array 对齐的布尔 mask
    idx = sel.indices(atom_array)            # 展平后的原子下标 (规则序 → 结构序)
    sub = sel.apply(atom_array)              # 选择后的 AtomArray
    sub = sel.apply(internal_coord)          # 选择后的 InternalCoord (见 apply 的约束)
    sel.to_csv("selection.csv", mode="selection", structure=atom_array)
    sel.run_editor(structure_file="x.pdb", rule_csv="rules.csv", selection_csv="out.csv")

命令行 (parser/runner 就在本模块, ``objects/cli.py`` 只做注册)::

    biorazer select -s x.pdb --rules rules.csv -o selection.csv --print --export
    biorazer select -s x.pdb --rules rules.csv          # 不开 --print 则开窗编辑

选择器吃得下 ``AtomArray`` 与 ``InternalCoord`` 两种目标: 两者暴露同样的
``ins_code`` / ``chain_id`` / ``res_id`` / ``atom_name`` / ``res_name``, 匹配一律按
PDB/auth 口径 (biotite 的 ``use_author_fields`` 默认值)。

规则表/选择表的 csv 由本模块读写; 文件 → 内存对象的读取走
:mod:`biorazer.structure.io` (``StructureFile_AtomArray``); 变成 mask / indices 的**转换器**
(``AtomArraySelector_AtomArrayMask`` / ``AtomArraySelector_AtomArrayIndices``) 在
:mod:`biorazer.structure.bridge.selector`; GUI 在同包的私有模块 ``_selector_gui``。

已知边界: 五元组 (ins_code, chain, resi, name, altloc) 不是全局唯一 —— 同一残基上多个原子
重名、或目标不带 altloc 标注时会出现同 key 的多个原子 (命中重复会在警告里报个数)。展平**默认
去重**: 同一个五元组只写第一次命中的那条规则 (``dedupe=False`` 则重复照写)。重复个数任何时候
都报出来, 不悄悄吞掉。
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter
from dataclasses import dataclass, field
from fnmatch import fnmatchcase
from pathlib import Path

import numpy as np

from biorazer.structure.objects.bt_atom_array import AtomArrayStack
from biorazer.structure.objects.internal_coords import InternalCoord

FIELDS: tuple[str, ...] = ("ins_code", "chain", "resi", "name", "altloc")

#: biotite spells "no alternate conformation" three ways depending on the
#: source: an empty PDB altLoc column is ``" "``, an mmCIF one ``"."``, a
#: missing value ``"?"``.  The field view normalises all of them, so a rule
#: saying "no altloc" has one spelling to write (the empty cell).
NULL_ALT = ("", " ", ".", "?")

#: 这两个字段的**空模式**不是"任意", 而是"没有" (无插入码 / 无 altloc) —— 否则
#: `,A,45,CA` 会悄悄连 `45A`、`A`-altloc 的原子一起命中。"任意"要写 `*`。
NOTHING_MEANS_EMPTY = ("ins_code", "altloc")


FIELD_ALIASES = {
    "ins_code": ("ins_code", "icode", "ins", "insertion_code", "insertion"),
    "chain": ("chain", "chain_id", "ch", "auth_asym_id"),
    "resi": ("resi", "res_id", "resid", "residue", "residue_number", "resseq", "auth_seq_id"),
    "name": ("name", "atom", "atom_name", "atom_id", "auth_atom_id"),
    "altloc": ("altloc", "altloc_id", "alt_id", "label_alt_id", "alt"),
}

INT_RANGE = re.compile(r"^(-?\d+)-(-?\d+)$")
SUFFIX_RANGE = re.compile(r"^(.*?)(\d+)-(.*?)(\d+)$")
PLAIN_RANGE = re.compile(r"^(.+)-(.+)$")
IS_INT = re.compile(r"^-?\d+$")


# ---------------------------------------------------------------- 模式匹配器

def _match_one(item: str, value: str, field: str) -> bool:
    """列表里的一项是否命中 value。"""
    if not item:
        return False
    if item.startswith("re:"):
        return re.fullmatch(item[3:], value) is not None
    if field in NOTHING_MEANS_EMPTY and item in ("-", "none"):
        return value == ""                                  # 显式写"无插入码"/"无 altloc"
    if "*" in item or "?" in item:
        return fnmatchcase(value, item)
    if field == "name":                                     # H1-H20 / CA-CB
        m = SUFFIX_RANGE.match(item)
        if m and m.group(1) == m.group(3):
            lo, hi = int(m.group(2)), int(m.group(4))
            vm = re.fullmatch(re.escape(m.group(1)) + r"(\d+)", value)
            return bool(vm) and lo <= int(vm.group(1)) <= hi
    m = INT_RANGE.match(item)                               # 1-10 / -3-5 (A-D 不匹配)
    if m and (field == "resi" or IS_INT.match(value)):
        return IS_INT.match(value) is not None and int(m.group(1)) <= int(value) <= int(m.group(2))
    m = PLAIN_RANGE.match(item)
    if m:                                                   # 单字符按字符序, 其余按字符串序
        lo, hi = m.group(1), m.group(2)
        if field in ("chain", "ins_code", "altloc") and len(lo) == len(hi) == 1:
            return len(value) == 1 and ord(lo) <= ord(value) <= ord(hi)
        if field == "name":
            return lo <= value <= hi
    return item == value


def match_field(pattern: str, value: str, field: str) -> bool:
    """字段模式 pattern 是否命中 value。field ∈ FIELDS。"""
    if not pattern:
        # 空 = 任意; ins_code / altloc 例外 (见 NOTHING_MEANS_EMPTY)
        return value == "" if field in NOTHING_MEANS_EMPTY else True
    if pattern.startswith("re:"):
        return re.fullmatch(pattern[3:], value) is not None
    return any(_match_one(i.strip(), value, field) for i in pattern.split(","))


def _is_literal(pattern: str) -> bool:
    """是不是一个可拿去查表/报错的字面值 (无通配/正则/列表/范围)。"""
    return bool(pattern) and not pattern.startswith("re:") and not any(c in pattern for c in "*?[,")


# ---------------------------------------------------------------- GUI 的「模式」↔ csv 字符串

MODES = ("任意", "字面值", "通配符", "列表", "范围", "正则")

SYNTAX_CHARS = ",*?[]"          # 这几个字符在模式语法里有含义, 字面值里出现就得转义


def _needs_escape(text: str) -> bool:
    """字面值会不会被当成语法 ("A,B" 是列表, "CA-CB" 是范围, "re:x" 是正则)。"""
    return (any(c in text for c in SYNTAX_CHARS) or text.startswith("re:")
            or bool(PLAIN_RANGE.match(text)))


def encode(field: str, mode: str, text: str) -> str:
    """(模式, 框里的输入) → 写进规则表的模式字符串。"""
    text = text.strip()
    if mode == "任意":
        return "*" if field in NOTHING_MEANS_EMPTY else ""    # 空串在这些字段另有含义
    if not text:
        return ""
    if mode == "正则":
        return "re:" + text
    if mode == "字面值" and _needs_escape(text):
        return "re:" + re.escape(text)                       # 免得 "A,B" 被当成列表
    return text


def decode(field: str, pattern: str):
    """规则表里的模式字符串 → (模式, 框里显示的输入)。手写的 csv 也认。"""
    p = pattern.strip()
    if not p:
        return ("字面值", "") if field in NOTHING_MEANS_EMPTY else ("任意", "")
    if p.startswith("re:"):
        return "正则", p[3:]
    if "," in p:
        return "列表", p
    if INT_RANGE.match(p) or PLAIN_RANGE.match(p):
        return "范围", p
    if "*" in p or "?" in p:
        return "通配符", p
    return "字面值", p


# ---------------------------------------------------------------- csv 读写

def find_columns(header) -> dict:
    """表头 → FIELDS 各字段的列号; 缺列直接报错, 不猜。"""
    idx = {}
    for i, h in enumerate(header):
        key = str(h).strip().lower().replace(" ", "_").replace("-", "_")
        for field, aliases in FIELD_ALIASES.items():
            if field not in idx and key in aliases:
                idx[field] = i
    missing = [f for f in FIELDS if f not in idx]
    if missing:
        raise ValueError(f"csv 缺少列 {missing}; 现有表头: {list(header)}")
    return idx


def _read_csv(path) -> tuple[list[list[str]], list[str]]:
    """→ (数据行, 表头)。"""
    with open(path, newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.reader(fh))
    if not rows:
        return [], list(FIELDS)
    return rows[1:], [h.strip() for h in rows[0]]


def _write_csv(path, header, rows) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


# ---------------------------------------------------------------- 匹配目标

class _TargetView:
    """选择所需的字段 + 展示用残基名, 全按 PDB/auth 口径。

    同时吃得下 ``AtomArray`` 与 ``InternalCoord``: 两者都暴露
    ``ins_code`` / ``chain_id`` / ``res_id`` / ``atom_name`` / ``res_name``。

    ``altloc`` 的值来自目标的 ``altloc_id``: ``AtomArray`` 上那是 biotite 的标注 —— 只有**按
    ``altloc="all"`` 读进来**的数组才有 (实测 biotite 1.6: 默认 ``altloc="first"`` 不带这个
    category, 见 :class:`biorazer.structure.io.protein.StructureFile_AtomArray`)。
    ``InternalCoord`` 上没有这个标注 (记录不带 altloc, 见
    :class:`~biorazer.structure.bridge.AtomArray_InternalCoord`), 于是一切目标都不带标注时
    ``has_altloc=False``, 所有原子的 altloc 记成 ``""`` (= 无 altloc)。"没有"的三种 sentinel
    (PDB 空列 ``" "`` / CIF ``"."`` / ``"?"``) 一律归一成 ``""`` (见 :data:`NULL_ALT`)。
    """

    def __init__(self, target):
        self.vals = {
            "ins_code": [str(x) for x in target.ins_code],
            "chain": [str(x) for x in target.chain_id],
            "resi": [str(x) for x in target.res_id],
            "name": [str(x) for x in target.atom_name],
        }
        self.n = len(self.vals["chain"])
        raw_alt = getattr(target, "altloc_id", None)
        self.has_altloc = raw_alt is not None
        if raw_alt is None:
            raw_alt = [""] * self.n
        self.vals["altloc"] = ["" if str(x) in NULL_ALT else str(x) for x in raw_alt]
        self.res = [str(x) for x in target.res_name]
        self.chains = sorted(set(self.vals["chain"]))
        self.ins_codes = sorted(x for x in set(self.vals["ins_code"]) if x)
        self.altlocs = sorted(x for x in set(self.vals["altloc"]) if x)
        keys = Counter(self.key(a) for a in range(self.n))
        self.dup_keys = {k for k, c in keys.items() if c > 1}

    def key(self, a: int):
        """(ins_code, chain, resi, name, altloc) —— 与 FIELDS 同序, 方便和 csv 行直接比。"""
        return tuple(self.vals[f][a] for f in FIELDS)

    def atom_str(self, a: int) -> str:
        """A/45B/CA —— 插入码直接跟在 resi 后面; 目标带 altloc 时再跟 ``(A)``。"""
        alt = self.vals["altloc"][a]
        return (f"{self.vals['chain'][a]}/{self.vals['resi'][a]}"
                f"{self.vals['ins_code'][a]}/{self.vals['name'][a]}"
                + (f"({alt})" if alt else ""))

    def summary(self) -> str:
        ins = ",".join(self.ins_codes) or "无"
        alt = (",".join(self.altlocs) or "无") if self.has_altloc else "无标注"
        return (f"{self.n} 原子, 链 {','.join(self.chains)}, 插入码 {ins}, "
                f"altloc {alt}, 重复五元组 {len(self.dup_keys)}")


def match_rule(patterns: dict, view: _TargetView) -> tuple[list[int], list[str]]:
    """一行规则 → (命中的原子下标列表, 警告列表)。"""
    warns = []
    if not view.has_altloc and patterns.get("altloc", "").strip() not in ("", "*"):
        warns.append("目标没有 altloc 标注 (InternalCoord 不带; 手工拼的 AtomArray 也不带; "
                     "读文件请用 StructureFile_AtomArray), altloc 模式无法判定")
    try:
        # ponytail: 每行 5×N 次匹配; 结构 >10^6 原子时再预编译谓词
        hits = [a for a in range(view.n)
                if all(match_field(patterns[f], view.vals[f][a], f) for f in FIELDS)]
    except re.error as exc:
        return [], warns + [f"正则错误: {exc}"]
    dup = sum(1 for a in hits if view.key(a) in view.dup_keys)
    if dup:
        warns.append(f"{dup} 个命中原子的五元组在目标里不唯一 (altloc 等)")
    if not hits:
        warns.append("0 命中")
        if _is_literal(patterns["chain"]) and patterns["chain"] not in view.chains:
            warns.append(f"chain {patterns['chain']!r} 不在目标里 "
                         f"(只有 {','.join(view.chains)})")
    return hits, warns


def _unreachable(ic: InternalCoord) -> list[int]:
    """从 anchor 出发在 bond 图 (``dihedra`` 的生长规则) 上够不到的原子下标。

    判据与 :meth:`~biorazer.structure.objects.InternalCoord.to_coords` 的同一个: 一个原子
    可放 ⟺ 它在 ``anchor`` 里, 或它有一条 ``dihedra`` 的四元组, 其三个父原子都可放。
    这里只走图, 不算坐标。
    """
    placed = set(ic.anchor)
    changed = True
    while changed:
        changed = False
        for i, j, k, l in ic.dihedra:
            if l not in placed and i in placed and j in placed and k in placed:
                placed.add(l)
                changed = True
    return sorted(set(range(len(ic))) - placed)


# ---------------------------------------------------------------- 选择器

@dataclass(slots=True, repr=False)
class AtomArraySelector:
    """一组选择式 (规则表), 可以匹配 ``AtomArray`` / ``InternalCoord``, 也可以读写 csv。

    Fields
    ------
    rules : list[list[str]]
        每行一条规则, 每格一条**模式**; 行内单元格与 ``header`` 对齐。字段列由
        :func:`find_columns` 按表头名 (别名见 ``FIELD_ALIASES``) 定位, 其余列是随行的
        批注, 展开时照抄到它展开出的每个原子行。
    header : list[str]
        列名, 默认就是 ``FIELDS``。

    ``repr=False`` 保留下面的短 repr (生成的 repr 会把每行规则都倒出来)。相等按字段值
    (两组规则 + 表头) 比较。
    """

    rules: list[list[str]] = field(default_factory=list)
    header: list[str] = field(default_factory=lambda: list(FIELDS))

    def __post_init__(self):
        self.header = [str(h).strip() for h in (self.header or FIELDS)]
        width = len(self.header)
        self.rules = [[str(c) for c in (row or [])] for row in (self.rules or [])]
        longer = [i for i, row in enumerate(self.rules) if len(row) > width]
        if longer:
            # 比表头还长的行 = csv 里那个多出来的逗号没加引号 (或列数真对不上)。
            # 补空好办, 多的单元格没有列名可挂, 只能报错 —— 静默丢掉就是改用户数据。
            raise ValueError(
                f"row(s) {longer} have more cells than the header ({width}): "
                f"{self.rules[longer[0]]!r}; 表头 {self.header}")
        # 短行补空, 后面按列号取值才不用每次判越界
        self.rules = [row + [""] * (width - len(row)) for row in self.rules]
        self.columns          # 表头缺列在这里就报错, 不等到匹配时

    def __repr__(self):
        return f"AtomArraySelector({len(self.rules)} rules, header={self.header})"

    @property
    def columns(self) -> dict:
        """字段名 → 列号 (缺列抛 ``ValueError``)。"""
        return find_columns(self.header)

    @property
    def patterns(self) -> list[dict]:
        """每行规则取成 ``{字段: 模式}`` (字段名已解析, 值已 strip)。"""
        cols = self.columns
        return [{f: row[cols[f]].strip() for f in FIELDS} for row in self.rules]

    # ---- csv ------------------------------------------------------------

    @classmethod
    def from_csv(cls, path, mode: str = "rule") -> "AtomArraySelector":
        """从 csv 构建选择器。

        Parameters
        ----------
        path : str or Path
            规则表 / 选择表路径。
        mode : {"rule", "selection"}
            ``"rule"`` 每格原样当**模式**读 (规则表); ``"selection"`` 每格是**具体值**,
            先经 :func:`encode` 按字面值转义再当模式读 —— 否则选择表里那些具体值
            (``A,B`` 是列表、``CA-CB`` 是范围、``re:x`` 是正则) 会被当语法误解。
            这样选择表也能当规则表打开 (见模块 docstring)。

        Returns
        -------
        AtomArraySelector
        """
        if mode not in ("rule", "selection"):
            raise ValueError(f"mode must be 'rule' or 'selection', got {mode!r}")
        rows, header = _read_csv(path)
        if mode == "selection":
            field_of = {i: f for f, i in find_columns(header).items()}
            rows = [[encode(field_of[i], "字面值", cell) if i in field_of else cell
                     for i, cell in enumerate(row)] for row in rows]
        return cls(rules=rows, header=header)

    def to_csv(self, path, mode: str = "rule", structure=None, dedupe: bool = True) -> None:
        """导出选择器。

        Parameters
        ----------
        path : str or Path
            输出路径。
        mode : {"rule", "selection"}
            ``"rule"`` 原样写回规则表 (与 :meth:`from_csv` 的 ``"rule"`` 往返一致);
            ``"selection"`` 按 ``structure`` 把每行规则展开成**每行一个原子**的选择表,
            字段列写具体值, 其余列照抄。``structure`` 必需 —— 它是展开的唯一依据。
        structure : AtomArray or InternalCoord or None
            展开规则表用的目标 (只有 ``mode="selection"`` 用得到)。
        dedupe : bool
            展开时去掉重复原子 (同一个五元组只留第一次命中的那条规则)。
        """
        if mode == "rule":
            _write_csv(path, self.header, self.rules)
            return
        if mode != "selection":
            raise ValueError(f"mode must be 'rule' or 'selection', got {mode!r}")
        if structure is None:
            raise ValueError("mode='selection' needs a structure: the selection table "
                             "is the rule table expanded against one")
        view, _hits, _keep, _warns, _dup = self._match(structure)
        _write_csv(path, self.header,
                   self.expanded_rows(view, self.pairs(structure, dedupe=dedupe)))

    # ---- 匹配 -----------------------------------------------------------

    @staticmethod
    def _one_model(target):
        """目标必须是**一个 model**: ``AtomArrayStack`` 的布尔下标切的是 model 轴,
        拿它当 AtomArray 会把选择悄悄变成"选 model"。

        Raises
        ------
        TypeError
            ``target`` 是 ``AtomArrayStack`` (先自己 ``stack[0]``)。
        """
        if isinstance(target, AtomArrayStack):
            raise TypeError(
                "AtomArraySelector works on a single model: this is an AtomArrayStack. "
                "Index the model you mean first (e.g. stack[0]) -- boolean indexing a "
                "stack selects models, not atoms.")
        return target

    def _match(self, target):
        """一行一行匹配。

        Returns
        -------
        tuple
            ``(view, hits, keep, warns, dup)`` -- ``hits[r]`` 是规则 ``r`` 命中的原子下标
            (按 target 顺序), ``keep[r]`` 只留其中五元组第一次出现的那几个, ``warns`` 是
            ``{行号: [警告]}``, ``dup`` 是命中里五元组重复的总个数。
        """
        view = _TargetView(self._one_model(target))
        hits, keep, warns, dup, seen = [], [], {}, 0, set()
        for r, patterns in enumerate(self.patterns):
            row_hits, row_warns = match_rule(patterns, view)
            if row_warns:
                warns[r] = row_warns
            first = []
            for a in row_hits:
                key = view.key(a)
                if key in seen:
                    dup += 1
                else:
                    seen.add(key)
                    first.append(a)
            hits.append(row_hits)
            keep.append(first)
        return view, hits, keep, warns, dup

    def resolve(self, target):
        """匹配每一行规则。

        Returns
        -------
        hits : list[list[int]]
            每行规则命中的原子下标 (按 target 顺序)。
        warns : dict[int, list[str]]
            ``{行号 (0 起): [警告]}`` —— 正则错误 / 0 命中 / 命中的五元组不唯一。
        dup : int
            命中里五元组重复的总个数 (按规则顺序算, 先出现的不计)。
        """
        _view, hits, _keep, warns, dup = self._match(target)
        return hits, warns, dup

    def pairs(self, target, dedupe: bool = True) -> list[tuple[int, int]]:
        """展平成 ``[(规则行号, 原子下标)]``, 顺序 = 规则顺序, 行内 = target 里的原子顺序。

        ``dedupe=True`` (默认) 时同一个五元组只保留第一次命中它的那条规则。
        """
        _view, hits, keep, _warns, _dup = self._match(target)
        chosen = keep if dedupe else hits       # 重复原子记在第一次命中它的规则上
        return [(r, a) for r, row in enumerate(chosen) for a in row]

    def indices(self, target, dedupe: bool = True) -> np.ndarray:
        """展平后的原子下标 (1D int 数组), 可用于直接索引 target。"""
        return np.array([a for _r, a in self.pairs(target, dedupe=dedupe)], dtype=int)

    def mask(self, target) -> np.ndarray:
        """与 target 对齐的布尔 mask (匹配顺序无关, 所以不用去重)。"""
        mask = np.zeros(len(target), dtype=bool)
        mask[self.indices(target, dedupe=False)] = True
        return mask

    def expanded_rows(self, view: _TargetView, pairs) -> list[list[str]]:
        """展平成选择表: 每行一个原子, 字段列写具体值, 规则行的其余列照抄。"""
        cols = self.columns
        out = []
        for r, a in pairs:
            cells = list(self.rules[r])
            for f, i in cols.items():
                cells[i] = view.vals[f][a]
            out.append(cells)
        return out

    # ---- 应用 -----------------------------------------------------------

    def apply(self, target):
        """返回选择后的目标。

        * ``AtomArray`` (或任何支持布尔下标的东西) → ``target[self.mask(target)]``;
        * ``InternalCoord`` → ``target[self.mask(target)]``, 且**结果里每个原子都必须还在
          bond 图里** —— ``InternalCoord`` 的子集会把跨出选择的 ``dihedra`` / ``bond_angles``
          / ``bond_distances`` 条目丢掉, 于是断口处的原子失去生长依据, 得到的是一个
          "长不出来"的图。这里在返回前按生长规则 (anchor + dihedra) 走一遍图, 有够不到的
          原子就抛 ``ValueError`` 并把它们按 ``atom_repr`` 列出来 —— 即**选择必须保住连通性**
          (要么选整个连通片段, 要么自己给 ``InternalCoord.anchor`` 补帧)。

        Parameters
        ----------
        target : AtomArray or InternalCoord
            匹配与切片的对象 (两者的注解字段同名同义)。

        Returns
        -------
        AtomArray or InternalCoord
            与 ``target`` 同类型的选择结果。
        """
        if not isinstance(target, InternalCoord):
            return target[self.mask(target)]
        sub = target[self.mask(target)]
        bare = _unreachable(sub)
        if bare:
            raise ValueError(
                "selection leaves atoms without a bond connection: "
                + ", ".join(sub.atom_repr(i) for i in bare)
                + " (their grow frame was cut off; select whole connected "
                  "fragments or give InternalCoord.anchor a frame)")
        return sub

    # ---- 编辑器 ---------------------------------------------------------

    def run_editor(self, structure_file=None, rule_csv=None, selection_csv=None,
                   dedupe: bool = True, on_ready=None):
        """打开原子选择表编辑器 (tkinter), 三个路径参数都可在 GUI 里选。

        Parameters
        ----------
        structure_file : str or Path or None
            PDB/CIF 结构 —— 规则展开的唯一依据; 留空则在 GUI 里 "打开结构…"。
        rule_csv : str or Path or None
            起始规则表; 留空则用本选择器已有的 ``rules``。
        selection_csv : str or Path or None
            导出选择表的默认路径 (留空: ``selection.csv``)。
        dedupe : bool
            展开时是否去重 (GUI 里是可勾选的初始值)。
        on_ready : callable or None
            自检钩子: 传了就在窗口起来后把 GUI 的命令/状态字典交给它, 然后自动关窗 ——
            供测试驱动全流程, 正常使用不传。

        Returns
        -------
        AtomArraySelector
            本选择器 (GUI 里保存/导出会就地更新它的 ``rules`` / ``header``)。
        """
        from biorazer.structure.objects._selector_gui import run_editor
        return run_editor(self, structure_file=structure_file, rule_csv=rule_csv,
                          selection_csv=selection_csv, dedupe=dedupe, on_ready=on_ready)


# ---------------------------------------------------------------- 不开窗的用法

def print_report(selector: AtomArraySelector, structure, dedupe: bool = True,
                 out=None) -> int:
    """打印每行规则的命中数/问题行; 给了 ``out`` 就顺手把选择表写过去。

    Returns
    -------
    int
        0 = 每行都有命中且无警告, 1 = 有问题行 (可直接当退出码)。
    """
    view = _TargetView(structure)
    hits, warns, dup = selector.resolve(structure)
    print(f"结构: {view.summary()}")
    print(f"规则表: {len(selector.rules)} 行, 列 {selector.header} -> {selector.columns}")
    for r, row in enumerate(selector.rules):
        shown = " ".join(view.atom_str(a) for a in hits[r][:12])
        print(f"行 {r + 1} {','.join(row):<36} → {len(hits[r]):>5} 原子  {shown}"
              + (f" … (+{len(hits[r]) - 12})" if len(hits[r]) > 12 else "")
              + "".join(f"\n      ⚠ {w}" for w in warns.get(r, [])))
    keep = selector.pairs(structure, dedupe=dedupe)
    print(f"--- 规则 {len(selector.rules)} 行 → 展平 {len(selector.indices(structure, dedupe=False))} "
          f"个原子, 重复 {dup}, {len(warns)} 行有问题; "
          f"选择表 {len(keep)} 个原子" + (" (已去重)" if dedupe else " (未去重)"))
    if out is not None:
        selector.to_csv(out, mode="selection", structure=structure, dedupe=dedupe)
        print(f"选择表 {out}: {len(keep)} 行 (每行一个原子)")
    return 1 if warns else 0


def _add_selector_parser(sub):
    """把 ``biorazer select`` 子命令挂到 argparse subparsers 上 (见 ``objects/cli.py``)。"""
    p = sub.add_parser(
        "select",
        help="规则表 → 每行一个原子的选择表 (不开窗打印/导出, 或打开编辑器)",
        description=(
            "用规则表把结构展开成**选择表** (严格每行一个原子), 或在 tkinter 编辑器里编辑规则表:\n"
            "规则表每行一条选择式, 5 列 ins_code/chain/resi/name/altloc, 每格是一条模式\n"
            "(空=任意; 字面值/通配符/列表/范围/正则; ins_code 与 altloc 的空=没有)。\n"
            "默认开窗; --print / --export 不开窗 (headless 重建)。"
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-s", "--structure", metavar="PDB|CIF",
                   help="结构文件 (按后缀自动识别); 规则展开的唯一依据, 留空则开窗后选")
    p.add_argument("--rules", metavar="CSV", default="rules.csv",
                   help="规则表: 每行一条选择式, 与编辑器的行一一对应")
    p.add_argument("-o", "--out", metavar="CSV", default="selection.csv",
                   help="选择表: 展平后的产物, 严格每行一个原子")
    p.add_argument("--no-dedupe", dest="dedupe", action="store_false",
                   help="选择表保留重复原子 (默认去重: 同五元组只写第一次命中的那条规则)")
    p.add_argument("--print", dest="do_print", action="store_true",
                   help="不开窗: 打印每行规则的命中数与问题行")
    p.add_argument("--export", action="store_true",
                   help="配合 --print: 顺手把选择表写到 -o")
    p.set_defaults(func=_run_selector)
    return p


def _run_selector(args) -> None:
    """``biorazer select`` 的执行体。"""
    from biorazer.structure.io import StructureFile_AtomArray

    selector = (AtomArraySelector.from_csv(args.rules) if Path(args.rules).exists()
                else AtomArraySelector())
    if args.do_print or args.export:
        if not args.structure:
            raise SystemExit("--print/--export 需要 -s/--structure")
        array = StructureFile_AtomArray(input_io=args.structure).read()
        print_report(selector, array, dedupe=args.dedupe,
                     out=args.out if args.export else None)
        return
    selector.run_editor(structure_file=args.structure, rule_csv=args.rules,
                        selection_csv=args.out, dedupe=args.dedupe)
