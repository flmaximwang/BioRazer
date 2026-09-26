"""Converters for an ``AtomArraySelection``: its two csv tables, and the per-atom
selection it makes on an ``AtomArray``.

Every class here is a biorazer ``A_B`` :class:`~biorazer.io.Converter` -- ``A`` is the
source, ``B`` the target.

**The two csv tables** (one selection syntax, two uses; see
:mod:`biorazer.structure.objects.selector`):

* :class:`RuleCsv_AtomArraySelection` / :class:`AtomArraySelection_RuleCsv` -- the
  **rule table**: one selection expression per row, each cell a pattern.  The
  round trip is lossless, ``0``-hit rows included: the rule table is the source of
  truth and only the expansion drops rows.
* :class:`SelectionCsv_AtomArraySelection` / :class:`AtomArraySelection_SelectionCsv`
  -- the **selection table**: exactly one atom per row, concrete values.  The
  reader escapes every field cell as a literal, so those concrete values (``A,B``,
  ``CA-CB``, ``re:x``) are not misread as syntax and a selection table can be
  opened as a rule table.  The writer needs the structure the table is expanded
  against -- it is the only thing the expansion can depend on.

These four are file-backed (``read()`` / ``write()``) because a csv table is a
file: the rule/selection table has no in-memory class of its own, unlike the
objects :mod:`biorazer.structure.bridge` normally bridges.  With that one
exception the module follows the package rule -- in-memory conversion, no files.

**The PyMOL selection text** -- :class:`PyMOLSelection_AtomArraySelection` /
:class:`AtomArraySelection_PyMOLSelection` move the same selection to and from a
chunk of text you can paste straight into PyMOL: one ``select <name>, ...`` line
whose terms are ``/model//chain/resi/name`` macros (the insertion code rides on ``resi``; the
macro has no altloc slot, so the text does not split conformers), one macro per matched atom,
joined with ``or``.  Nothing
attempts to translate a *pattern* (``*`` / ``1-10`` / ``re:``) into PyMOL's own
selection algebra -- the text is computer-generated, so a concrete macro per atom
is both easier and exact, and the whole "PyMOL cannot express this pattern"
problem never comes up.

The remaining PyMOL classes are compositions over the same middle layer, not a
second matching path:

* :class:`PyMOLSelection_SelectionCsv` -- text -> rules -> selection table (it
  still needs the structure the table is expanded against);
* :class:`SelectionCsv_PyMOLSelection` -- selection table -> rules -> text (no
  structure needed: every row is already a concrete atom);
* :class:`Mask_PyMOLSelection` / :class:`Indices_PyMOLSelection` -- a mask or a
  flat index array -> rules (one literal rule per selected atom) -> text; the
  values still come from the structure, so both take it as an argument.

**The per-atom selection** -- :class:`AtomArraySelection_AtomArrayMask` and
:class:`AtomArraySelection_AtomArrayIndices` implement the parent's
:meth:`~biorazer.io.Converter.convert` (one transform, no ``read()`` /
``write()``): the selection is held in ``input_io`` and ``convert()`` returns the
selection over an array passed per call.  The array is a parameter because only
an ``AtomArray`` says which atoms exist -- the shape the ``quads`` / ``anchor``
parameters of :class:`~biorazer.structure.bridge.atom_array.AtomArray_InternalCoord`
already use.

The two targets differ in ordering, which is why both exist:

* a **mask** is aligned with the array, so the rule order is invisible in it, and
* **indices** are flat (rule order first, then atom order within a rule), an order
  a mask cannot express.  ``dedupe=False`` keeps atoms matched by more than one
  rule once per match.
"""

from __future__ import annotations

import csv
import re

import numpy as np

from biorazer.io import Converter
from biorazer.structure.objects import AtomArray, AtomArraySelection
from biorazer.structure.objects.selector import (
    FIELDS,
    _TargetView,
    decode,
    encode,
    find_columns,
)


# ---------------------------------------------------------------- csv 文件

def _read_csv(path) -> tuple[list[list[str]], list[str]]:
    """``path`` → (数据行, 表头); 空文件当"没有数据行, 表头 = ``FIELDS``"。"""
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


class RuleCsv_AtomArraySelection(Converter):
    """Reads a **rule table** into an :class:`AtomArraySelection`.

    Every cell is a pattern, read as written (``*`` wildcard, ``,`` list,
    ``1-10`` range, ``re:`` regex, else a literal).  This is the file the
    tkinter editor edits, so this direction plus
    :class:`AtomArraySelection_RuleCsv` is the lossless round trip.

    Parameters
    ----------
    input_io : str or Path
        Rule table csv; the field columns are found by header name
        (:data:`~biorazer.structure.objects.selector.FIELD_ALIASES`), extra
        columns ride along as per-row notes.
    """

    def read(self) -> AtomArraySelection:
        """Read ``self.input_io``.

        Returns
        -------
        AtomArraySelection
            One rule per data row; missing field columns raise ``ValueError``.
        """
        rows, header = _read_csv(self.input_io)
        return AtomArraySelection(rules=rows, header=header)


class SelectionCsv_AtomArraySelection(Converter):
    """Reads a **selection table** (exactly one atom per row) into a selection.

    Each field cell is taken as a **literal value** and escaped through
    :func:`~biorazer.structure.objects.selector.encode` before it becomes a
    pattern -- otherwise the concrete values a selection table holds would be
    misread as syntax (``A,B`` a list, ``CA-CB`` a range, ``re:x`` a regex).
    That is what lets a selection table be reopened as a rule table: every
    concrete value then matches only itself.

    Parameters
    ----------
    input_io : str or Path
        Selection table csv.
    """

    def read(self) -> AtomArraySelection:
        """Read ``self.input_io``.

        Returns
        -------
        AtomArraySelection
            One single-atom rule per data row.
        """
        rows, header = _read_csv(self.input_io)
        field_of = {i: f for f, i in find_columns(header).items()}
        return AtomArraySelection(
            rules=[[encode(field_of[i], "字面值", cell) if i in field_of else cell
                    for i, cell in enumerate(row)] for row in rows],
            header=header)


class AtomArraySelection_RuleCsv(Converter):
    """Writes a selection back to its **rule table**, row for row.

    Parameters
    ----------
    output_io : str or Path
        Where the rule table goes.
    """

    def write(self, tmp: AtomArraySelection) -> None:
        """Write ``tmp`` to ``self.output_io``.

        Parameters
        ----------
        tmp : AtomArraySelection
            The selection whose ``header`` / ``rules`` are written verbatim --
            including the rules that match nothing.
        """
        _write_csv(self.output_io, tmp.header, tmp.rules)


class AtomArraySelection_SelectionCsv(Converter):
    """Expands a selection against a structure, one atom per row.

    Parameters
    ----------
    output_io : str or Path
        Where the selection table goes.
    """

    def write(self, tmp: AtomArraySelection, structure, dedupe: bool = True) -> None:
        """Write the expanded ``tmp`` to ``self.output_io``.

        Parameters
        ----------
        tmp : AtomArraySelection
            The rules to expand; its ``header`` is written as the table header
            and its extra columns are copied onto every row the rule expands to.
        structure : AtomArray or InternalCoord
            The target the expansion is made against -- the only thing the
            expansion can depend on, so it is an argument here rather than
            something the selection could hold.
        dedupe : bool
            Keep an atom matched by more than one rule only at its first match
            (the default), or once per matching rule.

        Raises
        ------
        ValueError
            ``structure`` is ``None``.
        """
        if structure is None:
            raise ValueError("a selection table is the rule table expanded against a "
                             "structure: pass one, there is nothing to expand otherwise")
        view, hits, keep, _warns, _dup = tmp._match(structure)
        chosen = keep if dedupe else hits
        _write_csv(self.output_io, tmp.header,
                   tmp.expanded_rows(view, [(r, a) for r, row in enumerate(chosen)
                                            for a in row]))


# ---------------------------------------------------------------- PyMOL 选择式文本

#: PyMOL 宏 ``/model//chain/resi/name`` 去掉前导 ``/`` 后按 ``/`` 切出的段数。
MACRO_PARTS = 5


def _pymol_macro(model: str, vals: dict) -> str:
    """一个原子的字段值 → PyMOL 宏 ``/model//A/45A/CA``。

    * 插入码没有自己的段, 直接跟在 ``resi`` 后面: 实测 ``/m//A/45A/CA`` 只命中带插入码的
      那个原子、``/m//A/45/CA`` 只命中不带的 —— 与本库 ``ins_code`` 的"空 = 没有"一致;
    * **宏不带 altloc** (接第 6 段直接 ``too many slashes``), 于是这段文本**不区分构象**:
      同名同号的多个构象在 PyMOL 里一起命中 (实测 ``/m//A/2/CB`` 命中 A、B 两份)。这是用户
      要的语义 —— 选择要么按残基/原子看整体, 要么由规则表里那格 altloc 去约束, 不在文本这层
      挑构象;
    * ``model`` 与 ``chain`` / ``resi`` / ``name`` 一样是必需值: 实测 ``//A/1/CA`` (空 model)
      命中 **0** 个原子, 给个空值出去就是一段选不中东西的文本;
    * ``ins_code`` 只能是单字符或空 —— 别让 ``none`` 这种模式哨兵值混进宏里。
    """
    if not model:
        raise ValueError(
            "PyMOL 宏的第一个 ``/`` 段就是对象名, 空着 (如 ``//A/1/CA``) 实测命中 0 个原子: "
            "传 model (cmd.load 时给的名字, 缺省是文件名)。")
    for field in ("chain", "resi", "name"):
        if not vals[field]:
            raise ValueError(f"原子 {vals} 的 {field} 是空的: PyMOL 宏省掉一格就等于把选择"
                             "放宽, 不写这种文本。")
    if len(vals["ins_code"]) > 1:
        raise ValueError(f"原子 {vals} 的 ins_code = {vals['ins_code']!r}: 只能是单字符或空")
    return f"/{model}//{vals['chain']}/{vals['resi']}{vals['ins_code']}/{vals['name']}"


def _pymol_macro_to_cells(term: str, lineno: int) -> list[str]:
    """PyMOL 宏 → 一行规则 (每格按字面值转义)。

    ``model`` / ``segi`` 在规则表里没有列: ``model`` 是 PyMOL 侧的对象名 (原子身份里没有它),
    ``segi`` 本库不认 —— 非空 ``segi`` 报错, 不静默丢掉。``altloc`` 宏里没有这一格, 而 PyMOL
    不带它就是**任意构象** —— 所以记成 ``*`` (任意), 不是"没有 altloc" (那会静默把选择收窄)。
    """
    macro, and_, suffix = term.strip().partition(" and ")
    if and_:
        raise ValueError(f"选择式第 {lineno} 行: {term.strip()!r} 里宏后面不该有别的东西 —— "
                         "宏是 `/model//chain/resi/name` (插入码跟在 resi 后), 不接关键字")
    parts = macro.strip().lstrip("/").split("/")
    if len(parts) != MACRO_PARTS:
        raise ValueError(
            f"选择式第 {lineno} 行: {term.strip()!r} 不是 PyMOL 宏 "
            f"(/model//chain/resi/name 共 {MACRO_PARTS} 段), 拿到 {len(parts)} 段; "
            "本转换器只读自己写出的宏形式 (规则式选择式走 RuleCsv_AtomArraySelection)")
    _model, segi, chain, resi_name, name = parts
    if segi:
        raise ValueError(f"选择式第 {lineno} 行: {term.strip()!r} 带 segi {segi!r} —— "
                         "规则表没有 segment 这一列, 不猜")
    m = re.fullmatch(r"(-?\d+)([A-Za-z]?)", resi_name)
    if not m:
        raise ValueError(f"选择式第 {lineno} 行: {term.strip()!r} 的 resi {resi_name!r} "
                         "不是具体残基号 (可带插入码, 如 45 或 45A)")
    # 前四格是具体值 (过 encode 转义); altloc 那一格是**模式** ``*`` = 任意构象 —— 不能过
    # encode, 它会把这个通配符当字面值转义成 ``re:\*`` (实测: 那就一个原子都命中不了)
    cells = [encode(f, "字面值", v) for f, v in
             zip(FIELDS[:4], (m.group(2), chain, m.group(1), name))]
    return cells + ["*"]


def _parse_pymol_text(text) -> AtomArraySelection:
    """PyMOL 选择式文本 → 规则表: 每个宏一条字面值规则 (纯原子的"或")。

    认 ``select <name>, <宏> or <宏> ...``; 空行、``#`` 注释、行尾 ``;`` (从 pml 里拷出来的
    那样) 都跳过。``select`` 名与 ``model`` 一样是 PyMOL 侧的东西, 不进规则表。
    """
    rows = []
    for lineno, line in enumerate(str(text).splitlines(), 1):
        line = line.partition("#")[0].strip().rstrip(";").strip()
        if not line:
            continue
        head, _comma, rest = line.partition(",")
        if head.strip().lower().startswith("select"):
            line = rest.strip()
            if not line:
                raise ValueError(f"选择式第 {lineno} 行: select 后面没有选择式")
        rows += [_pymol_macro_to_cells(term, lineno) for term in line.split(" or ")]
    return AtomArraySelection(rules=rows, header=list(FIELDS))


class AtomArraySelection_PyMOLSelection(Converter):
    """展开成一段可直接粘进 PyMOL 的选择式文本: 一个原子一个宏, ``or`` 串起来。

    ``convert()`` **返回** 文本 (落盘自己 ``Path(out).write_text(...)``)。刻意不把模式
    (通配/范围/正则) 翻译成 PyMOL 自己的选择代数: 文本是给机器生成、拿去粘贴的, 每个命中原子
    写成一个具体宏既省事又准 —— 顺带也就没有"PyMOL 表达不了这条模式"这回事。

    已知边界: 宏里没有 altloc 这一格, 所以文本**不区分构象** —— 同名同号的 A/B 两份一起命中
    (实测 ``/m//A/2/CB`` n=2)。这是要的语义 (同时选中多个构象); 反过来, 库这边按构象筛出来的
    原子 (规则表那格 altloc 写死 ``A``) 到了文本里也会放开成两份 —— 文本只回答"哪些原子",
    不回答"哪份构象"。
    """

    def convert(self, structure, model: str, name: str = "sel",
                dedupe: bool = True) -> str:
        """展开 ``self.input_io`` 里的选择。

        Parameters
        ----------
        structure : AtomArray or InternalCoord
            展开的依据 —— 只有它说得出哪些原子存在。
        model : str
            PyMOL 里这个结构的对象名 (``cmd.load`` 给的名字, 缺省是文件名)。
        name : str
            展开出的选择名; 一个选择一条命令, 所以只有这一个名字。
        dedupe : bool
            被多条规则命中的原子只写一次 (默认), 还是每条规则各写一次。

        Returns
        -------
        str
            ``select <name>, /model//A/1/CA or /model//A/2/CB ...``; 一个原子都没命中时是
            ``select <name>, none`` (PyMOL 的空选择)。

        Raises
        ------
        ValueError
            ``structure`` 是 ``None``; ``model`` 为空或字段给不出值 (见 :func:`_pymol_macro`)。
        """
        if structure is None:
            raise ValueError("a PyMOL text is the rule table expanded against a structure: "
                             "pass one, there is nothing to write otherwise")
        view, hits, keep, _warns, _dup = self.input_io._match(structure)
        vals = [{f: view.vals[f][a] for f in FIELDS}
                for row in (keep if dedupe else hits) for a in row]
        return (f"select {name}, "
                + (" or ".join(_pymol_macro(model, v) for v in vals) if vals else "none")
                + "\n")


class PyMOLSelection_AtomArraySelection(Converter):
    """读一段 PyMOL 选择式文本 → 一个原子一条规则的 :class:`AtomArraySelection`。

    只认 :class:`AtomArraySelection_PyMOLSelection` 写出的形状 (逗号 / ``or`` 分隔的宏);
    人工写的规则式选择式 (``chain A and resi 1-10``) 不在其中 —— 那套语法由
    :class:`RuleCsv_AtomArraySelection` 负责, 猜一半比不解析更坏。

    每条规则一个原子, 所以 ``model`` / ``select`` 名 / ``segi`` 不进规则; ``altloc`` 宏里没有
    这一格, 而 PyMOL 不带它就是**任意构象** —— 记成 ``*`` (任意), 不是"没有 altloc" (那会静默
    把选择收窄)。
    """

    def convert(self, text) -> AtomArraySelection:
        """解析 ``text`` (多行也行, 每行的 ``select`` 都算进去)。

        Parameters
        ----------
        text : str
            PyMOL 选择式文本。

        Returns
        -------
        AtomArraySelection
            宏顺序 → 规则顺序, 每个宏一条字面值规则。
        """
        return _parse_pymol_text(text)


def _literal_selection(structure, indices) -> AtomArraySelection:
    """一组原子下标 → 每个原子一条字面值规则的规则表 (mask / indices 转换的中间层)。

    具体值一律过 :func:`encode` 当字面值 (值里带 ``,`` / ``*`` 之类不会被读成语法); 下标去重
    并按结构顺序排 —— 那正是 PyMOL 列选择时的顺序。

    Raises
    ------
    ValueError
        下标越界 (mask 是另一个结构的, 或者 mask 比结构长)。
    TypeError
        ``structure`` 是 ``AtomArrayStack`` (同一道门: 布尔下标在它上面切的是 model 轴)。
    """
    view = _TargetView(AtomArraySelection._one_model(structure))
    idx = sorted({int(i) for i in np.ravel(indices)})
    bad = [i for i in idx if not 0 <= i < view.n]
    if bad:
        raise ValueError(f"atom indices out of range for this structure ({view.n} atoms): "
                         f"{bad[:5]} —— mask/indices 是逐原子的, 别拿另一个结构的来")
    return AtomArraySelection(
        rules=[[encode(f, "字面值", view.vals[f][a]) for f in FIELDS] for a in idx],
        header=list(FIELDS))


def _indices_to_pymol(indices, structure, model: str, name: str) -> str:
    """下标 → 规则表 → 文本 (两个类共用的一条链, 没有第二条匹配路径)。"""
    selection = _literal_selection(structure, indices)
    return AtomArraySelection_PyMOLSelection(input_io=selection).convert(structure, model, name)


class Mask_PyMOLSelection(Converter):
    """布尔 mask → PyMOL 文本 (中间层是规则表: 每个 ``True`` 原子一条字面值规则)。

    mask 是 :mod:`biorazer.structure.selection.mask` 一族函数的产物 (与结构逐原子对齐),
    拿它生成文本要结构: 宏里的 ``chain`` / ``resi`` / ``name`` / ``alt`` 只有结构知道。
    """

    def convert(self, structure, model: str, name: str = "sel") -> str:
        """把 ``self.input_io`` 里的 mask 展开成 PyMOL 文本。

        Parameters
        ----------
        structure : AtomArray or InternalCoord
            mask 对齐的那个结构 (也是值的来源)。
        model : str
            PyMOL 里的对象名 (同上游)。
        name : str
            ``select`` 出来的选择名。

        Returns
        -------
        str
            ``select <name>, /model//A/1/CA or ...``; 全 False 时是 ``select <name>, none``。

        Raises
        ------
        ValueError
            mask 长度与结构原子数不一致。
        """
        structure = AtomArraySelection._one_model(structure)
        mask = np.asarray(self.input_io)
        if mask.shape != (len(structure),):
            raise ValueError(f"mask shape {mask.shape} 对不上结构的 {len(structure)} 个原子: "
                             "mask 是逐原子对齐的, 换结构要重新算")
        return _indices_to_pymol(np.flatnonzero(mask), structure, model, name)


class Indices_PyMOLSelection(Converter):
    """原子下标 (1D int 数组) → PyMOL 文本 (中间层同样是规则表)。

    与 :class:`Mask_PyMOLSelection` 只差输入那一端: 展平的下标 (``selection.index.*`` /
    :class:`AtomArraySelection_AtomArrayIndices` 的产物) 照样要结构 —— 下标只说第几个原子。
    重复下标只写一次, 顺序无意义 (PyMOL 按结构序列出选择)。
    """

    def convert(self, structure, model: str, name: str = "sel") -> str:
        """把 ``self.input_io`` 里的下标展开成 PyMOL 文本 (参数同上游)。

        Parameters
        ----------
        structure : AtomArray or InternalCoord
            下标所指的结构 (也是值的来源)。
        model : str
            PyMOL 里的对象名。
        name : str
            ``select`` 出来的选择名。

        Returns
        -------
        str
            ``select <name>, /model//A/1/CA or ...``; 空下标时是 ``select <name>, none``。

        Raises
        ------
        ValueError
            下标越界 (不是这个结构的下标)。
        """
        return _indices_to_pymol(self.input_io, structure, model, name)


class PyMOLSelection_SelectionCsv(Converter):
    """PyMOL 文本 → 选择表 (严格每行一个原子)。

    中间层就是规则表: 文本先解析成 :class:`AtomArraySelection`, 再交给
    :class:`AtomArraySelection_SelectionCsv` 按结构展开 —— 没有第二条匹配路径。
    """

    def write(self, text, structure, dedupe: bool = True) -> None:
        """把 ``text`` 展开成选择表写到 ``self.output_io``。

        Parameters
        ----------
        text : str
            一段 PyMOL 选择式文本 (见 :class:`PyMOLSelection_AtomArraySelection`)。
        structure : AtomArray or InternalCoord
            展开的依据。
        dedupe : bool
            同一个原子被多个宏写到时只留第一次 (默认)。
        """
        selection = PyMOLSelection_AtomArraySelection().convert(text)
        AtomArraySelection_SelectionCsv(output_io=self.output_io).write(
            selection, structure, dedupe=dedupe)


class SelectionCsv_PyMOLSelection(Converter):
    """选择表 → PyMOL 文本 (中间层同样是规则表)。

    选择表每格是**具体值**, 拿它当 :class:`SelectionCsv_AtomArraySelection` 打开 (具体值按
    字面值转义), 再交给 :class:`AtomArraySelection_PyMOLSelection` 逐个原子写宏 —— 不需要
    结构: 一个原子一条规则, 宏就写完了, 没什么可展开的。

    Raises
    ------
    ValueError
        某格不是具体值 (通配/范围/列表/正则混进了选择表): 那种值写进宏就成了别的意思。
    """

    def read(self, model: str, name: str = "sel") -> str:
        """读 ``self.input_io``, 返回粘进 PyMOL 的文本。

        Parameters
        ----------
        model : str
            PyMOL 里的对象名 (同上游)。
        name : str
            ``select`` 出来的选择名。

        Returns
        -------
        str
            ``select <name>, /model//A/1/CA or ...``; 空表是 ``select <name>, none``。
        """
        selection = SelectionCsv_AtomArraySelection(input_io=self.input_io).read()
        cols = selection.columns
        macros = []
        for row in selection.rules:
            vals = {f: row[cols[f]].strip() for f in FIELDS}
            for field in FIELDS:
                if decode(field, vals[field])[0] != "字面值":
                    raise ValueError(
                        f"选择表 {field} 列有一格 {vals[field]!r} 不是具体值 (通配/范围/列表/"
                        "正则): 写进 PyMOL 宏会变成别的意思, 不写。")
            macros.append(_pymol_macro(model, vals))
        return f"select {name}, " + (" or ".join(macros) if macros else "none") + "\n"


class AtomArraySelection_AtomArrayMask(Converter):
    """Builds a boolean mask from an :class:`AtomArraySelection`.

    The mask is aligned with the input array (``True`` at the atoms any rule
    matches), so it can be used directly as ``atom_array[mask]`` and composed with
    the helpers in :mod:`biorazer.structure.selection.mask`.
    """

    def convert(self, atom_array: AtomArray) -> np.ndarray:
        """Select on ``atom_array`` with ``self.input_io``.

        Parameters
        ----------
        atom_array : AtomArray
            The array the patterns are matched against, and the array the
            returned mask is aligned with.

        Returns
        -------
        numpy.ndarray
            1D boolean mask of shape ``atom_array.shape``, ``True`` at every atom
            matched by at least one rule.
        """
        return self.input_io.mask(atom_array)


class AtomArraySelection_AtomArrayIndices(Converter):
    """Builds flat atom indices from an :class:`AtomArraySelection`.

    The order is the one a mask cannot express: rule by rule, and within a rule
    the atom order of the array.  It is the order
    :class:`AtomArraySelection_SelectionCsv` writes into a selection table, so
    indices and the exported table agree row by row.
    """

    def convert(self, atom_array: AtomArray, dedupe: bool = True) -> np.ndarray:
        """Select on ``atom_array`` with ``self.input_io``.

        Parameters
        ----------
        atom_array : AtomArray
            The array the patterns are matched against.
        dedupe : bool
            Keep an atom matched by more than one rule only at its first match
            (the default), or once per matching rule.

        Returns
        -------
        numpy.ndarray
            1D ``int`` array of atom indices.
        """
        return self.input_io.indices(atom_array, dedupe=dedupe)
