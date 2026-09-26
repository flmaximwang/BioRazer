# -*- coding: utf-8 -*-
"""Tests for :mod:`biorazer.structure.objects.selector` 和它的 bridge 转换器。

被测的是"规则表 → 每行一个原子的选择表"这条链路:

* 匹配语法 (通配符/列表/范围/正则/字面值) 与 ``ins_code`` / ``altloc`` 的"空 = 没有"特例;
* GUI 的「模式」↔ csv 字符串 (encode/decode) 往返;
* 行操作 (add_rule / remove_rule / move_rule = GUI 的 加行 / 删行 / 上移下移);
* 规则表 → 选择表 → 规则表 的往返 (四个 csv 转换器), 以及**选择表当规则表打开**时
  每个具体值只命中它自己;
* 展平顺序 (规则序 → 原子序) 与去重语义, mask / indices 两个转换器;
* ``apply`` 对 ``AtomArray`` 与 ``InternalCoord`` 两种目标的判据 (后者要求选择保住连通性);
* GUI 编辑器全流程 (建表/提示/高亮/展开/回写), 没有显示环境时 skip。

夹具不读任何外部结构文件: ``atom_array`` 用一个**手工打包的小数组** (3 个原子, 其中一对同名
原子靠 altloc 区分), ``altloc`` 相关的用例用 ``tests.test_mutation._ideal_chain`` 长一条理想
主链写成 PDB 再经 :class:`biorazer.structure.io.StructureFile_AtomArray` 读回 (``altloc_id`` 标注只有
这条读入路径才有); ``InternalCoord`` 由仓库自己的 bridge 从理想主链转来。
"""

import csv
import pathlib

import numpy as np
import pytest

from biorazer.structure.bridge import (
    AtomArray_InternalCoord,
    AtomArraySelection_AtomArrayMask,
    AtomArraySelection_AtomArrayIndices,
    AtomArraySelection_RuleCsv,
    AtomArraySelection_SelectionCsv,
    RuleCsv_AtomArraySelection,
    SelectionCsv_AtomArraySelection,
)
from biorazer.structure.io import StructureFile_AtomArray
from biorazer.structure.objects import AtomArray, AtomArraySelection
from biorazer.structure.objects.selector import (
    FIELDS,
    decode,
    encode,
    find_columns,
    match_field,
)

from .test_mutation import _ideal_chain

# --------------------------------------------------------------------------
# fixtures / helpers
# --------------------------------------------------------------------------

#: 匹配语法用例: (模式, 值, 字段, 期望命中)
SYNTAX_CASES = [
    ("", "A", "chain", True), ("", "", "ins_code", True), ("", "A", "ins_code", False),
    ("", "", "altloc", True), ("", "A", "altloc", False),
    ("*", "A", "chain", True), ("*", "", "ins_code", True), ("-", "", "ins_code", True),
    ("-", "A", "ins_code", False), ("none", "", "ins_code", True),
    ("-", "", "altloc", True), ("none", "A", "altloc", False), ("*", "B", "altloc", True),
    ("A,B", "B", "chain", True), ("A,B", "C", "chain", False),
    ("A-D", "C", "chain", True), ("A-D", "E", "chain", False), ("A-D", "CA", "chain", False),
    ("A-C", "B", "altloc", True), ("A-C", "D", "altloc", False),
    ("1-10", "10", "resi", True), ("1-10", "11", "resi", False), ("1-10", "1", "resi", True),
    ("-5--1", "-3", "resi", True), ("-5--1", "0", "resi", False), ("-3-5", "-1", "resi", True),
    ("1,5-9,20", "7", "resi", True), ("1,5-9,20", "4", "resi", False),
    ("1,5-9,20", "20", "resi", True),
    ("H1-H20", "H17", "name", True), ("H1-H20", "H21", "name", False),
    ("H1-H20", "CA", "name", False), ("H1-H20", "H", "name", False),
    ("C1-C5", "C3", "name", True), ("C1-C5", "CB", "name", False),
    ("CA-CB", "CA", "name", True), ("CA-CB", "CB", "name", True), ("CA-CB", "CC", "name", False),
    ("C*", "CB", "name", True), ("C*", "N", "name", False), ("?", "A", "chain", True),
    ("re:C[AB]", "CA", "name", True), ("re:C[AB]", "CG", "name", False),
    ("re:C[A-Z]{1,3}", "CGD", "name", True), ("re:C[A-Z]{1,3}", "CGD2", "name", False),
    ("re:1[0-9]", "12", "resi", True), ("re:1[0-9]", "22", "resi", False),
    ("CA", "CA", "name", True), ("CA", "CAA", "name", False),
    ("a", "A", "chain", False),                            # 区分大小写
]


def _tiny_array():
    """手工打包的 3 原子数组: A/1/CA, A/2/CB, B/1A/N (无 altloc 标注)。"""
    arr = AtomArray(3)
    arr.coord = np.zeros((3, 3), float)
    arr.atom_name = np.array(["CA", "CB", "N"], dtype="U4")
    arr.res_id = np.array([1, 2, 1], dtype=np.int32)
    arr.res_name = np.array(["GLY", "GLY", "ALA"], dtype="U3")
    arr.chain_id = np.array(["A", "A", "B"], dtype="U4")
    arr.element = np.array(["C", "C", "N"], dtype="U2")
    arr.ins_code = np.array(["", "", "A"], dtype="U4")
    arr.hetero = np.array([False] * 3)
    return arr


@pytest.fixture
def chain_ic():
    """理想主链 (4 个残基, 只含 N/CA/C/O) 的 :class:`InternalCoord`。"""
    return AtomArray_InternalCoord(_ideal_chain(4, "A")).convert()


# --------------------------------------------------------------------------
# 匹配语法
# --------------------------------------------------------------------------

@pytest.mark.parametrize("pattern, value, field, expect", SYNTAX_CASES)
def test_match_field_table(pattern, value, field, expect):
    """语法表逐条命中判断 (含 ins_code / altloc 的"空 = 没有"特例)。"""
    assert match_field(pattern, value, field) is expect


ENCODE_CASES = [
    # (字段, 模式, 框里输入, 值, 期望命中)
    ("ins_code", "任意", "", "", True), ("ins_code", "任意", "", "A", True),
    ("ins_code", "字面值", "", "", True), ("ins_code", "字面值", "", "A", False),
    ("altloc", "任意", "", "", True), ("altloc", "任意", "", "B", True),
    ("altloc", "字面值", "", "", True), ("altloc", "字面值", "", "B", False),
    ("chain", "任意", "", "Z", True), ("chain", "字面值", "A", "A", True),
    ("chain", "字面值", "A,B", "A,B", True),               # 含语法字符 → 转义成正则, 仍是字面值
    ("chain", "字面值", "A,B", "A", False),
    ("chain", "列表", "A,B", "B", True), ("chain", "列表", "A,B", "C", False),
    ("name", "字面值", "CA-CB", "CC", False),              # 字面值不当范围用 ("CC" 在范围里会命中)
    ("resi", "范围", "1-10", "7", True), ("resi", "范围", "1-10", "11", False),
    ("name", "通配符", "C*", "CB", True), ("name", "通配符", "C*", "N", False),
    ("name", "正则", "C[AB]", "CA", True), ("name", "正则", "C[AB]", "CG", False),
    ("altloc", "列表", "A,B", "B", True), ("altloc", "范围", "A-C", "B", True),
]


@pytest.mark.parametrize("field, mode, text, value, expect", ENCODE_CASES)
def test_encode_then_match(field, mode, text, value, expect):
    """GUI 模式 → csv 字符串 → 命中判断。"""
    assert match_field(encode(field, mode, text), value, field) is expect


@pytest.mark.parametrize("pattern", ["", "*", "re:C[AB]", "A,B", "1-10", "H1-H20", "CA",
                                     "-5", "none"])
@pytest.mark.parametrize("field", FIELDS)
def test_encode_decode_roundtrip(field, pattern):
    """手写/回写的 csv 字符串原样保留 (encode ∘ decode = id)。"""
    assert encode(field, *decode(field, pattern)) == pattern


def test_fields_include_altloc():
    """altloc 是第 5 个匹配字段, 表头别名认 biotite 的 altloc_id。"""
    assert FIELDS == ("ins_code", "chain", "resi", "name", "altloc")
    assert find_columns(["ins_code", "chain_id", "res_id", "atom_name", "altloc_id"]) == {
        f: i for i, f in enumerate(FIELDS)}
    with pytest.raises(ValueError):
        find_columns(["chain", "resi", "name", "altloc"])     # 缺 ins_code 就报错, 不猜
    with pytest.raises(ValueError):
        find_columns(["ins_code", "chain", "resi", "name"])   # 缺 altloc 也一样


# --------------------------------------------------------------------------
# dataclass 契约
# --------------------------------------------------------------------------

def test_dataclass_contract():
    """None → 空; 短行补到表头宽; 表头缺列当场报错; 相等按值比。"""
    empty = AtomArraySelection()
    assert empty.header == list(FIELDS) and empty.rules == []

    sel = AtomArraySelection(rules=[["", "A"]], header=list(FIELDS))
    assert sel.rules == [["", "A", "", "", ""]]
    assert sel == AtomArraySelection(rules=[["", "A"]])
    assert repr(sel) == f"AtomArraySelection(1 rules, header={list(FIELDS)})"

    with pytest.raises(ValueError):
        AtomArraySelection(rules=[["", "A", "1", "CA"]], header=["ins_code", "chain", "resi"])

    # 比表头还长的行 = csv 里那个逗号没加引号, 报错而不是静默丢单元格
    with pytest.raises(ValueError) as excinfo:
        AtomArraySelection(rules=[["", "A", "1", "CA", "", "批注", "多出来的"]])
    assert "more cells than the header" in str(excinfo.value)


def test_tkinter_is_not_imported_with_the_selector():
    """编辑器是惰性 import 的: ``import objects`` 不该把 tkinter 拉进来 (无显示环境照样能用库)。"""
    import os
    import subprocess
    import sys

    repo_root = pathlib.Path(__file__).resolve().parents[1]
    code = ("import sys, biorazer.structure.objects as O;"
            "print('tkinter' in sys.modules);"
            "print(callable(O.AtomArraySelection))")
    env = dict(os.environ)
    env["PYTHONPATH"] = str(repo_root) + os.pathsep + env.get("PYTHONPATH", "")
    proc = subprocess.run([sys.executable, "-c", code], cwd=repo_root, env=env,
                          capture_output=True, text=True, check=True)
    assert proc.stdout.split() == ["False", "True"], proc.stdout


# --------------------------------------------------------------------------
# 匹配 / 展平 / mask
# --------------------------------------------------------------------------

def test_resolve_reports_hits_and_warnings():
    """每行命中数 + 三类警告: 0 命中 / 链不在目标里 / 正则错误 / 五元组重复。"""
    arr = _tiny_array()
    sel = AtomArraySelection(rules=[
        ["", "A", "1-2", "C*", ""],                      # 命中 0, 1
        ["A", "B", "1", "N", ""],                        # 命中 2 (ins_code A)
        ["", "Z", "1", "CA", ""],                        # 0 命中 + 链不存在
        ["", "A", "1", "re:[", ""],                      # 正则错误
        ["", "A", "1", "CA", "A"],                       # 目标没有 altloc 标注
    ])
    hits, warns, dup = sel.resolve(arr)
    assert hits == [[0, 1], [2], [], [], []]
    assert warns[2] == ["0 命中", "chain 'Z' 不在目标里 (只有 A,B)"]
    assert warns[3][0].startswith("正则错误")
    assert warns[4] == ["目标没有 altloc 标注 (InternalCoord 不带; 手工拼的 AtomArray 也不带; "
                        "读文件请用 StructureFile_AtomArray), altloc 模式无法判定", "0 命中"]
    assert dup == 0


def test_flat_order_dedupe_and_mask():
    """展平顺序 = 规则序 → 原子序; 去重只留第一次命中的那条规则; mask 与顺序无关。"""
    arr = _tiny_array()
    sel = AtomArraySelection(rules=[
        ["", "A", "1-2", "C*", ""],                      # 0, 1
        ["", "A", "1", "CA", ""],                        # 0 (重复)
        ["A", "B", "1", "N", ""],                        # 2 (那个原子带插入码 A)
    ])
    assert sel.pairs(arr) == [(0, 0), (0, 1), (2, 2)]
    assert sel.pairs(arr, dedupe=False) == [(0, 0), (0, 1), (1, 0), (2, 2)]
    assert sel.indices(arr).tolist() == [0, 1, 2]
    assert sel.indices(arr, dedupe=False).tolist() == [0, 1, 0, 2]
    assert sel.mask(arr).tolist() == [True, True, True]
    _, _, dup = sel.resolve(arr)
    assert dup == 1


def test_mask_and_indices_from_the_bridge():
    """两个 A_B 转换器 = 对象上同名查询的一层 "方向" 命名。"""
    arr = _tiny_array()
    sel = AtomArraySelection(rules=[["", "A", "1", "CA", ""], ["A", "B", "1", "N", ""]])
    mask = AtomArraySelection_AtomArrayMask(input_io=sel).convert(arr)
    idx = AtomArraySelection_AtomArrayIndices(input_io=sel).convert(arr)
    assert mask.tolist() == sel.mask(arr).tolist() == [True, False, True]
    assert idx.tolist() == sel.indices(arr).tolist() == [0, 2]
    assert (arr[mask].atom_name == arr[idx].atom_name).all()


# --------------------------------------------------------------------------
# csv
# --------------------------------------------------------------------------

def test_csv_roundtrip_and_selection_table(tmp_path):
    """规则表往返一致; 选择表严格每行一个原子 (具体值 + 去重), 且能当规则表打开。"""
    arr = _tiny_array()
    sel = AtomArraySelection(rules=[
        ["", "A", "1-2", "C*", "*", "规则一"],
        ["", "A", "1", "CA", "*", "规则二(与规则一重复)"],
        ["", "Z", "1", "CA", "*", "坏链"],
    ], header=list(FIELDS) + ["note"])

    rules_path = tmp_path / "rules.csv"
    AtomArraySelection_RuleCsv(output_io=rules_path).write(sel)
    assert RuleCsv_AtomArraySelection(input_io=rules_path).read() == sel
    # 0 命中的行也原样回写 —— 规则表是唯一真源, 展开时才丢
    assert [r[-1] for r in
            RuleCsv_AtomArraySelection(input_io=rules_path).read().rules] == [
        "规则一", "规则二(与规则一重复)", "坏链"]

    out_path = tmp_path / "selection.csv"
    AtomArraySelection_SelectionCsv(output_io=out_path).write(sel, arr)
    rows, header = list(csv.reader(out_path.open())), None
    header, body = rows[0], rows[1:]
    assert header == list(FIELDS) + ["note"]
    assert [r[:5] for r in body] == [["", "A", "1", "CA", ""], ["", "A", "2", "CB", ""]]
    assert [r[-1] for r in body] == ["规则一", "规则一"]      # 重复的那条被去重丢掉

    AtomArraySelection_SelectionCsv(output_io=out_path).write(sel, arr, dedupe=False)
    body = list(csv.reader(out_path.open()))[1:]
    assert len(body) == 3 and body[2][-1] == "规则二(与规则一重复)"

    # 选择表当规则表打开: 每个具体值只命中它自己
    back = SelectionCsv_AtomArraySelection(input_io=out_path).read()
    hits, warns, dup = back.resolve(arr)
    assert [h for h in hits] == [[0], [1], [0]]
    assert dup == 1                                          # 第三行与第一行同一个原子

    with pytest.raises(ValueError):
        # 展开需要目标: 选择表就是规则表按某个结构展平的结果, 没结构没法展
        AtomArraySelection_SelectionCsv(output_io=tmp_path / "x.csv").write(sel, None)


def test_row_ops_mirror_the_editor():
    """add_rule / remove_rule / move_rule = GUI 的 加行 / 删行 / 上移下移。"""
    arr = _tiny_array()
    sel = AtomArraySelection(header=list(FIELDS) + ["note"])
    assert sel.add_rule() == 0 and sel.rules == [[""] * 6]   # 加一整行空的 (命中一切)
    # ins_code 的空 = "没有插入码", 要"任意"得写 "*" (目标里那个 N 带插入码 A)
    assert sel.add_rule({"chain": "B", "name": "N", "ins_code": "*"}) == 1
    # 列名认表头名与字段别名; 值是规则表里那格字符串, 不是 GUI 的 (模式, 输入) 对
    assert sel.add_rule({"atom_name": "CA", "auth_asym_id": "A", "note": "A 链 CA"},
                        index=0) == 0
    assert sel.rules == [["", "A", "", "CA", "", "A 链 CA"],
                         [""] * 6,
                         ["*", "B", "", "N", "", ""]]
    assert sel.indices(arr, dedupe=False).tolist() == [0, 0, 1, 2]   # 展平顺序 = 规则序

    # 上移/下移: 动的是"哪个原子由哪条规则写进选择表"
    assert sel.move_rule(0, 1) == 1 and sel.move_rule(2, 1) == 2     # 越界不动
    assert sel.rules[1][:5] == ["", "A", "", "CA", ""]
    assert sel.indices(arr, dedupe=False).tolist() == [0, 1, 0, 2]

    assert sel.remove_rule(0) == [""] * 6                            # 删行返回被删的那行
    assert len(sel.rules) == 2 and sel.remove_rule(-1)[-2:] == ["", ""]   # 负号从末尾数
    with pytest.raises(IndexError):
        sel.remove_rule(5)
    with pytest.raises(ValueError) as excinfo:
        sel.add_rule({"no_such_column": "x"})
    assert "unknown column" in str(excinfo.value)


def test_selection_csv_read_escapes_literals(tmp_path):
    """选择表里的具体值含语法字符时按字面值转义, 不会被当列表/范围/正则。"""
    path = tmp_path / "sel.csv"
    path.write_text("ins_code,chain,resi,name,altloc\n"
                    ",A,1,\"CA,CB\",\n"
                    ",A,1,\"re:x\",\n"
                    ",A,1,\"CA-CB\",\n")
    sel = SelectionCsv_AtomArraySelection(input_io=path).read()
    assert [row[:5] for row in sel.rules] == [
        ["", "A", "1", "re:CA,CB", ""],
        ["", "A", "1", "re:re:x", ""],
        ["", "A", "1", "re:CA\\-CB", ""]]
    # 原样当规则表读则不是字面值: "CA,CB" 是列表 (命中 CA 或 CB)
    plain = AtomArraySelection(rules=[["", "A", "1-2", "CA,CB", ""]])
    arr = _tiny_array()
    assert plain.indices(arr).tolist() == [0, 1]
    assert sel.indices(arr).tolist() == []


# --------------------------------------------------------------------------
# altloc
# --------------------------------------------------------------------------

def test_altloc_field_matches_the_read_path(altloc_pdb):
    """altloc 的值来自 biotite 的 altloc_id (只有 altloc="all" 读入才有)。"""
    arr = StructureFile_AtomArray(input_io=altloc_pdb).read()
    assert len(arr) == 13
    expected_alt = [" "] * 13
    expected_alt[5], expected_alt[12] = "A", "B"        # 第 2 残基 CA (在中间) + 追加的那份
    assert arr.altloc_id.tolist() == expected_alt
    dup_ca = (arr.res_id == 2) & (arr.atom_name == "CA")
    assert dup_ca.sum() == 2 and set(arr.altloc_id[dup_ca]) == {"A", "B"}

    def ca(indices):
        return sorted(i for i in indices if arr.atom_name[i] == "CA")

    both = AtomArraySelection(rules=[["", "A", "1-3", "CA", "*"]])
    only_a = AtomArraySelection(rules=[["", "A", "1-3", "CA", "A"]])
    only_b = AtomArraySelection(rules=[["", "A", "1-3", "CA", "B"]])
    none_alt = AtomArraySelection(rules=[["", "A", "1-3", "CA", ""]])
    assert len(ca(both.indices(arr))) == 4                   # 3 个残基 + 重复的那个
    assert ca(only_a.indices(arr)) == [5] and ca(only_b.indices(arr)) == [12]
    assert [arr.res_id[i] for i in ca(only_a.indices(arr))] == [2]
    assert (ca(only_a.indices(arr)) + ca(only_b.indices(arr))
            == sorted(set(ca(both.indices(arr))) - set(ca(none_alt.indices(arr)))))
    # 带 altloc 的规则在无标注的目标上报 "无法判定"
    _, warns, _ = only_a.resolve(_tiny_array())
    assert warns[0] == ["目标没有 altloc 标注 (InternalCoord 不带; 手工拼的 AtomArray 也不带; "
                        "读文件请用 StructureFile_AtomArray), altloc 模式无法判定", "0 命中"]


# --------------------------------------------------------------------------
# apply
# --------------------------------------------------------------------------

def test_apply_atomarray_is_masked_indexing():
    """apply(AtomArray) = 布尔下标切片 (含坐标与注解)。"""
    arr = _tiny_array()
    sel = AtomArraySelection(rules=[["", "A", "2", "*", ""]])
    sub = sel.apply(arr)
    assert isinstance(sub, AtomArray)
    assert sub.atom_name.tolist() == ["CB"] and sub.res_id.tolist() == [2]
    assert (sub.coord == arr[sel.mask(arr)].coord).all()
    assert (sub.atom_name == arr[sel.mask(arr)].atom_name).all()


def test_apply_internalcoord_keeps_coordinates(chain_ic):
    """整链选择: 结果原子数一致, 坐标逐原子与输入相同 (对拍, 不是单测)。"""
    full = AtomArraySelection(rules=[["", "A", "*", "*", ""]])
    sub = full.apply(chain_ic)
    assert len(sub) == len(chain_ic) == 16
    before, after = chain_ic.to_coords(), sub.to_coords()
    assert max(float(np.linalg.norm(np.asarray(after[i], float) - np.asarray(before[i], float)))
               for i in range(len(sub))) < 1e-9
    assert [sub.atom_repr(i) for i in range(len(sub))] == \
           [chain_ic.atom_repr(i) for i in range(len(chain_ic))]

    # 残基 1-2 是一个自足的片段 (anchor 在残基 1), 8 个原子, 仍能长出来
    head = AtomArraySelection(rules=[["", "A", "1-2", "*", ""]]).apply(chain_ic)
    assert len(head) == 8 and len(head.to_coords()) == 8
    xyz = chain_ic.to_coords()
    assert all(np.allclose(head.to_coords()[i], xyz[i]) for i in range(8))


def test_apply_internalcoord_rejects_a_cut_frame(chain_ic):
    """断了 anchor 的选择: 原子还在图里但长不出来 → ValueError 列出裸露原子。"""
    sel = AtomArraySelection(rules=[["", "A", "2-4", "*", ""]])
    with pytest.raises(ValueError) as excinfo:
        sel.apply(chain_ic)
    msg = str(excinfo.value)
    assert "without a bond connection" in msg
    assert "A:2:ALA:N" in msg and "A:4:ALA:O" in msg
    # 同一选择在 AtomArray 上照常切片 —— 约束只属于 IC 的生长图
    assert len(sel.apply(_ideal_chain(4, "A"))) == 12


def test_atom_array_stack_is_rejected():
    """2 维目标 (AtomArrayStack) 直接 TypeError —— 布尔下标在它上面切的是 model 轴。"""
    from biotite.structure import AtomArrayStack

    arr = _tiny_array()
    stacked = AtomArrayStack(2, len(arr))                   # (depth, length): coord 是 2 维
    stacked.coord = np.stack([arr.coord, arr.coord + 1.0])
    for category in ("atom_name", "chain_id", "res_id", "res_name", "element", "ins_code",
                     "hetero"):
        setattr(stacked, category, getattr(arr, category))
    sel = AtomArraySelection(rules=[["", "A", "*", "*", ""]])
    with pytest.raises(TypeError) as excinfo:
        sel.apply(stacked)
    assert "single model" in str(excinfo.value)
    assert len(sel.apply(stacked[0])) == 2                   # 取一个 model 就照常能用


def test_apply_leaves_the_input_alone(chain_ic):
    """apply 不就地改输入 (IC 子集是新建的对象, 连通表重新编号)。"""
    full = AtomArraySelection(rules=[["", "A", "1-3", "*", ""]])
    sub = full.apply(chain_ic)
    assert len(chain_ic) == 16 and len(sub) == 12
    assert max(chain_ic.dihedra) == max(q for q in chain_ic.dihedra)
    assert all(i < 12 for q in sub.dihedra for i in q)


# --------------------------------------------------------------------------
# GUI 编辑器 (没有显示环境就 skip)
# --------------------------------------------------------------------------

def _need_display():
    tk = pytest.importorskip("tkinter")
    try:
        root = tk.Tk()
    except tk.TclError as exc:                               # 无 DISPLAY 的 CI
        pytest.skip(f"no display for tkinter: {exc}")
    root.destroy()
    return tk


def test_editor_full_flow(tmp_path, altloc_pdb):
    """编辑器全流程: 建表 / 灰字提示 / 激活行 / 问题行 / 展开 / 回写两个 csv / 行操作。

    结构 = :func:`altloc_pdb` (3 残基主链 + 第 2 残基 CA 的 A/B 两个构象, 13 原子),
    规则 5 行: 命中 / 0 命中 / 命中 / 坏链 / 与第 1 行重复。
    """
    _need_display()
    from biorazer.structure.objects._selector_gui import (
        BAND_ACTIVE,
        BAND_PROBLEM,
        BG_PROBLEM,
    )

    rules_path = tmp_path / "rules.csv"
    rules_path.write_text(
        "ins_code,chain,resi,name,altloc,note\n"
        ",A,1-3,C*,*,规则一\n"
        "A,A,2,CA,,\"规则二(该残基没有插入码, 0 命中)\"\n"   # 批注含逗号 → csv 要引号
        ",A,1,N,*,规则三(N 不在 C* 里)\n"
        ",Z,1,CA,*,规则四(坏链)\n"
        ",A,2,CA,A,规则五(与规则一重复)\n")
    out_path = tmp_path / "selection.csv"

    sel = AtomArraySelection()
    seen = {}

    def probe(api):
        state, selector = api["state"], api["selector"]
        assert len(selector.rules) == 5 and len(selector.header) == 6
        assert [len(row) for row in state["cells"]] == [6] * 5

        # 窗口钉在亮色外观上 (macOS 默认 auto = 跟随系统深色)
        root = api["root"]
        if root.tk.call("tk", "windowingsystem") == "aqua":
            seen["appearance"] = root.tk.call("tk::unsupported::MacWindowStyle",
                                              "appearance", root._w)

        seen["start"] = (state["focus"], [n.cget("text") for n in state["rownums"]])
        # 灰字提示: 空框给"该模式下该怎么写", 换模式跟着换, 聚焦清掉
        cell = state["cells"][1][4]                          # 第 2 行的 altloc 格 (模式空)
        seen["hints"] = [cell.mode.get(), cell.entry.get(), cell.get()]
        cell.mode.set("任意")
        cell.refresh_hint()
        seen["hints_any"] = [cell.entry.get(), cell.get()]
        cell.mode.set("字面值")
        cell.refresh_hint()
        cell._focus_in()
        seen["cleared"] = cell.text.get()
        cell._focus_out()
        seen["active_after_hint"] = [i for i, f in enumerate(state["rowframes"])
                                     if f.cget("bg") == BAND_ACTIVE]
        state["cells"][0][0]._focus_in()                     # 先点回第 1 行 (它的格子不红)

        api["check_all"]()                                   # 解析全部 → 红带/红底
        seen["problem_bands"] = [i for i, f in enumerate(state["rowframes"])
                                 if f.cget("bg") == BAND_PROBLEM]
        seen["problem_cells"] = sum(1 for c in state["cells"] for x in c
                                    if x.entry.cget("bg") == BG_PROBLEM)
        seen["active"] = [i for i, f in enumerate(state["rowframes"])
                          if f.cget("bg") == BAND_ACTIVE]
        seen["rownums"] = [n.cget("text") for n in state["rownums"]]
        seen["status"] = api["status_var"].get()

        state["cells"][4][0]._focus_in()                     # 点第 5 行 → 激活行跟着换
        seen["active_after_click"] = [i for i, f in enumerate(state["rowframes"])
                                      if f.cget("bg") == BAND_ACTIVE]
        seen["rownums_after_click"] = [n.cget("text") for n in state["rownums"]]

        api["export"](str(out_path))                          # 显式路径: 不弹确认框
        seen["export"] = (list(csv.reader(open(out_path)))[1:], api["out_var"].get())

        api["dedupe_var"].set(False)                          # 关掉去重 → 重复照写
        api["export"](str(tmp_path / "raw.csv"))
        seen["raw"] = list(csv.reader(open(tmp_path / "raw.csv")))[1:]
        api["dedupe_var"].set(True)

        api["save_rules"](str(tmp_path / "rules_roundtrip.csv"))   # 规则表回写 (唯一真源)
        seen["rules_back"] = list(csv.reader(open(tmp_path / "rules_roundtrip.csv")))

        n0 = len(selector.rules)                              # 行操作
        api["add_row"]()
        api["move_row"](-1)
        api["del_row"]()
        seen["row_ops"] = (n0, len(selector.rules))
        seen["cells_after"] = [len(r) for r in state["cells"]]

    sel.run_editor(structure_file=str(altloc_pdb), rule_csv=str(rules_path),
                   selection_csv=str(out_path), on_ready=probe)

    assert seen["hints"] == ["字面值", "留空=无 altloc", ""], seen["hints"]
    assert seen["appearance"] == "aqua", seen["appearance"]   # 固定亮色, 不跟系统深色
    assert seen["hints_any"] == ["*=任意 altloc", "*"], seen["hints_any"]
    assert seen["cleared"] == ""
    assert seen["start"] == (0, ["▶ 1", "2", "3", "4", "5"]), seen["start"]
    # 点过第 2 行的格子 → 激活行跟到那一行 (点任意一格即选中该行)
    assert seen["active_after_hint"] == [1]
    assert seen["active"] == [0] and seen["rownums"][0] == "▶ 1"
    assert seen["active_after_click"] == [4]
    assert seen["rownums_after_click"] == ["1", "2", "3", "4", "▶ 5"]
    # 激活行是蓝带, 会盖掉它自己的红带: 所以第 2 行的红带要在它不被激活时才看得到
    assert seen["problem_bands"] == [1, 3]                   # 0 命中 / 坏链
    # 底色刷在输入框上 (下拉框不刷): 每个问题行的 6 列各有一个 entry
    assert seen["problem_cells"] == 2 * len(FIELDS) + 2 * 1
    assert "2 行有问题" in seen["status"], seen["status"]

    rows, out_var = seen["export"]                           # 选择表: 每行一个原子, 已去重
    assert len(rows) == 8 and all(len(r) == 6 for r in rows)
    assert len({tuple(r[:5]) for r in rows}) == len(rows), "去重后不该有重复原子"
    notes = [r[-1] for r in rows]
    assert notes.count("规则一") == 7 and notes[-1] == "规则三(N 不在 C* 里)"
    assert set(r[4] for r in rows if r[3] == "CA") == {"A", "B", ""}   # 构象列是具体值
    assert str(out_path) in out_var and "8 个原子" in out_var

    assert len(seen["raw"]) == 9, seen["raw"]                # 关掉去重后重复照写
    assert [r[-1] for r in seen["raw"]].count("规则五(与规则一重复)") == 1

    rules_back = seen["rules_back"]                          # 规则表回写: 行/列都不变
    assert rules_back[0] == ["ins_code", "chain", "resi", "name", "altloc", "note"]
    assert [r[-1] for r in rules_back[1:]] == [
        "规则一", "规则二(该残基没有插入码, 0 命中)", "规则三(N 不在 C* 里)",
        "规则四(坏链)", "规则五(与规则一重复)"]
    assert seen["row_ops"] == (5, 5) and seen["cells_after"] == [6] * 5


# --------------------------------------------------------------------------
# 选择器吃 InternalCoord 目标 / bridge 拒绝带 altloc 的数组
# --------------------------------------------------------------------------

class TestSelectorOnInternalCoord:
    """IC 的 record 不带 altloc, 所以 IC 目标与手工数组一样是"没有这个 category"。"""

    def test_altloc_rule_on_ic_cannot_tell(self):
        """带 altloc 约束的规则在 IC 上报"无法判定", 而不是静默当命中。"""
        ic = AtomArray_InternalCoord(input_io=_tiny_array()).convert()
        sel = AtomArraySelection(rules=[["", "A", "1", "CA", "A"]])
        hits, warns, _dup = sel.resolve(ic)
        assert hits == [[]]
        assert warns[0] == ["目标没有 altloc 标注 (InternalCoord 不带; 手工拼的 AtomArray 也不带; "
                            "读文件请用 StructureFile_AtomArray), altloc 模式无法判定", "0 命中"]
        # 不带 altloc 约束的规则在 IC 上照常命中, tag 也不带构象后缀
        hits2 = AtomArraySelection(rules=[["", "A", "1", "CA", ""]]).resolve(ic)[0]
        assert hits2 == [[0]] and ic.atom_repr(hits2[0][0]) == "A:1:GLY:CA"


class TestBridgeRefusesAltloc:
    """IC 是"一个原子名一格"的生长树 —— 带 altloc 的数组不给建 (会静默丢副本)。"""

    def test_array_with_altloc_is_refused(self, altloc_pdb):
        arr = StructureFile_AtomArray(input_io=altloc_pdb).read()
        with pytest.raises(ValueError) as exc:
            AtomArray_InternalCoord(input_io=arr).convert()
        msg = str(exc.value)
        assert "carries alternate conformations (2 of 13 atoms" in msg
        assert "A/2/ALA/CA altloc 'A'" in msg          # 点名第一个带标签的原子
        assert "InternalCoord has no altloc" in msg

    def test_array_without_altloc_converts(self):
        """没有 altloc 标注的数组照常建 IC —— guard 只挡带标签的。"""
        from tests.test_mutation import _ideal_chain
        arr = _ideal_chain(3, "A")
        assert not hasattr(arr, "altloc_id")
        ic = AtomArray_InternalCoord(input_io=arr).convert()
        assert len(ic.atoms) == len(arr) == 12
        assert not hasattr(ic, "altloc_id")


# --------------------------------------------------------------------------
# CLI 子命令 (biorazer select)
# --------------------------------------------------------------------------

class TestSelectorCli:
    """``biorazer select`` 的注册与 headless 运行 (parser/runner 在 selector 模块)。"""

    @staticmethod
    def _parser():
        import argparse

        from biorazer.structure.objects import cli as objects_cli

        parser = argparse.ArgumentParser()
        sub = parser.add_subparsers(dest="command")
        objects_cli.register_subcommand(sub)
        return parser, sub

    def test_registered(self):
        """子命令挂上去了, func 指向 runner。"""
        from biorazer.structure.objects import selector

        _, sub = self._parser()
        assert "select" in sub.choices
        assert sub.choices["select"].get_default("func") is selector._run_selector

    def test_headless_run(self, tmp_path, altloc_pdb, capsys):
        """--print --export: 逐行报命中数, 并把选择表写到 -o (不建窗)。"""
        parser, _ = self._parser()
        rules = tmp_path / "rules.csv"
        rules.write_text("ins_code,chain,resi,name,altloc,note\n"
                         ",A,1-3,C*,*,规则一\n"
                         ",Z,1,CA,*,坏链\n")
        out = tmp_path / "selection.csv"
        args = parser.parse_args(["select", "-s", str(altloc_pdb), "--rules", str(rules),
                                  "-o", str(out), "--print", "--export"])
        args.func(args)

        printed = capsys.readouterr().out
        assert "altloc A,B" in printed and "重复五元组 0" in printed      # 结构摘要
        assert "行 1" in printed and "7 原子" in printed
        assert "⚠ 0 命中" in printed and "chain 'Z' 不在目标里 (只有 A)" in printed
        assert "--- 规则 2 行 → 展平 7 个原子, 重复 0, 1 行有问题" in printed

        rows = list(csv.reader(open(out)))
        assert rows[0] == list(FIELDS) + ["note"]
        assert len(rows) == 8 and rows[1] == ["", "A", "1", "CA", "", "规则一"]

    def test_print_needs_a_structure(self):
        """--print 没有 -s 直接报错 (没有结构就没法展开)。"""
        parser, _ = self._parser()
        args = parser.parse_args(["select", "--print"])
        with pytest.raises(SystemExit):
            args.func(args)

    def test_top_level_cli_registers_it(self):
        """biorazer.cli 的门面把 select 一起注册 (漏一行就静默少一个子命令)。"""
        import argparse

        from biorazer import cli as top_cli

        parser = argparse.ArgumentParser()
        sub = parser.add_subparsers(dest="command")
        top_cli.access_cli.register_subcommand(sub)
        top_cli.colabfold_cli.register_subcommand(sub)
        top_cli.plot_cli.register_subcommand(sub)
        top_cli.report_cli.register_subcommand(sub)
        top_cli.selector_cli.register_subcommand(sub)
        assert "select" in sub.choices
