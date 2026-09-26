# -*- coding: utf-8 -*-
"""原子选择编辑器的窗口 (tkinter) —— :mod:`~biorazer.structure.objects.selector` 的私有模块。

入口只有一个: :meth:`~biorazer.structure.objects.selector.AtomArraySelection.run_editor`,
它惰性 import 本模块并把选择实例交进来 (tkinter 因此不进 ``import objects`` 的路径)。
编辑器**就地**改传进来的选择 (``rules`` / ``header``), 边上给两个 csv:
起始规则表 (:class:`~biorazer.structure.bridge.RuleCsv_AtomArraySelection`) 与导出的选择表
(:class:`~biorazer.structure.bridge.AtomArraySelection_SelectionCsv`)。结构走库自己的读入口
(:class:`biorazer.structure.io.StructureFile_AtomArray`, 按后缀走 ``Pdb_AtomArray`` / ``Cif_AtomArray``)。
行操作 (加行/删行/上移下移) 直接调 :class:`~biorazer.structure.objects.selector.AtomArraySelection`
的同名方法, 不另写一套。

GUI 每个字段是「模式下拉 + 输入框」; 输入框为空时框里显示该模式下该怎么写 (灰字提示,
聚焦即清)。表格配色: 有问题行红、激活行淡蓝。

窗口**固定亮色**: macOS 上 Tk 默认跟随系统外观 (``appearance=auto``) 而默认配色的控件
用的是语义颜色, 系统深色下整窗跟着变黑 —— 见 :func:`force_light_appearance`。
"""

from __future__ import annotations

from pathlib import Path

from biorazer.structure.bridge import (
    AtomArraySelection_RuleCsv,
    AtomArraySelection_SelectionCsv,
    RuleCsv_AtomArraySelection,
)
from biorazer.structure.io import StructureFile_AtomArray
from biorazer.structure.objects.selector import (
    MODES,
    AtomArraySelection,
    _TargetView,
    decode,
    encode,
    match_rule,
)

# 表格配色: 激活行 / 有问题行 (红) / 平常
BG_IDLE, BG_PROBLEM, BG_ACTIVE = "white", "#ffd9d9", "#e8f0fe"
BAND_IDLE, BAND_ACTIVE, BAND_PROBLEM = "#e0e0e0", "#1a73e8", "#e8a0a0"
FG_ROW_IDLE, FG_ROW_ACTIVE = "#888888", "#1a73e8"

# 空框里显示的提示 = 该字段在该模式下该怎么写
HINTS = {
    "ins_code": {"任意": "*=任意插入码", "字面值": "留空=无插入码", "通配符": "*",
                 "列表": "A,B", "范围": "A-D", "正则": "[AB]"},
    "chain": {"任意": "留空=任意链", "字面值": "A", "通配符": "A*", "列表": "A,B",
              "范围": "A-D", "正则": "[AB]"},
    "resi": {"任意": "留空=任意号", "字面值": "45", "通配符": "4*", "列表": "1,5-9",
             "范围": "1-10", "正则": "1[0-9]"},
    "name": {"任意": "留空=任意原子", "字面值": "CA", "通配符": "C*", "列表": "CA,CB",
             "范围": "H1-H20", "正则": "C[AB]"},
    "altloc": {"任意": "*=任意 altloc", "字面值": "留空=无 altloc", "通配符": "*",
               "列表": "A,B", "范围": "A-D", "正则": "[AB]"},
}


class FieldCell:
    """规则表里一个字段: [模式下拉] [输入框]。空框显示该模式下该怎么写的提示。"""

    def __init__(self, parent, field, pattern, on_change):
        import tkinter as tk
        from tkinter import ttk

        self.field, self.on_change, self._focused, self._hint = field, on_change, False, False
        mode, text = decode(field, pattern)
        self.mode, self.text = tk.StringVar(value=mode), tk.StringVar(value=text)
        self.combo = ttk.Combobox(parent, textvariable=self.mode, values=MODES,
                                  width=6, state="readonly")
        self.entry = tk.Entry(parent, textvariable=self.text, width=15)
        self.widgets = [self.combo, self.entry]
        self.combo.bind("<<ComboboxSelected>>", lambda _e: (self.refresh_hint(), self.on_change()))
        self.combo.bind("<FocusIn>", lambda _e: self.on_change())   # 点下拉也算选中该行
        self.entry.bind("<KeyRelease>", lambda _e: self.on_change())
        self.entry.bind("<FocusIn>", self._focus_in)
        self.entry.bind("<FocusOut>", self._focus_out)
        self.refresh_hint()

    def _focus_in(self, _e=None):
        self._focused = True
        if self._hint:                                          # 灰字是提示, 不是输入
            self.text.set("")
            self._hint = False
            self.entry.configure(fg="black")
        self.on_change()

    def _focus_out(self, _e=None):
        self._focused = False
        self.refresh_hint()

    def refresh_hint(self):
        if self._focused:
            return
        if self._hint or not self.text.get().strip():           # 空框 / 框里现在是提示 → 写本模式的提示
            self.text.set(HINTS[self.field][self.mode.get()])
            self._hint = True
            self.entry.configure(fg="grey")

    def get(self) -> str:
        """该字段写进规则表的字符串 (提示文字不算输入)。"""
        return encode(self.field, self.mode.get(), "" if self._hint else self.text.get())

    def paint(self, problem: bool, active: bool):
        """背景色: 有问题 → 红; 激活行 → 淡蓝; 否则白。"""
        self.entry.configure(bg=BG_PROBLEM if problem else (BG_ACTIVE if active else BG_IDLE))


class PlainCell:
    """字段列之外的列: 原样显示, 原样回写。"""

    def __init__(self, parent, text, on_change):
        import tkinter as tk

        self.text = tk.StringVar(value=text)
        self.entry = tk.Entry(parent, textvariable=self.text, width=16)
        self.widgets = [self.entry]
        self.entry.bind("<KeyRelease>", lambda _e: on_change())
        self.entry.bind("<FocusIn>", lambda _e: on_change())

    def get(self) -> str:
        return self.text.get()

    def paint(self, problem: bool, active: bool):
        self.entry.configure(bg=BG_PROBLEM if problem else (BG_ACTIVE if active else BG_IDLE))


def force_light_appearance(window):
    """把窗口外观钉成亮色 (macOS 的 NSAppearance), 返回设成了什么; 别的平台返回 ``None``。

    Tk 在 macOS 上默认 ``appearance=auto`` (= 跟随系统), 而默认配色的控件用的是
    ``systemWindowBackgroundColor`` / ``systemTextColor`` 这类**语义颜色** —— 系统深色下
    整个窗口 (连下拉框) 一起变黑, 编辑器里那些写死的白/红/蓝底色反倒成了异类。
    这里把窗口外观钉成 ``aqua`` (= 亮色), 子控件继承同一个 NSAppearance, 配色一并回来。

    实测 (Tk 8.6.13): ``tk::unsupported::MacWindowStyle appearance .`` 默认回答 ``auto``,
    设成 ``aqua`` 后读回 ``aqua``; 老 Tk 没这个子命令, 那就保持系统外观 (不报错)。
    """
    if window.tk.call("tk", "windowingsystem") != "aqua":
        return None
    try:
        window.tk.call("tk::unsupported::MacWindowStyle", "appearance", window._w, "aqua")
    except window.tk.TclError:
        return None
    return "aqua"


def run_editor(selector: AtomArraySelection, structure_file=None, rule_csv=None,
               selection_csv=None, dedupe: bool = True, on_ready=None):
    """打开编辑器窗口; 见 :meth:`AtomArraySelection.run_editor`。

    Parameters
    ----------
    selector : AtomArraySelection
        就地编辑的选择。
    structure_file, rule_csv, selection_csv : str or Path or None
        起始结构 / 规则表 / 导出路径 (都能在 GUI 里改)。
    dedupe : bool
        "展开时去重" 勾选框的初始值。
    on_ready : callable or None
        自检钩子: 传了就在窗口起来后把 GUI 的 ``api`` 字典交给它 (``api`` 里有全部命令
        与 state), 然后自动关窗 —— 供测试驱动全流程, 正常使用不传。
    """
    import tkinter as tk
    import tkinter.font as tkfont
    from tkinter import filedialog, messagebox

    state = {"structure": None, "view": None, "columns": None, "cells": [],
             "bad": {}, "rowframes": [], "rownums": [], "focus": 0, "job": None}
    state["rules_path"] = Path(rule_csv) if rule_csv else None
    state["out_path"] = Path(selection_csv) if selection_csv else Path("selection.csv")

    root = tk.Tk()
    root.title("原子选择表生成器 · 规则表 → 每行一个原子")
    root.geometry("1180x760+80+80")                             # 固定位置: 好截图/好找
    root.lift()
    root.attributes("-topmost", True)                           # 起窗先亮出来, 1.5 秒后交还
    root.after(1500, lambda: root.attributes("-topmost", False))

    struct_var = tk.StringVar(value="结构: (未加载)")
    rules_var = tk.StringVar(value=f"规则表: {state['rules_path'] or '(本选择器自带的 rules)'}")
    out_var = tk.StringVar(value=f"选择表(输出): {state['out_path']}")
    row_var = tk.StringVar(value="当前行: — (点任意一格即选中该行)")
    status_var = tk.StringVar(value="")
    dedupe_var = tk.BooleanVar(value=dedupe)                    # 展开时去掉重复原子 (同五元组只留第一次)
    plain_font = tkfont.nametofont("TkDefaultFont")
    bold_font = plain_font.copy()
    bold_font.configure(weight="bold")

    def sync(r):
        selector.rules[r] = [cell.get() for cell in state["cells"][r]]

    def sync_all():
        for r in range(len(state["cells"])):
            sync(r)

    def row_patterns(r):
        cols = state["columns"]
        return {f: selector.rules[r][cols[f]].strip() for f in cols}

    def ensure_visible(rowframe):
        """激活行滚进视野 (行多的时候最容易"看不见选的是哪行")。"""
        canvas.update_idletasks()
        y, h = rowframe.winfo_y(), rowframe.winfo_height()
        top = canvas.canvasy(0)
        if y < top or y + h > top + canvas.winfo_height():
            canvas.yview_moveto(max(0.0, (y - h) / max(1, inner.winfo_height())))

    def set_active(r):
        """激活行: ▶ 行号 + 蓝色带 + 输入框淡蓝; 有问题的行始终是红带/红底。"""
        r = max(0, min(r, len(state["rowframes"]) - 1)) if state["rowframes"] else 0
        state["focus"] = r
        for i, frame in enumerate(state["rowframes"]):
            active, problem = i == r, state["bad"].get(i, False)
            band = BAND_ACTIVE if active else (BAND_PROBLEM if problem else BAND_IDLE)
            frame.configure(bg=band, highlightbackground=band)
            state["rownums"][i].configure(text=f"▶ {i + 1}" if active else str(i + 1),
                                          fg=FG_ROW_ACTIVE if active else FG_ROW_IDLE,
                                          font=bold_font if active else plain_font)
            for cell in state["cells"][i]:
                cell.paint(problem=problem, active=active)
        if state["rowframes"]:
            row_var.set(f"当前行 {r + 1} 的原子序列 (按结构文件顺序):")
            ensure_visible(state["rowframes"][r])

    def revalidate():
        """重算每行有没有问题 (缓存进 state['bad']) 并按缓存重画整表。"""
        state["bad"] = ({r: bool(match_rule(row_patterns(r), state["view"])[1])
                         for r in range(len(selector.rules))}
                        if state["view"] is not None and state["columns"] is not None else {})
        sync_all()
        set_active(state["focus"])

    def expand_now():
        """表里当前的规则 → 展平结果; 没有结构/表时给空。"""
        if state["view"] is None or state["columns"] is None:
            return [], {}, 0, []
        sync_all()
        _view, hits, keep, warns, dup = selector._match(state["structure"])
        flat = [(r, a) for r, row in enumerate(hits) for a in row]
        uniq = [(r, a) for r, row in enumerate(keep) for a in row]
        return flat, warns, dup, uniq

    def refresh_status():
        flat, warns, dup, uniq = expand_now()
        keep = len(uniq) if dedupe_var.get() else len(flat)
        bits = [f"规则 {len(selector.rules)} 行"]
        if state["view"] is not None:
            bits += [f"结构 {state['view'].summary()}",
                     f"展平 {len(flat)} → 选择表 {keep} 个原子",
                     (f"去重丢掉 {dup}" if dedupe_var.get() else f"重复 {dup} (未去重)") if dup else "无重复",
                     f"{len(warns)} 行有问题" if warns else "每行都有命中"]
        else:
            bits.append("无结构")
        status_var.set(" · ".join(bits))

    def show_preview(r):
        preview.configure(state="normal")
        preview.delete("1.0", "end")
        if state["view"] is not None and state["columns"] is not None and r < len(selector.rules):
            hits, warns = match_rule(row_patterns(r), state["view"])
            out = [f"行 {r + 1}: {len(hits)} 个原子"]
            out += [f"  ⚠ {w}" for w in warns]
            if hits:
                shown = " ".join(state["view"].atom_str(a) for a in hits[:60])
                out.append("  " + shown + (f" … (+{len(hits) - 60})" if len(hits) > 60 else ""))
            preview.insert("1.0", "\n".join(out))
        preview.configure(state="disabled")

    def on_edit(r):
        sync(r)
        if state["view"] is not None and state["columns"] is not None:
            state["bad"][r] = bool(match_rule(row_patterns(r), state["view"])[1])
        set_active(r)
        if state["job"]:
            root.after_cancel(state["job"])
        state["job"] = root.after(250, lambda: (show_preview(r), refresh_status()))

    def rebuild():
        for w in inner.winfo_children():
            w.destroy()
        state["cells"], state["rowframes"], state["rownums"] = [], [], []
        col2field = {i: f for f, i in (state["columns"] or {}).items()}
        # 表头: 第 0 列是行号; 每个字段占 2 列 (模式下拉 + 输入框), 其余列 1 列
        tk.Label(inner, text="#", width=4, anchor="e").grid(row=0, column=0, sticky="e")
        at, gcol = {}, 1
        for c, h in enumerate(selector.header):
            at[c] = gcol
            span = 2 if c in col2field else 1
            tk.Label(inner, text=h, anchor="w").grid(row=0, column=gcol, columnspan=span,
                                                     sticky="w", padx=3)
            gcol += span
        for r, row in enumerate(selector.rules):
            band = tk.Frame(inner, bd=0, highlightthickness=2, bg=BAND_IDLE, highlightbackground=BAND_IDLE)
            band.grid(row=r + 1, column=0, columnspan=gcol, sticky="we", padx=1, pady=1)
            num = tk.Label(band, text=str(r + 1), width=4, anchor="e", fg=FG_ROW_IDLE, font=plain_font)
            num.grid(row=0, column=0, sticky="e", padx=2, pady=1)
            cells = []
            for c in range(len(selector.header)):
                value = row[c] if c < len(row) else ""
                cell = (FieldCell(band, col2field[c], value, lambda r=r: on_edit(r))
                        if c in col2field else PlainCell(band, value, lambda r=r: on_edit(r)))
                for i, w in enumerate(cell.widgets):
                    w.grid(row=0, column=at[c] + i, sticky="we", padx=2, pady=1)
                cells.append(cell)
            state["cells"].append(cells)
            state["rowframes"].append(band)
            state["rownums"].append(num)
        inner.update_idletasks()
        canvas.configure(scrollregion=canvas.bbox("all"))
        revalidate()

    def load_struct(path):
        if not path:
            return
        try:
            state["structure"] = StructureFile_AtomArray(input_io=path).read()
            state["view"] = _TargetView(state["structure"])
        except Exception as exc:                                # noqa: BLE001 - GUI 里报给人看
            messagebox.showerror("读结构失败", str(exc))
            return
        struct_var.set(f"结构: {path} ({state['view'].summary()})")
        revalidate()
        refresh_status()

    def load_rules(path):
        try:
            fresh = RuleCsv_AtomArraySelection(input_io=path).read()
        except (OSError, ValueError) as exc:
            messagebox.showerror("读规则表失败", str(exc))
            return
        selector.header, selector.rules = fresh.header, fresh.rules
        state["rules_path"] = Path(path)
        state["columns"] = selector.columns
        state["focus"] = 0
        rules_var.set(f"规则表: {path} ({len(selector.rules)} 行)")
        rebuild()
        refresh_status()
        show_preview(0)

    def save_rules(path=None):
        sync_all()
        path = path or filedialog.asksaveasfilename(
            defaultextension=".csv",
            initialfile=(state["rules_path"] or Path("rules.csv")).name)
        if not path:
            return
        AtomArraySelection_RuleCsv(output_io=path).write(selector)
        state["rules_path"] = Path(path)
        rules_var.set(f"规则表: {path} (已保存 {len(selector.rules)} 行规则)")

    def export(path=None):
        if state["structure"] is None:
            messagebox.showwarning("没有结构", "先打开 PDB/CIF —— 选择表要按结构展开")
            return
        flat, warns, dup, uniq = expand_now()
        keep = uniq if dedupe_var.get() else flat
        tail = (f"去重后 {len(keep)} 个原子" if dedupe_var.get()
                else f"未去重 {len(keep)} 个原子 (其中重复 {dup})")
        if warns and (path is None):                            # 0 命中的规则会凭空消失, 先问一声
            bad = "; ".join(f"行 {r + 1}: {'; '.join(w)}" for r, w in sorted(warns.items())[:6])
            if not messagebox.askokcancel(
                    "这些规则命中不到原子",
                    f"{len(warns)} 行有问题, 它们不会出现在选择表里:\n{bad}\n\n仍要导出 {tail} 吗?"):
                return
        path = path or filedialog.asksaveasfilename(
            defaultextension=".csv", initialfile=state["out_path"].name)
        if not path:
            return
        AtomArraySelection_SelectionCsv(output_io=path).write(
            selector, state["structure"], dedupe=dedupe_var.get())
        state["out_path"] = Path(path)
        out_var.set(f"选择表(输出): {path} ({len(keep)} 个原子"
                    + (f", 去重丢掉 {len(flat) - len(keep)}, 重复 {dup}" if dedupe_var.get() and dup
                       else f", 含重复 {dup}" if dup else "") + ")")
        refresh_status()

    def check_all():
        sync_all()
        if state["view"] is None or state["columns"] is None:
            return
        bad, lines = [], []
        for r in range(len(selector.rules)):
            _hits, warns = match_rule(row_patterns(r), state["view"])
            state["bad"][r] = bool(warns)
            if warns:
                bad.append(r)
                lines.append(f"行 {r + 1} [{','.join(selector.rules[r])}] ⚠ {'; '.join(warns)}")
        set_active(state["focus"])                              # 按新结果重画 (红带/红底)
        flat, _warns, dup, uniq = expand_now()
        keep = len(uniq) if dedupe_var.get() else len(flat)
        head = (f"规则 {len(selector.rules)} 行: {len(bad)} 行有问题; "
                f"展平 {len(flat)} → 选择表 {keep} 个原子"
                + (f", 去重丢掉 {dup}" if dedupe_var.get() and dup else f", 重复 {dup}" if dup else ""))
        preview.configure(state="normal")
        preview.delete("1.0", "end")
        preview.insert("1.0", "\n".join([head] + lines[:40]
                                        + ([f"… 其余 {len(bad) - 40} 行见红色高亮"] if len(bad) > 40 else [])))
        preview.configure(state="disabled")
        refresh_status()

    def add_row():
        sync_all()
        selector.add_rule()
        rebuild()
        refresh_status()

    def del_row():
        if state["focus"] < len(selector.rules):
            sync_all()
            selector.remove_rule(state["focus"])
            state["focus"] = max(0, state["focus"] - 1)
            rebuild()
            refresh_status()

    def move_row(delta):
        sync_all()
        state["focus"] = selector.move_rule(state["focus"], delta)
        rebuild()

    # ---- 部件
    top = tk.Frame(root)
    top.pack(fill="x", padx=8, pady=(6, 2))
    tk.Label(top, textvariable=struct_var, anchor="w").pack(anchor="w")
    tk.Label(top, textvariable=rules_var, anchor="w").pack(anchor="w")
    tk.Label(top, textvariable=out_var, anchor="w").pack(anchor="w")

    bar = tk.Frame(root)
    bar.pack(fill="x", padx=8, pady=(4, 0))
    for text, cmd in (("打开结构…", lambda: load_struct(filedialog.askopenfilename(
                          filetypes=[("structure", "*.pdb *.cif *.mmcif"), ("all", "*")]))),
                      ("打开规则表…", lambda: load_rules(filedialog.askopenfilename(
                          filetypes=[("csv", "*.csv"), ("all", "*")]))),
                      ("保存规则表", lambda: save_rules()),
                      ("规则表另存为…", lambda: save_rules(filedialog.asksaveasfilename(
                          defaultextension=".csv",
                          initialfile=(state["rules_path"] or Path("rules.csv")).name))),
                      ("导出选择表", lambda: export()),
                      ("选择表另存为…", lambda: export(filedialog.asksaveasfilename(
                          defaultextension=".csv", initialfile=state["out_path"].name)))):
        tk.Button(bar, text=text, command=cmd).pack(side="left", padx=2)

    bar2 = tk.Frame(root)
    bar2.pack(fill="x", padx=8, pady=2)
    for text, cmd in (("解析全部", check_all), ("添加行", add_row), ("删除行", del_row),
                      ("上移", lambda: move_row(-1)), ("下移", lambda: move_row(1))):
        tk.Button(bar2, text=text, command=cmd).pack(side="left", padx=2)
    tk.Checkbutton(bar2, text="展开时去重 (同一原子只留第一次出现)", variable=dedupe_var,
                   command=refresh_status).pack(side="left", padx=10)
    tk.Label(bar2, text="每行一条规则: 左边下拉选「怎么解释」, 空框灰字就是该模式的写法。",
             anchor="w", fg="#555555").pack(side="left", padx=6)

    mid = tk.Frame(root)
    mid.pack(fill="both", expand=True, padx=8)
    canvas = tk.Canvas(mid, highlightthickness=0)
    vsb = tk.Scrollbar(mid, orient="vertical", command=canvas.yview)
    canvas.configure(yscrollcommand=vsb.set)
    vsb.pack(side="right", fill="y")
    canvas.pack(side="left", fill="both", expand=True)
    inner = tk.Frame(canvas)
    canvas.create_window((0, 0), window=inner, anchor="nw")
    inner.bind("<Configure>", lambda _e: canvas.configure(scrollregion=canvas.bbox("all")))
    canvas.bind_all("<MouseWheel>", lambda e: canvas.yview_scroll(-e.delta, "units"))

    tk.Label(root, textvariable=row_var, anchor="w").pack(fill="x", padx=8, pady=(6, 0))
    preview = tk.Text(root, height=9, wrap="word", state="disabled", background="#f6f6f6")
    preview.pack(fill="both", padx=8, pady=(2, 2))
    tk.Label(root, textvariable=status_var, anchor="w").pack(fill="x", padx=8, pady=(0, 6))
    root.bind("<Command-s>", lambda _e: save_rules())
    root.bind("<Command-e>", lambda _e: export())

    def start():
        # 表头在构造时已经验过 (AtomArraySelection.__post_init__), 这里只取一次列号缓存
        state["columns"] = selector.columns
        if rule_csv:
            load_rules(rule_csv)
        else:
            rebuild()
        if structure_file:
            load_struct(structure_file)
        refresh_status()

    failures = []
    if on_ready is not None:                                    # 自检钩子: 驱动完自动关窗
        def _probe():
            try:
                on_ready({"root": root, "state": state, "selector": selector,
                          "struct_var": struct_var, "rules_var": rules_var, "out_var": out_var,
                          "row_var": row_var, "status_var": status_var, "dedupe_var": dedupe_var,
                          "start": start, "rebuild": rebuild, "check_all": check_all,
                          "on_edit": on_edit, "show_preview": show_preview,
                          "refresh_status": refresh_status, "save_rules": save_rules,
                          "export": export, "load_rules": load_rules, "load_struct": load_struct,
                          "add_row": add_row, "del_row": del_row, "move_row": move_row})
            except AssertionError as exc:
                failures.append(exc)
            finally:
                root.destroy()                                  # 断言失败也不留窗给人
        root.after(300, _probe)

    start()
    # 外观要等窗口建好再钉: Tk 在没有 NSWindow 时设不上 (会警告 "Failed to read
    # appearance name"), 而 title/geometry/lift/topmost 又可能把设过的外观打回 auto。
    root.update_idletasks()
    force_light_appearance(root)
    root.mainloop()
    if failures:
        raise failures[0]
    return selector
