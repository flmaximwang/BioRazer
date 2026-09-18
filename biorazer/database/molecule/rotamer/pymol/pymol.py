# -*- coding: utf-8 -*-
"""PyMOL 自带的 Dunbrack rotamer 库读取器。

数据集来源与诚实性说明
────────────────────────
PyMOL **自带**一份 Dunbrack rotamer 库, 无需另行下载。它随 PyMOL 一起安装,
位于 ``$PYMOL_DATA/chempy/sidechains/`` (macOS 打包版为
``/Applications/PyMOL.app/Contents/share/pymol/data/chempy/sidechains/``)::

    sc_bb_ind.pkl    19 KB   backbone-independent (每个残基一套离散 rotamer)
    sc_bb_dep.pkl   1.4 MB   backbone-dependent (按 phi/psi 10 deg 分箱)
    sc_library.pkl   28 KB   更早的简表
    sc_bb_ind.py / sc_bb_dep.py / generate1.py / generate2.py   生成脚本

**命名陷阱**: 这个目录叫 ``data/chempy/``, 与 PyMOL 的 python 代码模块
``site-packages/chempy/`` (分子数据模型, 36 个 .py, 无任何数据库) **同名但
完全无关**。rotamer 库只在 ``data/chempy/sidechains/`` 下, 不在代码模块里。

消费方是 ``pymol/wizard/mutagenesis.py``, 它读::

    os.environ['PYMOL_DATA'] + "/chempy/sidechains/sc_bb_ind.pkl"
    os.environ['PYMOL_DATA'] + "/chempy/sidechains/sc_bb_dep.pkl"

**实测内容** (HIS, 6VY1 分析时核实):

* ``sc_bb_ind``: 9 个 rotamer, 带 ``FREQ`` 与 chi 四元组 key, 概率和为 1。
* ``sc_bb_dep``: 169 个 (phi,psi) 箱, 合计 1141 条 HIS rotamer。

**与 Dunbrack 2002 的关系** (实测核对): ``sc_bb_dep`` 是
:func:`biorazer.database.molecule.rotamer.rosetta.read_rosetta_text` 读出的
``bbdep02.May.sortlib`` 按概率降序的**子集** —— 同一 (phi,psi) 箱内, PyMOL
保留前若干条高概率 rotamer, 低概率条目被生成脚本丢弃。实测 169 个共同箱中
168 个完全一致, 剩余 1 个也只是 PyMOL 少 3 条低概率条目 (PyMOL 的 6 条在
Dun02 中全部能找到且概率一致)。

**与 ``by_residue.SIDECHAIN_CHI`` 的关系**: PyMOL pickle 的 key 就是 chi 原子
名四元组, 与 :data:`...by_residue.SIDECHAIN_CHI` 完全一致 (同源 Dunbrack),
可直接对接。

.. warning::
   **key 顺序陷阱**: PyMOL pickle 的 dict key **顺序不等于 chi1..chiN 顺序**。
   实测 7 个残基被打乱 —— 例如 GLU 的 key 顺序是 ``(chi2, chi3, chi1)``、
   ILE 的 chi2 用了合并原子名 ``CD1+CD``。若按 key 顺序直接取 chi, 会把 chi
   值挂到错误的 chi 编号上。本模块的 reader 统一按
   :data:`...by_residue.SIDECHAIN_CHI` 的规范顺序重排 (见 :func:`_order_chi`),
   并把 ``CD1+CD`` 规范化成 ``CD1``。返回的 ``chi_quads`` 始终是 chi1..chiN 顺序。
"""

from __future__ import annotations

import pickle
from pathlib import Path

from ..rotamer import RotamerLibrary, RotamerRecord

try:  # 规范 chi 顺序 (chi1..chiN); 见下方 _order_chi 说明
    from ...bond.dihedral.protein.by_residue import SIDECHAIN_CHI as _SIDECHAIN_CHI
except ImportError:  # pragma: no cover - 仅在包被拆散时发生
    _SIDECHAIN_CHI = {}

__all__ = [
    "PYMOL_ROTAMER_DIR",
    "PYMOL_SC_BB_IND",
    "PYMOL_SC_BB_DEP",
    "PYMOL_SC_LIBRARY",
    "read_pymol_ind",
    "read_pymol_dep",
    "read_pymol_library",
]

#: PyMOL rotamer 库所在目录 (macOS 打包版默认路径)
PYMOL_ROTAMER_DIR = Path(
    "/Applications/PyMOL.app/Contents/share/pymol/data/chempy/sidechains"
)
PYMOL_SC_BB_IND = PYMOL_ROTAMER_DIR / "sc_bb_ind.pkl"
PYMOL_SC_BB_DEP = PYMOL_ROTAMER_DIR / "sc_bb_dep.pkl"
PYMOL_SC_LIBRARY = PYMOL_ROTAMER_DIR / "sc_library.pkl"

#: PyMOL 生成脚本对 ILE 的 chi2 用了合并名 ``CD1+CD``, 需映射到规范原子名 ``CD1``
_ATOM_NAME_FIX = {"CD1+CD": "CD1"}


def _normalise_quad(quad: tuple) -> tuple[str, str, str, str]:
    return tuple(_ATOM_NAME_FIX.get(a, a) for a in quad)  # type: ignore[return-value]


def _chi_quads(rec: dict) -> tuple[tuple[str, str, str, str], ...]:
    """从一条 PyMOL rotamer 记录取出 chi 原子名四元组。

    .. warning::
       PyMOL pickle 的 **dict key 顺序不等于 chi1..chiN 顺序**! 实测 7 个残基
       (ARG/GLN/GLU/ILE/LYS/MET/PRO) 的 key 顺序被打乱, 例如 GLU 的 key 顺序
       是 (chi2, chi3, chi1)、ILE 的 chi2 用了合并名 ``CD1+CD``。
       因此这里只做原子名规范化, **顺序**由 :func:`_order_chi` 按
       :data:`...by_residue.SIDECHAIN_CHI` 的规范定义重排。
    """
    return tuple(_normalise_quad(k) for k in rec if k != "FREQ")


def _order_chi(rec: dict, resname: str) -> tuple[tuple[float, ...], tuple[tuple[str, str, str, str], ...]]:
    """按规范 chi1..chiN 顺序取出 chi 角度, 返回 ``(chi, quads)``。

    规范顺序取自 :data:`...by_residue.SIDECHAIN_CHI` (同源 Dunbrack)。找不到
    该残基的规范定义时退回文件内 key 顺序, 并如实返回该顺序作为 quads。
    """
    canon = _SIDECHAIN_CHI.get(resname)
    if canon:
        quads = tuple(tuple(q) for q in canon)
        vals = []
        for q in quads:
            key = next((k for k in rec if k != "FREQ" and _normalise_quad(k) == q), None)
            vals.append(float(rec[key]) if key is not None else float("nan"))
        return tuple(vals), quads  # type: ignore[return-value]
    quads = _chi_quads(rec)
    return tuple(float(rec[q]) for q in quads), quads


def read_pymol_ind(path: Path | str = PYMOL_SC_BB_IND) -> dict[str, RotamerLibrary]:
    """读 PyMOL ``sc_bb_ind.pkl`` (backbone-independent)。

    Parameters
    ----------
    path : pickle 路径, 默认 PyMOL 安装目录下的 ``sc_bb_ind.pkl``

    Returns
    -------
    ``{RES: RotamerLibrary}``; 每条 record 的 ``chi`` 按文件中 key 顺序
    (即 chi1, chi2, ...), ``chi_quads`` 给出对应原子名四元组, ``phi``/``psi``
    为 ``None`` (该库不区分骨架)。

    Notes
    -----
    实测: 18 种残基; HIS 9 个 rotamer, 概率和 = 1。
    """
    raw = pickle.load(open(path, "rb"))
    out: dict[str, RotamerLibrary] = {}
    for res, recs in raw.items():
        lib = RotamerLibrary(resname=res, source=f"pymol_sc_bb_ind:{path}")
        for rec in recs:
            chi, quads = _order_chi(rec, res)
            if lib.chi_quads is None:
                lib.chi_quads = quads
            lib.records.append(
                RotamerRecord(
                    resname=res,
                    chi=chi,
                    probability=float(rec["FREQ"]),
                    phi=None,
                    psi=None,
                    source=lib.source,
                    chi_quads=quads,
                )
            )
        out[res] = lib
    return out


def read_pymol_dep(path: Path | str = PYMOL_SC_BB_DEP) -> dict[str, RotamerLibrary]:
    """读 PyMOL ``sc_bb_dep.pkl`` (backbone-dependent)。

    Parameters
    ----------
    path : pickle 路径, 默认 PyMOL 安装目录下的 ``sc_bb_dep.pkl``

    Returns
    -------
    ``{RES: RotamerLibrary}``, 每条 record 带 ``phi``/``psi`` (10 deg 箱)。

    Notes
    -----
    实测: 3569 条 (RES, phi, psi) 箱记录; HIS 1141 条 / 169 箱。
    该库是 Dunbrack 2002 按概率降序的子集, 见模块 docstring。
    """
    raw = pickle.load(open(path, "rb"))
    out: dict[str, RotamerLibrary] = {}
    for key, recs in raw.items():
        res, phi, psi = key[0], float(key[1]), float(key[2])
        if res not in out:
            out[res] = RotamerLibrary(
                resname=res, source=f"pymol_sc_bb_dep:{path}"
            )
        lib = out[res]
        for rec in recs:
            chi, quads = _order_chi(rec, res)
            if lib.chi_quads is None:
                lib.chi_quads = quads
            lib.records.append(
                RotamerRecord(
                    resname=res,
                    chi=chi,
                    probability=float(rec["FREQ"]),
                    phi=phi,
                    psi=psi,
                    source=lib.source,
                    chi_quads=quads,
                )
            )
    return out


def read_pymol_library(path: Path | str = PYMOL_SC_LIBRARY) -> dict[str, RotamerLibrary]:
    """读 PyMOL ``sc_library.pkl`` (更早的简表, 每残基一套 rotamer)。

    结构与 ``sc_bb_ind`` 相同 (``{RES: [{FREQ, quad: deg, ...}]}``),
    但条目更少、不含完整侧链二面角, 仅作兼容用途。
    """
    raw = pickle.load(open(path, "rb"))
    out: dict[str, RotamerLibrary] = {}
    for res, recs in raw.items():
        lib = RotamerLibrary(resname=res, source=f"pymol_sc_library:{path}")
        for rec in recs:
            if not isinstance(rec, dict) or "FREQ" not in rec:
                continue
            chi, quads = _order_chi(rec, res)
            if lib.chi_quads is None:
                lib.chi_quads = quads
            lib.records.append(
                RotamerRecord(
                    resname=res,
                    chi=chi,
                    probability=float(rec["FREQ"]),
                    source=lib.source,
                    chi_quads=quads,
                )
            )
        out[res] = lib
    return out
