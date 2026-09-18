# -*- coding: utf-8 -*-
"""Rosetta 内附 / Dunbrack 官方 rotamer 库读取器。

数据集来源与诚实性说明
────────────────────────
Rosetta 随发行版附带多套 rotamer 库, 全部位于
``$ROSETTA3_DB/rotamer/`` (macOS 实测:
``/Applications/Rosetta/rosetta.source.release-408/main/database/rotamer/``)::

    bbdep02.May.sortlib                               Dunbrack 2002 合并库 (466829 行)
    bbdep02.May.sortlib-correct.12.2010               2010 年修正版
    bbdep02.May.sortlib.Dunbrack02.lib.bin            Rosetta 预编译 bin
    shapovalov/StpDwn_*-*-*/<res>.bbdep.rotamers.lib.gz   Shapovalov 2010 逐残基库
    shapovalov/Suggested.Default.txt                  官方推荐 = StpDwn_5-5-5
    ExtendedOpt1-5/<res>.bbdep.rotamers.lib.gz        另一套 (扩展模式)
    beta_nov2016/                                     beta_nov2016 打分函数用
    corrections_conway2016/                           2016 修正
    ncaa_rotlibs/ peptoid_rotlibs/ dna/               非天然/核酸

**列语义** (对照 Rosetta 源码
``core/pack/dunbrack/RotamericSingleResidueDunbrackLibraryParser.cc`` 的
``read_file()`` 注释 a–q, 非推测; N = 主链扭转数列数, 通常 2 = phi,psi)::

    a.    three-letter-code
    b,c.  phi, psi ... (共 N 列)
    d.    count
    e-h.  r1..r4            rotamer well 索引 (源码注释: "only used for warnings")
    i.    probability
    i'.   -log(Probabil)    仅 Shapovalov 文件有 (夹在 probability 与 chimeans 之间)
    j-m.  chimean1..chimean4
    n-q.  chisd1..chisd4

总列数: **15+N = Dun02** (无 -lnP), **16+N = Shapovalov** (有 -lnP)。判定逻辑
见源码 ``check_for_extra_column()``::

    if ( counter == 11 + num_mainchain_torsions_ + max_possible_chis_ - shift ) return false;
    runtime_assert( counter == 12 + num_mainchain_torsions_ + max_possible_chis_ - shift, ... )

即 15+N / 16+N 两种合法列数 (max_possible_chis_ = 4)。本模块按此自动识别,
不依赖文件名猜测。

**r1..r4 列的含义**: 源码把它们读进 ``rotwells_``, 后续用于
``determine_rotamer_well_order()`` 把 rotamer 重排到 Rosetta 的
``[0, 360)`` 约定顺序。``0`` 表示该 chi 在此 rotamer 中未定义 (例如 HIS 只有
2 个 chi, r3/r4 恒为 0)。

**实测内容** (6VY1 分析时核实):

* ``bbdep02.May.sortlib``: 466829 行 / 18 种残基 / 每行 17 列 (= 15+2)。
  HIS 12321 行 = 1369 箱 × 9 rotamer。
* ``shapovalov/StpDwn_5-5-5/his.bbdep.rotamers.lib.gz``: 49284 行 / 每行 18 列
  (= 16+2); 注释头给出 ``Number of chi angles = 2``、
  ``bins for each discrete chi angle = [3, 12]``、``Number of rotamers = 36``、
  ``phi/psi step = 10 deg``。HIS 1369 箱 × 36 rotamer。

**与 PyMOL 的关系**: PyMOL 自带的 ``sc_bb_dep.pkl`` 是 ``bbdep02.May.sortlib``
按概率降序的子集 (实测 169 个共同箱中 168 个完全一致)。
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

from ..rotamer import RotamerLibrary, RotamerRecord

__all__ = [
    "ROSETTA_ROTAMER_DIR",
    "DEFAULT_STEPDOWN",
    "N_ROTWELL_COLS",
    "N_MAINCHAIN",
    "MAX_POSSIBLE_CHIS",
    "default_bbdep02_path",
    "default_shapovalov_path",
    "read_rosetta_text",
    "read_shapovalov",
    "read_bbdep02",
]

#: Rosetta ``database/rotamer`` 根目录 (macOS 实测路径)
ROSETTA_ROTAMER_DIR = Path(
    "/Applications/Rosetta/rosetta.source.release-408/main/database/rotamer"
)

#: Rosetta 官方推荐的步降版本 (``shapovalov/Suggested.Default.txt`` 内容)
DEFAULT_STEPDOWN = "StpDwn_5-5-5"

#: 文本库中 r1..r4 的列数 (源码注释 e-h)
N_ROTWELL_COLS = 4
#: 主链扭转数列数 (phi, psi)
N_MAINCHAIN = 2
#: Rosetta 文本库的 chi 列数上限 (源码 ``max_possible_chis_``)
MAX_POSSIBLE_CHIS = 4


def default_bbdep02_path(root: Path | str = ROSETTA_ROTAMER_DIR) -> Path:
    """Dunbrack 2002 合并库 ``bbdep02.May.sortlib`` 的路径。"""
    return Path(root) / "bbdep02.May.sortlib"


def default_shapovalov_path(resname: str, stepdown: str = DEFAULT_STEPDOWN,
                            root: Path | str = ROSETTA_ROTAMER_DIR) -> Path:
    """``shapovalov/<stepdown>/<res>.bbdep.rotamers.lib.gz`` 的路径。

    Parameters
    ----------
    resname : 三字母残基名, 如 ``"HIS"`` (大小写不敏感)
    stepdown : 步降目录名, 默认 ``StpDwn_5-5-5`` (Rosetta 推荐值);
        可选 ``StpDwn_0-0-0`` / ``StpDwn_2-2-2`` / ``StpDwn_10-10-10`` /
        ``StpDwn_20-20-20`` / ``StpDwn_25-25-25``
    root : Rosetta ``database/rotamer`` 根目录
    """
    return Path(root) / "shapovalov" / stepdown / f"{resname.lower()}.bbdep.rotamers.lib.gz"


def _is_shapovalov(n_cols: int, n_mainchain: int = N_MAINCHAIN,
                   max_chis: int = MAX_POSSIBLE_CHIS) -> bool:
    """按列数判定是否 Shapovalov 格式 (多一列 ``-log(P)``)。

    对照源码 ``check_for_extra_column()``: 15+N = Dun02, 16+N = Shapovalov。
    """
    return n_cols == 12 + n_mainchain + max_chis


def _parse_header(text: str) -> dict:
    """从 ``#`` 注释头提取元信息 (chi 数、bin 数、步长等)。"""
    info: dict = {}
    for line in text.splitlines():
        if not line.startswith("#"):
            continue
        body = line.lstrip("#").strip()
        for label, key in (
            ("Number of chi angles (degrees of freedom)", "n_chi"),
            ("Number of chi angles treated as discrete", "n_chi_discrete"),
            ("Number of chi angles treated as continuous", "n_chi_continuous"),
            ("Number of rotamers for discrete chi angles", "n_rotamers"),
            ("TotalDatapointsNum", "n_datapoints"),
            ("phi step, deg", "phi_step"),
            ("psi step, deg", "psi_step"),
            ("Residue type", "resname"),
        ):
            if body.startswith(label):
                val = body[len(label):].strip()
                try:
                    info[key] = int(val)
                except ValueError:
                    try:
                        info[key] = float(val)
                    except ValueError:
                        info[key] = val
        m = re.match(r"Number of bins for each discrete chi angle\s*\[(.*)\]", body)
        if m:
            info["chi_bins"] = [int(x) for x in m.group(1).split(",")]
    return info


def read_rosetta_text(path: Path | str, resname: str | None = None,
                      n_mainchain: int = N_MAINCHAIN,
                      max_chis: int = MAX_POSSIBLE_CHIS
                      ) -> RotamerLibrary | dict[str, RotamerLibrary]:
    """读 Rosetta / Dunbrack 空白分隔文本 rotamer 库 (可含多残基)。

    适用 ``bbdep02.May.sortlib``、``shapovalov/StpDwn_*/*.lib.gz``、
    ``ExtendedOpt1-5/*.lib.gz`` 等。``.gz`` 后缀自动 gunzip。

    Parameters
    ----------
    path : 文件路径 (支持 ``.gz``)
    resname : 只保留该残基 (``None`` = 返回 ``{RES: RotamerLibrary}``)
    n_mainchain : 主链扭转数列数, 默认 2 (phi, psi)
    max_chis : 文件中的 chi 列数, 默认 4 (Rosetta 文本库固定 4 列, 未定义的为 0)

    Returns
    -------
    ``resname`` 为 ``None`` 时返回 ``{RES: RotamerLibrary}``, 否则返回单个
    :class:`RotamerLibrary`。

    Notes
    -----
    列切分严格按源码 ``read_file()`` 的 a–q 顺序, 依赖空白分隔而非固定列宽。
    Shapovalov 文件由列数自动识别 (多一列 ``-log(P)``), 不靠文件名猜测。
    ``#`` 注释头的元信息放在返回对象的 ``header`` 属性里。
    """
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", errors="replace") as fh:
        text = fh.read()

    header = _parse_header(text)
    body = [l for l in text.splitlines() if l.strip() and not l.startswith("#")]
    if not body:
        raise ValueError(f"no data lines in {path}")

    out: dict[str, RotamerLibrary] = {}
    for line in body:
        f = line.split()
        shap = _is_shapovalov(len(f), n_mainchain, max_chis)
        exp = (16 if shap else 15) + n_mainchain
        if len(f) != exp:
            raise ValueError(
                f"{path.name}: expected {exp} columns (shapovalov={shap}), got {len(f)}"
            )
        res = f[0]
        if resname and res != resname:
            continue
        i = 1
        bb = [float(x) for x in f[i:i + n_mainchain]]
        i += n_mainchain
        count = int(f[i]); i += 1
        rotwell = tuple(int(x) for x in f[i:i + N_ROTWELL_COLS]); i += N_ROTWELL_COLS
        prob = float(f[i]); i += 1
        minusln = float("nan")
        if shap:
            minusln = float(f[i]); i += 1
        chimean = tuple(float(x) for x in f[i:i + max_chis]); i += max_chis
        chistd = tuple(float(x) for x in f[i:i + max_chis]); i += max_chis

        if res not in out:
            out[res] = RotamerLibrary(
                resname=res, source=path.name, n_mainchain=n_mainchain
            )
        out[res].records.append(
            RotamerRecord(
                resname=res,
                chi=chimean,
                probability=prob,
                chi_std=chistd,
                rotwell=rotwell,
                count=count,
                minus_log_prob=minusln,
                phi=bb[0] if n_mainchain >= 1 else None,
                psi=bb[1] if n_mainchain >= 2 else None,
                source=out[res].source,
            )
        )
    for lib in out.values():
        lib.header = header
    if resname:
        if resname not in out:
            raise KeyError(f"{resname} not found in {path}")
        return out[resname]
    return out


def read_shapovalov(resname: str, stepdown: str = DEFAULT_STEPDOWN,
                    root: Path | str = ROSETTA_ROTAMER_DIR,
                    path: Path | str | None = None) -> RotamerLibrary:
    """读 Rosetta 内附的 Shapovalov 2010 骨架依赖库 (单残基一个 gz 文件)。

    Parameters
    ----------
    resname : 三字母残基名 (大小写不敏感)
    stepdown : ``StpDwn_5-5-5`` (Rosetta 推荐) / ``StpDwn_0-0-0`` 等
    root : Rosetta ``database/rotamer`` 根目录
    path : 直接指定文件 (给了就忽略上面三项)

    Returns
    -------
    :class:`RotamerLibrary`; ``header`` 含 chi 数/bin 数/步长等元信息。

    Notes
    -----
    实测 HIS: 49284 条 = 1369 个 (phi,psi) 箱 × 36 rotamer; ``chi_bins = [3, 12]``。
    """
    p = Path(path) if path else default_shapovalov_path(resname, stepdown, root)
    lib = read_rosetta_text(p, resname=resname.upper())
    assert isinstance(lib, RotamerLibrary)
    return lib


def read_bbdep02(root: Path | str = ROSETTA_ROTAMER_DIR,
                 resname: str | None = None
                 ) -> RotamerLibrary | dict[str, RotamerLibrary]:
    """读 Dunbrack 2002 合并库 ``bbdep02.May.sortlib`` (17 列 = 15+2, 无 -lnP)。

    Parameters
    ----------
    root : Rosetta ``database/rotamer`` 根目录
    resname : 只保留该残基 (``None`` = 返回全部 18 种)

    Notes
    -----
    实测: 466829 行 / 18 种残基; HIS 12321 行 = 1369 箱 × 9 rotamer。
    该库是 PyMOL ``sc_bb_dep.pkl`` 的超集。
    """
    return read_rosetta_text(default_bbdep02_path(root), resname=resname)
