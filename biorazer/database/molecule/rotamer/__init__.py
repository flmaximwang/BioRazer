# -*- coding: utf-8 -*-
"""External side-chain rotamer library readers, split by source.

``biorazer.database.molecule.rotamer`` 提供**外部** rotamer 库的解析器。
与 :mod:`biorazer.database.molecule.bond.dihedral.protein.by_residue` 的关系:

* ``by_residue`` 内嵌的是**小体积**的 chi 定义与 canonical/命名 rotamer 框架
  (Rosetta ICOOR 理想值 + Dunbrack bin 中心), 够用但**不含**骨架依赖数值表。
* 本包读的是**外部大库**: PyMOL 自带的 Dunbrack pickle 与 Rosetta 内附的
  Dunbrack 2002 / Shapovalov 2010 文本库 (逐 phi/psi 箱的均值/方差/概率)。
  这些库体积大、许可各异, 故**不 vendored**, 按本机安装路径读取。

子包按来源划分
--------------

:mod:`.rosetta`
    Rosetta 内附 / Dunbrack 官方文本库。
    入口: :func:`~.rosetta.read_rosetta_text` (通用)、
    :func:`~.rosetta.read_shapovalov` (Shapovalov 2010 单残基 gz)、
    :func:`~.rosetta.read_bbdep02` (Dunbrack 2002 合并库)。
    列语义对照 Rosetta 源码 ``core/pack/dunbrack/
    RotamericSingleResidueDunbrackLibraryParser.cc`` 的 ``read_file()``
    注释 a–q 实测核实。

:mod:`.pymol`
    PyMOL 自带的 Dunbrack pickle。
    入口: :func:`~.pymol.read_pymol_ind` (backbone-independent)、
    :func:`~.pymol.read_pymol_dep` (backbone-dependent)。

公共记录类型 (跨来源统一)
-------------------------

:class:`RotamerRecord` / :class:`RotamerLibrary` 定义在 :mod:`.rotamer`,
两个子包的 reader 都返回它, 下游无需关心来源。所有角度单位 **度 (degree)**。

典型用法::

    from biorazer.database.molecule.rotamer import pymol, rosetta

    # PyMOL 自带库: 按该位点的 phi/psi 取概率最高的 rotamer
    dep = pymol.read_pymol_dep()                 # {RES: RotamerLibrary}
    top = dep["HIS"].top(1, phi=-60, psi=-40)
    print(top[0].chi)                            # (chi1, chi2) 度数

    # Rosetta 内附 Shapovalov 2010 (默认 StpDwn_5-5-5)
    his = rosetta.read_shapovalov("HIS")
    print(len(his), his.bins[:3])                # 49284 条, 1369 个 (phi,psi) 箱
    print(his.header["chi_bins"])                # [3, 12]

    # Dunbrack 2002 合并库
    bbdep = rosetta.read_bbdep02(resname="HIS")
"""

from .rotamer import (  # noqa: F401
    RotamerRecord,
    RotamerLibrary,
    N_MAINCHAIN,
    DEFAULT_BIN_GRID,
)
from . import rosetta  # noqa: F401
from . import pymol  # noqa: F401

__all__ = [
    # 公共记录类型
    "RotamerRecord",
    "RotamerLibrary",
    "N_MAINCHAIN",
    "DEFAULT_BIN_GRID",
    # 来源子包
    "rosetta",
    "pymol",
]
