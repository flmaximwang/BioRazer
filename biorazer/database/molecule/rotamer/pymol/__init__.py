# -*- coding: utf-8 -*-
"""PyMOL 自带的 Dunbrack rotamer 库读取器。

见 :mod:`.pymol` 的模块 docstring 了解库文件位置、``data/chempy`` 与代码模块
``site-packages/chempy`` 的命名陷阱、以及与 Dunbrack 2002 的实测关系。
"""

from .pymol import (  # noqa: F401
    PYMOL_ROTAMER_DIR,
    PYMOL_SC_BB_IND,
    PYMOL_SC_BB_DEP,
    PYMOL_SC_LIBRARY,
    read_pymol_ind,
    read_pymol_dep,
    read_pymol_library,
)

__all__ = [
    "PYMOL_ROTAMER_DIR",
    "PYMOL_SC_BB_IND",
    "PYMOL_SC_BB_DEP",
    "PYMOL_SC_LIBRARY",
    "read_pymol_ind",
    "read_pymol_dep",
    "read_pymol_library",
]
