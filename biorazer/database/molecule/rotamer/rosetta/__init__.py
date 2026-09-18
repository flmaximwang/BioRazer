# -*- coding: utf-8 -*-
"""Rosetta 内附 / Dunbrack 官方 rotamer 库读取器。

见 :mod:`.rosetta` 的模块 docstring 了解库文件清单、列语义 (对照 Rosetta
源码) 与实测内容。
"""

from .rosetta import (  # noqa: F401
    ROSETTA_ROTAMER_DIR,
    DEFAULT_STEPDOWN,
    N_ROTWELL_COLS,
    N_MAINCHAIN,
    MAX_POSSIBLE_CHIS,
    default_bbdep02_path,
    default_shapovalov_path,
    read_rosetta_text,
    read_shapovalov,
    read_bbdep02,
)

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
