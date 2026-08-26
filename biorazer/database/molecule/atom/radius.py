# -*- coding: utf-8 -*-
"""Van der Waals radii (Å) by element.

数据集来源与诚实性说明
────────────────────────
本模块记录元素**范德华半径 (Å)**。原 ``biorazer/database/atom.py`` 只存
6 个元素的点值 (H/C/N/O/P/S), 未记录出处 —— 这些数值与 biotite 的
``vdw_radii=\"Single\"`` 集 (Bondi 1964 派生) 一致, 但本库无法回溯到
单一论文, 故 ``source`` 标注为 ``\"consensus\"`` 并如实说明。

每条记录为 ``{mean, std, lb, up, source}``; 未查到 spread 的字段为
``np.nan``。
"""

import numpy as np

#: 元素范德华半径 (Å)。``mean`` 即半径值; 无样本 spread, 故 std/lb/up 为 nan。
ATOM_RADIUS = {
    "H": {"mean": np.float32(0), "std": np.nan, "lb": np.nan, "up": np.nan, "source": "consensus"},
    "C": {"mean": np.float32(1.70), "std": np.nan, "lb": np.nan, "up": np.nan, "source": "consensus"},
    "N": {"mean": np.float32(1.60), "std": np.nan, "lb": np.nan, "up": np.nan, "source": "consensus"},
    "O": {"mean": np.float32(1.50), "std": np.nan, "lb": np.nan, "up": np.nan, "source": "consensus"},
    "P": {"mean": np.float32(1.80), "std": np.nan, "lb": np.nan, "up": np.nan, "source": "consensus"},
    "S": {"mean": np.float32(1.80), "std": np.nan, "lb": np.nan, "up": np.nan, "source": "consensus"},
}

#: 兼容旧名 (``biorazer/database/atom.py`` 的 ``vdw_dict``/``vdw_radii``)。
vdw_dict = {elm: float(rec["mean"]) for elm, rec in ATOM_RADIUS.items()}


def vdw_radii(element):
    """Van der Waals radius (Å) of an element (vectorized over arrays)."""
    return vdw_dict[element]


vdw_radii = np.vectorize(vdw_radii)
