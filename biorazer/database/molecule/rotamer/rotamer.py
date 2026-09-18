# -*- coding: utf-8 -*-
"""rotamer 库的公共记录类型。

``biorazer.database.molecule.rotamer`` 下有两个来源子包:

* :mod:`.rosetta` -- Rosetta 内附 / Dunbrack 官方文本库
  (``bbdep02.May.sortlib``、``shapovalov/StpDwn_*``)。
* :mod:`.pymol`   -- PyMOL 自带的 Dunbrack pickle
  (``data/chempy/sidechains/sc_bb_ind.pkl`` / ``sc_bb_dep.pkl``)。

两者都把数据归一化成这里定义的 :class:`RotamerRecord` /
:class:`RotamerLibrary`, 因此下游代码无需关心来源。

所有角度单位 **度 (degree)**。
"""

from __future__ import annotations

from dataclasses import dataclass, field

__all__ = ["RotamerRecord", "RotamerLibrary", "N_MAINCHAIN", "DEFAULT_BIN_GRID"]

#: 默认主链扭转数列数 (phi, psi)
N_MAINCHAIN = 2
#: backbone-dependent 库的 (phi, psi) 箱步长 (度); PyMOL 与 Shapovalov 均为 10
DEFAULT_BIN_GRID = 10.0


@dataclass
class RotamerRecord:
    """一条 rotamer 记录 (跨库统一字段)。

    Attributes
    ----------
    resname : 三字母残基名
    chi : chi 角度均值 (度); Rosetta 文本库固定 4 列, 未定义的 chi 为 ``0.0``
        (对应 ``rotwell == 0``)
    probability : 该 rotamer 在所属 (phi,psi) 箱内的概率
    chi_std : chi 标准差 (度); 缺失为 ``()`` 或 ``nan``
    rotwell : Rosetta 的 rotamer well 索引 (r1..r4); 非 Rosetta 源为 ``None``
    count : 该 rotamer 的观测计数 (源码列 d); 非 Rosetta 源为 ``None``
    minus_log_prob : Shapovalov 库的 ``-log(P)`` 列 (源码列 i'); 无则为 ``nan``
    phi, psi : 所属主链箱 (度); backbone-independent 源为 ``None``
    source : 数据来源标签 (文件路径或库名)
    chi_quads : chi 对应的原子名四元组 (仅 PyMOL pickle 源提供, 其余为 ``None``)
    """

    resname: str
    chi: tuple[float, ...]
    probability: float
    chi_std: tuple[float, ...] = ()
    rotwell: tuple[int, ...] | None = None
    count: int | None = None
    minus_log_prob: float = float("nan")
    phi: float | None = None
    psi: float | None = None
    source: str = ""
    chi_quads: tuple[tuple[str, str, str, str], ...] | None = None


@dataclass
class RotamerLibrary:
    """一个残基的 rotamer 库 (可含多个 (phi,psi) 箱)。

    ``records`` 保留文件顺序 (Shapovalov 按概率降序)。
    """

    resname: str
    records: list[RotamerRecord] = field(default_factory=list)
    source: str = ""
    chi_quads: tuple[tuple[str, str, str, str], ...] | None = None
    n_mainchain: int = N_MAINCHAIN
    #: 解析自 ``#`` 注释头的元信息 (chi 数/bin 数/步长等); PyMOL 源为空
    header: dict = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.records)

    def __iter__(self):
        return iter(self.records)

    @property
    def bins(self) -> list[tuple[float, float]]:
        """出现过的 (phi, psi) 箱 (去重, 保持文件顺序)。"""
        seen, out = set(), []
        for r in self.records:
            if r.phi is None:
                continue
            key = (r.phi, r.psi)
            if key not in seen:
                seen.add(key)
                out.append(key)
        return out

    def for_backbone(self, phi: float, psi: float,
                     grid: float = DEFAULT_BIN_GRID) -> list[RotamerRecord]:
        """取最接近的 (phi, psi) 箱的 rotamer, 按概率降序。

        箱以 ``grid`` (默认 10 度, 库的固有分辨率) 为步长, 输入就近取整到箱中心;
        backbone-independent 库 (``phi``/``psi`` 为 ``None``) 直接返回全部记录。
        该箱无数据时退回全局按概率降序。
        """
        if not self.records or self.records[0].phi is None:
            return sorted(self.records, key=lambda r: -r.probability)
        tp = max(-180.0, min(180.0, round(phi / grid) * grid))
        ts = max(-180.0, min(180.0, round(psi / grid) * grid))
        cand = [r for r in self.records if r.phi == tp and r.psi == ts]
        if not cand:
            return sorted(self.records, key=lambda r: -r.probability)
        return sorted(cand, key=lambda r: -r.probability)

    def top(self, n: int = 1, phi: float | None = None,
            psi: float | None = None) -> list[RotamerRecord]:
        """概率最高的前 ``n`` 条 (给定 ``phi``/``psi`` 则限定在该箱内)。"""
        recs = (self.for_backbone(phi, psi)
                if phi is not None and psi is not None
                else sorted(self.records, key=lambda r: -r.probability))
        return recs[:n]
