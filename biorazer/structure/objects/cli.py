"""``biorazer select`` 子命令门面: 选择器 (规则表 → 选择表) 的 CLI 注册处。

parser/runner 在同目录的 :mod:`~biorazer.structure.objects.selector` (单一 argparse 来源),
本模块只负责组合注册, 并向上层模块 (``biorazer.cli``) 暴露 :func:`register_subcommand`。
"""

from .selector import _add_selector_parser

__all__ = ["register_subcommand"]


def register_subcommand(sub) -> None:
    """在 argparse subparsers 上注册选择器子命令"""
    _add_selector_parser(sub)
