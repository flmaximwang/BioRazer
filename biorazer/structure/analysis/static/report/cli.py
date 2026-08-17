"""biorazer structure report 子命令门面: report-interface-dsasa (后续 report 子命令也注册于此)。

各子命令的 parser/runner 位于同目录的 report 模块 (contact.py), 本模块负责
组合注册, 并向上层模块 (biorazer.cli) 暴露 register_subcommand 与 report API。
"""
from .contact import _add_interface_dsasa_parser, report_interface_dSASA

__all__ = [
    "register_subcommand",
    "report_interface_dSASA",
]


def register_subcommand(sub) -> None:
    """在 argparse subparsers 上注册 structure report 相关子命令"""
    _add_interface_dsasa_parser(sub)
