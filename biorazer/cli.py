"""biorazer CLI — 顶层多子命令入口"""
import argparse
import sys

from biorazer.access import cli as access_cli
from biorazer.access.server.colabfold_msa import cli as colabfold_cli
from biorazer.sequence.analysis.alignment.plot import cli as plot_cli
from biorazer.structure.analysis.static.report import cli as report_cli
from biorazer.structure.objects import cli as selector_cli


def _force_utf8_console() -> None:
    """让中文帮助/错误信息在非 UTF-8 控制台上也能写出去。

    Windows 上 stdout 被重定向或捕获时用的是 locale 代码页 (cp1252/cp437/...),
    中文帮助会直接抛 UnicodeEncodeError 把 CLI 弄崩 —— 实测 CI 的 windows 腿
    `python -m biorazer.cli --help` returncode=1, 本机用 PYTHONIOENCODING=cp1252
    可复现。交互式控制台走 UTF-16 API, encoding 本来就是 utf-8, 不受影响;
    这里只兜底被重定向/被捕获的流。
    """
    for stream in (sys.stdout, sys.stderr):
        reconfigure = getattr(stream, "reconfigure", None)
        if reconfigure is None:
            continue
        encoding = (stream.encoding or "").lower().replace("-", "")
        if encoding and encoding != "utf8":
            reconfigure(encoding="utf-8", errors="replace")


def main():
    _force_utf8_console()
    parser = argparse.ArgumentParser(
        prog="biorazer",
        description="BioRazer 生物信息分析平台",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", help="子命令")

    # 注册子命令
    access_cli.register_subcommand(sub)
    colabfold_cli.register_subcommand(sub)
    plot_cli.register_subcommand(sub)
    report_cli.register_subcommand(sub)
    selector_cli.register_subcommand(sub)

    args = parser.parse_args()
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
