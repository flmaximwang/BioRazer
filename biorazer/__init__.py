"""BioRazer: a platform for analyzing various biological information.

这是一个**常规包** (有本文件), 而不是原来的 PEP 420 命名空间包。差别是实测过的:
命名空间包会把 ``sys.path`` 上**每一个**同名目录都并进 ``__path__``, 于是第二份
checkout (worktree、副本、陈旧的 editable 路径) 可以**悄悄只供应一部分模块**;
而且扫描时一个"带 ``__init__.py``"的同名目录会**整个顶掉**先出现的命名空间
portion —— 发现器遇到没有 ``__init__.py`` 的目录只记成 portion 继续往后找, 直到
碰上第一个常规包才停。有了本文件, 这个名字只属于一个目录, ``biorazer.__file__``
也会指出来是哪一份。
"""

from importlib.metadata import PackageNotFoundError, version

try:  # 版本只有一个来源: pyproject.toml (经安装后的 metadata 暴露)
    __version__ = version("biorazer")
except PackageNotFoundError:  # 直接在源码树上 import, 没装过
    __version__ = "0.0.0.dev0"
