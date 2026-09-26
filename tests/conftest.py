"""pytest fixtures for the biorazer test suite."""

import tempfile
import os
import numpy as np
import pytest


@pytest.fixture
def tmp_msa_dir():
    """提供临时目录用于 MSA 输出测试，自动清理。"""
    with tempfile.TemporaryDirectory(prefix="biorazer_test_msa_") as d:
        yield d


@pytest.fixture
def sample_a3m_files():
    """创建两个临时的 A3M 文件用于 merge_a3m 测试。"""
    files = []
    contents = [
        (">seq1\nMTSENLYFQG\n", "a3m_1.a3m"),
        (">seq2\nWPKL\n", "a3m_2.a3m"),
    ]
    for content, name in contents:
        f = tempfile.NamedTemporaryFile(
            mode="w", suffix=f"_{name}", delete=False
        )
        f.write(content)
        f.close()
        files.append(f.name)
    yield files
    for p in files:
        os.unlink(p)


@pytest.fixture
def altloc_pdb(tmp_path):
    """一个带 altloc 的真实结构文件: 理想主链 → PDB (第 17 列手写构象字母)。

    第 2 号残基的 ``CA`` 有两个构象 (A/B), 其余原子的 altLoc 列留空 —— 即
    ``StructureFile_AtomArray`` 读回来带 ``altloc_id`` 标注的那种文件。写入口本身不输出
    ``altloc_id`` (实测 biotite 1.6), 所以构象字母直接写进文件的第 17 列。
    """
    from biorazer.structure.io.protein import AtomArray_Pdb
    from tests.test_mutation import _ideal_chain

    arr = _ideal_chain(3, "A")
    dup = int(np.where((arr.res_id == 2) & (arr.atom_name == "CA"))[0][0])
    arr = arr[np.concatenate([np.arange(len(arr)), [dup]])]
    path = tmp_path / "altloc.pdb"
    AtomArray_Pdb(input_io=None, output_io=str(path)).write(arr)

    alt = [" "] * len(arr)
    alt[dup], alt[-1] = "A", "B"
    out, seen = [], 0
    for line in path.read_text().splitlines():
        if line.startswith("ATOM"):
            line = line[:16] + alt[seen] + line[17:]
            seen += 1
        out.append(line)
    assert seen == len(arr), (seen, len(arr))               # 原子行顺序 = 数组顺序
    path.write_text("\n".join(out) + "\n")
    return path
