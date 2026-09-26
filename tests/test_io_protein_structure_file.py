"""Tests for :class:`biorazer.structure.io.protein.StructureFile_AtomArray`.

The suffix dispatch is the whole class: each case pins that the reader the
suffix selected is the one that ran (逐原子与直读一致), and that an unknown
suffix fails loudly instead of guessing a format.
"""

import numpy as np
import pytest

from biorazer.structure.io import Cif_AtomArray, Pdb_AtomArray, StructureFile_AtomArray
from biorazer.structure.io.protein import AtomArray_Cif, AtomArray_Pdb
from biorazer.structure.objects import AtomArray


def _ala() -> AtomArray:
    """One ALA residue's backbone, i.e. enough columns for a real file round trip."""
    array = AtomArray(4)
    array.atom_name = np.array(["N", "CA", "C", "O"])
    array.res_name = np.array(["ALA"] * 4)
    array.chain_id = np.array(["A"] * 4)
    array.res_id = np.array([1] * 4)
    array.ins_code = np.array([""] * 4)
    array.element = np.array(["N", "C", "C", "O"])
    array.coord = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [2.0, 1.4, 0.0], [3.0, 2.0, 0.0]])
    return array


def _same_atoms(left: AtomArray, right: AtomArray) -> bool:
    return (len(left) == len(right)
            and left.atom_name.tolist() == right.atom_name.tolist()
            and left.res_id.tolist() == right.res_id.tolist()
            and np.allclose(left.coord, right.coord))


@pytest.fixture
def pdb_path(tmp_path):
    path = tmp_path / "ala.pdb"
    AtomArray_Pdb(output_io=path).write(_ala())
    return path


@pytest.mark.parametrize("suffix,reader,writer", [
    ("pdb", Pdb_AtomArray, AtomArray_Pdb),
    ("cif", Cif_AtomArray, AtomArray_Cif),
    ("mmcif", Cif_AtomArray, AtomArray_Cif),        # .mmcif 也走 cif 分支
])
def test_read_dispatches_on_the_suffix(tmp_path, suffix, reader, writer):
    """.pdb → Pdb_AtomArray, .cif/.mmcif → Cif_AtomArray, 结果与直读一致。"""
    path = tmp_path / f"ala.{suffix}"
    writer(output_io=path).write(_ala())

    got = StructureFile_AtomArray(input_io=path).read()
    assert _same_atoms(got, reader(input_io=path).read(altloc="all"))
    assert _same_atoms(got, _ala())


def test_altloc_defaults_to_all_and_can_be_overridden(pdb_path):
    """默认 altloc="all" (选择器要它), 调用方仍可按 biotite 的语义覆盖。"""
    assert _same_atoms(StructureFile_AtomArray(input_io=pdb_path).read(),
                       Pdb_AtomArray(input_io=pdb_path).read(altloc="all"))
    # 这份文件没有 altloc, first 与 all 结果相同 —— 这里钉的是参数确实传下去了
    assert _same_atoms(StructureFile_AtomArray(input_io=pdb_path).read(altloc="first"),
                       Pdb_AtomArray(input_io=pdb_path).read(altloc="first"))


def test_unsupported_suffix_raises(pdb_path, tmp_path):
    """认不出的后缀直接报错, 报错里带上路径 —— 不猜格式。"""
    path = tmp_path / "ala.xyz"
    path.write_text(pdb_path.read_text())
    with pytest.raises(ValueError) as excinfo:
        StructureFile_AtomArray(input_io=path).read()
    assert "unsupported structure format" in str(excinfo.value)
    assert "ala.xyz" in str(excinfo.value)
