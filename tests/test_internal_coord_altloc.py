"""Tests for the altLoc label that :class:`InternalCoordAtom` carries.

An ``InternalCoord`` is built from an ``AtomArray`` (and written back to one);
the label only exists on an array read with ``altloc="all"``, so both
directions have to agree on one spelling for "no alternate conformation" --
otherwise the round trip silently loses the annotation the selector filters on.
"""

import numpy as np

from biorazer.structure.bridge import AtomArray_InternalCoord, InternalCoord_AtomArray
from biorazer.structure.io import StructureFile_AtomArray
from biorazer.structure.objects.internal_coords import NULL_ALT

from .test_mutation import _ideal_chain


class TestInternalCoordAltloc:
    """``InternalCoordAtom.altloc_id``: 建 IC 时带上, 写回时恢复。"""

    def test_records_carry_the_label_from_the_array(self, altloc_pdb):
        """AtomArray → IC: 每条 record 的 altloc_id 与源数组逐原子一致 (空列归一成空串)。"""
        arr = StructureFile_AtomArray(input_io=altloc_pdb).read()
        ic = AtomArray_InternalCoord(input_io=arr).convert()

        assert len(ic.atoms) == len(arr) == 13
        assert ic.altloc_id.tolist() == [
            "" if str(x) in NULL_ALT else str(x) for x in arr.altloc_id]
        assert ic.altloc_id.tolist().count("A") == 1        # 第 2 号残基的 CA(A)
        assert ic.altloc_id.tolist().count("B") == 1        # 与它同名的 CA(B)

    def test_read_path_without_the_category_gives_empty_labels(self):
        """没按 altloc="all" 读的数组没有这个 category —— record 一律空, 不报错。"""
        ic = AtomArray_InternalCoord(_ideal_chain(4, "A")).convert()
        assert ic.altloc_id.tolist() == [""] * len(ic.atoms)

    def test_write_path_restores_the_annotation(self):
        """IC → AtomArray: 只要有一条 record 带标签就建 altloc_id 标注, 全空则不建。"""
        aa = _ideal_chain(4, "A")
        aa.set_annotation("altloc_id", np.array(["A"] + [""] * (len(aa) - 1)))
        back = InternalCoord_AtomArray(
            input_io=AtomArray_InternalCoord(input_io=aa).convert()).convert()

        assert "altloc_id" in back.get_annotation_categories()
        assert back.altloc_id[0] == "A" and set(back.altloc_id[1:]) == {""}

        plain = InternalCoord_AtomArray(
            input_io=AtomArray_InternalCoord(_ideal_chain(4, "A")).convert()).convert()
        assert "altloc_id" not in plain.get_annotation_categories()

    def test_mutation_keeps_the_backbone_label(self):
        """改侧链时骨架 record 会被**重建** —— altloc 得跟着走, 否则静默丢标注。"""
        from biorazer.structure.manipulation.internal_coord.mutation import mutate

        aa = _ideal_chain(4, "A")
        aa.set_annotation("altloc_id", np.array(["", "A"] + [""] * (len(aa) - 2)))
        ic = AtomArray_InternalCoord(input_io=aa).convert()
        labelled = ic.atoms[1]                     # 第 1 个残基的 CA
        assert (labelled.name, labelled.altloc_id) == ("CA", "A")

        out = mutate(ic, ["A1V"], chain="A")        # ALA -> VAL: 侧链重建
        assert len(out.atoms) > len(ic.atoms), "VAL 的侧链原子比理想骨架多"
        ca = [a for a in out.atoms if a.res_id == 1 and a.name == "CA"]
        assert len(ca) == 1 and ca[0].altloc_id == "A"
        assert ca[0] is not labelled, "骨架 record 是新对象 (输入那份不被就地改写)"
        # 新加的侧链原子来自模板: 不带 altloc
        assert all(a.altloc_id == "" for a in out.atoms if a.name in ("CB", "CG1", "CG2"))
