"""Converters for an ``AtomArraySelection``: its two csv tables, and the per-atom
selection it makes on an ``AtomArray``.

Every class here is a biorazer ``A_B`` :class:`~biorazer.io.Converter` -- ``A`` is the
source, ``B`` the target.

**The two csv tables** (one selection syntax, two uses; see
:mod:`biorazer.structure.objects.selector`):

* :class:`RuleCsv_AtomArraySelection` / :class:`AtomArraySelection_RuleCsv` -- the
  **rule table**: one selection expression per row, each cell a pattern.  The
  round trip is lossless, ``0``-hit rows included: the rule table is the source of
  truth and only the expansion drops rows.
* :class:`SelectionCsv_AtomArraySelection` / :class:`AtomArraySelection_SelectionCsv`
  -- the **selection table**: exactly one atom per row, concrete values.  The
  reader escapes every field cell as a literal, so those concrete values (``A,B``,
  ``CA-CB``, ``re:x``) are not misread as syntax and a selection table can be
  opened as a rule table.  The writer needs the structure the table is expanded
  against -- it is the only thing the expansion can depend on.

These four are file-backed (``read()`` / ``write()``) because a csv table is a
file: the rule/selection table has no in-memory class of its own, unlike the
objects :mod:`biorazer.structure.bridge` normally bridges.  With that one
exception the module follows the package rule -- in-memory conversion, no files.

**The per-atom selection** -- :class:`AtomArraySelection_AtomArrayMask` and
:class:`AtomArraySelection_AtomArrayIndices` implement the parent's
:meth:`~biorazer.io.Converter.convert` (one transform, no ``read()`` /
``write()``): the selection is held in ``input_io`` and ``convert()`` returns the
selection over an array passed per call.  The array is a parameter because only
an ``AtomArray`` says which atoms exist -- the shape the ``quads`` / ``anchor``
parameters of :class:`~biorazer.structure.bridge.atom_array.AtomArray_InternalCoord`
already use.

The two targets differ in ordering, which is why both exist:

* a **mask** is aligned with the array, so the rule order is invisible in it, and
* **indices** are flat (rule order first, then atom order within a rule), an order
  a mask cannot express.  ``dedupe=False`` keeps atoms matched by more than one
  rule once per match.
"""

from __future__ import annotations

import csv

import numpy as np

from biorazer.io import Converter
from biorazer.structure.objects import AtomArray, AtomArraySelection
from biorazer.structure.objects.selector import FIELDS, encode, find_columns


# ---------------------------------------------------------------- csv 文件

def _read_csv(path) -> tuple[list[list[str]], list[str]]:
    """``path`` → (数据行, 表头); 空文件当"没有数据行, 表头 = ``FIELDS``"。"""
    with open(path, newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.reader(fh))
    if not rows:
        return [], list(FIELDS)
    return rows[1:], [h.strip() for h in rows[0]]


def _write_csv(path, header, rows) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


class RuleCsv_AtomArraySelection(Converter):
    """Reads a **rule table** into an :class:`AtomArraySelection`.

    Every cell is a pattern, read as written (``*`` wildcard, ``,`` list,
    ``1-10`` range, ``re:`` regex, else a literal).  This is the file the
    tkinter editor edits, so this direction plus
    :class:`AtomArraySelection_RuleCsv` is the lossless round trip.

    Parameters
    ----------
    input_io : str or Path
        Rule table csv; the field columns are found by header name
        (:data:`~biorazer.structure.objects.selector.FIELD_ALIASES`), extra
        columns ride along as per-row notes.
    """

    def read(self) -> AtomArraySelection:
        """Read ``self.input_io``.

        Returns
        -------
        AtomArraySelection
            One rule per data row; missing field columns raise ``ValueError``.
        """
        rows, header = _read_csv(self.input_io)
        return AtomArraySelection(rules=rows, header=header)


class SelectionCsv_AtomArraySelection(Converter):
    """Reads a **selection table** (exactly one atom per row) into a selection.

    Each field cell is taken as a **literal value** and escaped through
    :func:`~biorazer.structure.objects.selector.encode` before it becomes a
    pattern -- otherwise the concrete values a selection table holds would be
    misread as syntax (``A,B`` a list, ``CA-CB`` a range, ``re:x`` a regex).
    That is what lets a selection table be reopened as a rule table: every
    concrete value then matches only itself.

    Parameters
    ----------
    input_io : str or Path
        Selection table csv.
    """

    def read(self) -> AtomArraySelection:
        """Read ``self.input_io``.

        Returns
        -------
        AtomArraySelection
            One single-atom rule per data row.
        """
        rows, header = _read_csv(self.input_io)
        field_of = {i: f for f, i in find_columns(header).items()}
        return AtomArraySelection(
            rules=[[encode(field_of[i], "字面值", cell) if i in field_of else cell
                    for i, cell in enumerate(row)] for row in rows],
            header=header)


class AtomArraySelection_RuleCsv(Converter):
    """Writes a selection back to its **rule table**, row for row.

    Parameters
    ----------
    output_io : str or Path
        Where the rule table goes.
    """

    def write(self, tmp: AtomArraySelection) -> None:
        """Write ``tmp`` to ``self.output_io``.

        Parameters
        ----------
        tmp : AtomArraySelection
            The selection whose ``header`` / ``rules`` are written verbatim --
            including the rules that match nothing.
        """
        _write_csv(self.output_io, tmp.header, tmp.rules)


class AtomArraySelection_SelectionCsv(Converter):
    """Expands a selection against a structure, one atom per row.

    Parameters
    ----------
    output_io : str or Path
        Where the selection table goes.
    """

    def write(self, tmp: AtomArraySelection, structure, dedupe: bool = True) -> None:
        """Write the expanded ``tmp`` to ``self.output_io``.

        Parameters
        ----------
        tmp : AtomArraySelection
            The rules to expand; its ``header`` is written as the table header
            and its extra columns are copied onto every row the rule expands to.
        structure : AtomArray or InternalCoord
            The target the expansion is made against -- the only thing the
            expansion can depend on, so it is an argument here rather than
            something the selection could hold.
        dedupe : bool
            Keep an atom matched by more than one rule only at its first match
            (the default), or once per matching rule.

        Raises
        ------
        ValueError
            ``structure`` is ``None``.
        """
        if structure is None:
            raise ValueError("a selection table is the rule table expanded against a "
                             "structure: pass one, there is nothing to expand otherwise")
        view, hits, keep, _warns, _dup = tmp._match(structure)
        chosen = keep if dedupe else hits
        _write_csv(self.output_io, tmp.header,
                   tmp.expanded_rows(view, [(r, a) for r, row in enumerate(chosen)
                                            for a in row]))


class AtomArraySelection_AtomArrayMask(Converter):
    """Builds a boolean mask from an :class:`AtomArraySelection`.

    The mask is aligned with the input array (``True`` at the atoms any rule
    matches), so it can be used directly as ``atom_array[mask]`` and composed with
    the helpers in :mod:`biorazer.structure.selection.mask`.
    """

    def convert(self, atom_array: AtomArray) -> np.ndarray:
        """Select on ``atom_array`` with ``self.input_io``.

        Parameters
        ----------
        atom_array : AtomArray
            The array the patterns are matched against, and the array the
            returned mask is aligned with.

        Returns
        -------
        numpy.ndarray
            1D boolean mask of shape ``atom_array.shape``, ``True`` at every atom
            matched by at least one rule.
        """
        return self.input_io.mask(atom_array)


class AtomArraySelection_AtomArrayIndices(Converter):
    """Builds flat atom indices from an :class:`AtomArraySelection`.

    The order is the one a mask cannot express: rule by rule, and within a rule
    the atom order of the array.  It is the order
    :class:`AtomArraySelection_SelectionCsv` writes into a selection table, so
    indices and the exported table agree row by row.
    """

    def convert(self, atom_array: AtomArray, dedupe: bool = True) -> np.ndarray:
        """Select on ``atom_array`` with ``self.input_io``.

        Parameters
        ----------
        atom_array : AtomArray
            The array the patterns are matched against.
        dedupe : bool
            Keep an atom matched by more than one rule only at its first match
            (the default), or once per matching rule.

        Returns
        -------
        numpy.ndarray
            1D ``int`` array of atom indices.
        """
        return self.input_io.indices(atom_array, dedupe=dedupe)
