# -*- coding: utf-8 -*-
"""Generate the literal ``IC_<Res>_<SS>[_<rotamer>]`` named templates.

Rewrites **each per-residue module** to a header that imports
:func:`build_template` and a literal block of named-template constants built
directly from it, e.g. ``gly.py``::

    from biorazer.database.molecule.icoor.protein.template._builder import build_template

    IC_Gly_HelixAlpha = build_template("GLY", "alpha-helix", "canonical")
    ...
    __all__ = ["IC_Gly_HelixAlpha", ...]

Every named template is a real, visible constant living next to its residue,
built on import with no intermediate ``{ss: {rotamer}}`` dict and no runtime
magic.  The package ``__init__`` re-exports them via a per-residue
``from .gly import *``.  The block is regenerated in place (anything below the
generation marker is replaced), so re-running is idempotent.

Run (from the repo root)::

    /opt/envs/BioRazer/bin/python scripts/generate_internal_coord_template_named.py
"""
from __future__ import annotations

import pathlib
import sys

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

PKG = REPO_ROOT / "biorazer/database/molecule/icoor/protein/template"

#: Secondary-structure class -> identifier-safe token (mirrors ``_naming.SS_TOKEN``).
_SS_TOKEN = {
    "alpha-helix": "HelixAlpha",
    "3-10-helix": "Helix310",
    "pi-helix": "HelixPi",
    "polyproline-II": "HelixPPII",
    "beta-strand": "Strand",
    "parallel-beta-strand": "StrandParallel",
    "antiparallel-beta-strand": "StrandAntiParallel",
    "beta-bridge": "Bridge",
    "turn": "Turn",
    "bend": "Bend",
    "coil": "Coil",
    "cis-peptide-bond": "CisPeptide",
}
_ROTAMER_TOKEN = {"g-": "gminus", "t": "trans", "g+": "gplus"}


def _template_token(resn, ss, rotamer="canonical"):
    """Identifier-safe ``IC_<Res>_<SS>[_<rotamer>]`` token (mirrors
    ``_naming.template_token``)."""
    tok = f"IC_{resn.title()}_{_SS_TOKEN[ss]}"
    if rotamer != "canonical":
        tok += "_" + "_".join(_ROTAMER_TOKEN[p] for p in rotamer.split("/"))
    return tok

GEN_MARKER = "# === GENERATED IC_* named templates - do not edit below ==="

_HEADER = '''# -*- coding: utf-8 -*-
"""Ideal {secondary-structure x rotamer} {names[aa]} {var} ``InternalCoord`` templates.

Each template is the ideal heavy-atom topology for {aa} anchored at
{{N, CA, C}}, with the side chain grown to ``rotamer`` and the
``phi``/``psi``/``omega`` of ``ss`` carried on the instance.  See
:mod:`biorazer.database.molecule.icoor.protein.template` and
:mod:`._builder` for the build rule and the per-conformer caveat.
"""

from biorazer.database.molecule.icoor.protein.template._builder import build_template
'''


def _generate_block(resn):
    """The named-template block (everything below GEN_MARKER) for one residue."""
    tokens = []
    lines = [GEN_MARKER]
    lines.append(
        "# Regenerate with: /opt/envs/BioRazer/bin/python "
        "scripts/generate_internal_coord_template_named.py")
    lines.append("")
    for ss in _SS_CLASSES:
        for rot in _rotamer_names(resn):
            tok = _template_token(resn, ss, rot)
            tokens.append(tok)
            lines.append(f"{tok} = build_template({resn!r}, {ss!r}, {rot!r})")
    lines.append("")
    lines.append("__all__ = [")
    for tok in sorted(set(tokens)):
        lines.append(f"    {tok!r},")
    lines.append("]")
    return "\n".join(lines) + "\n"


def _rotamer_names(resn):
    # mirror biorazer..._builder.rotamer_names without importing the package
    # (the package __init__ imports every residue module -> circular import
    # while we're mid-regeneration).  Build the list from the common-chi axes.
    n = _COMMON_CHI_AXES.get(resn, 0)
    if n == 0:
        return ["canonical"]
    if n == 1:
        return ["canonical"] + list(_ROT_BIN)
    return ["canonical"] + [f"{a}/{b}" for a in _ROT_BIN for b in _ROT_BIN]


def main():
    total = 0
    for resn in _RESIDUES:
        path = PKG / f"{resn.lower()}.py"
        text = _HEADER
        block = _generate_block(resn)
        path.write_text(text + "\n" + block, encoding="utf-8")
        total += len(_rotamer_names(resn)) * len(_SS_CLASSES)
    print(f"wrote IC_* blocks into {len(_RESIDUES)} residue modules "
          f"({total} named templates)")


# SS classes (keys of SS_BB_TORSION_ANGLE), rotamer bin centers, and the
# standard amino-acid residues, kept in sync with
# biorazer.database.molecule.icoor.protein.  A residue's named templates
# enumerate all SS x common-rotamer combos, mirroring the build_template sweep.
_RESIDUES = (
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
)
_SS_CLASSES = (
    "alpha-helix", "3-10-helix", "pi-helix", "polyproline-II", "beta-strand",
    "parallel-beta-strand", "antiparallel-beta-strand", "beta-bridge", "turn",
    "bend", "coil", "cis-peptide-bond",
)
_ROT_BIN = ("g-", "t", "g+")
_COMMON_CHI_AXES = {
    "ALA": 0, "GLY": 0, "PRO": 0,
    "CYS": 1, "SER": 1, "THR": 1, "VAL": 1,
    "ASN": 2, "ASP": 2, "GLN": 2, "GLU": 2, "HIS": 2, "ILE": 2, "LEU": 2,
    "LYS": 2, "MET": 2, "PHE": 2, "TRP": 2, "TYR": 2, "ARG": 2,
}


if __name__ == "__main__":
    main()