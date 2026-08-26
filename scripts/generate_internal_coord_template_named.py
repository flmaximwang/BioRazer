# -*- coding: utf-8 -*-
"""Generate the literal ``IC_<Res>_<SS>[_<rotamer>]`` named templates.

Appends, **into each per-residue module**, a block of literal constants that
alias its own ``<RESN>_TEMPLATES`` entries, e.g. ``gly.py``::

    IC_Gly_HelixAlpha = GLY_TEMPLATES['alpha-helix']['canonical']
    ...
    __all__ = ["IC_Gly_HelixAlpha", ...]

so every named template is a real, visible constant living next to its residue
-- no runtime magic.  The package ``__init__`` re-exports them via a per-residue
``from .gly import *``.  The block is regenerated in place (anything below the
generation marker is replaced), so re-running is idempotent.

Run (from the repo root)::

    /opt/envs/BioRazer/bin/python scripts/generate_internal_coord_template_named.py
"""
from __future__ import annotations

import importlib
import pathlib

from biorazer.database.internal_coord_template import RESIDUES
from biorazer.database.internal_coord_template._naming import template_token

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
PKG = REPO_ROOT / "biorazer/database/internal_coord_template"

GEN_MARKER = "# === GENERATED IC_* named templates - do not edit below ==="


def _generate_block(resn, templates):
    tokens = []
    lines = [GEN_MARKER]
    lines.append(
        "# Regenerate with: /opt/envs/BioRazer/bin/python "
        "scripts/generate_internal_coord_template_named.py")
    lines.append("")
    var = f"{resn}_TEMPLATES"
    for ss, rot_dict in templates.items():
        for rot in rot_dict:
            tok = template_token(resn, ss, rot)
            tokens.append(tok)
            lines.append(f"{tok} = {var}[{ss!r}][{rot!r}]")
    lines.append("")
    lines.append("__all__ = [")
    for tok in sorted(set(tokens)):
        lines.append(f"    {tok!r},")
    lines.append("]")
    return "\n".join(lines) + "\n"


def main():
    total = 0
    for resn in RESIDUES:
        path = PKG / f"{resn.lower()}.py"
        mod = importlib.import_module(
            f"biorazer.database.internal_coord_template.{resn.lower()}")
        templates = getattr(mod, f"{resn}_TEMPLATES")
        text = path.read_text(encoding="utf-8")
        idx = text.find(GEN_MARKER)
        if idx != -1:
            text = text[:idx].rstrip() + "\n"
        else:
            text = text.rstrip() + "\n"
        block = _generate_block(resn, templates)
        path.write_text(text + "\n" + block, encoding="utf-8")
        total += sum(len(rots) for rots in templates.values())
    print(f"wrote IC_* blocks into {len(RESIDUES)} residue modules "
          f"({total} named templates)")


if __name__ == "__main__":
    main()