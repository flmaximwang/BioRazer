# -*- coding: utf-8 -*-
"""Ideal {secondary-structure x rotamer} {names[aa]} {var} ``InternalCoord`` templates.

Each template is the ideal heavy-atom topology for {aa} anchored at
{{N, CA, C}}, with the side chain grown to ``rotamer`` and the
``phi``/``psi``/``omega`` of ``ss`` carried on the instance.  See
:mod:`biorazer.database.molecule.icoor.protein.template` and
:mod:`._builder` for the build rule and the per-conformer caveat.
"""

from biorazer.database.molecule.icoor.protein.template._builder import build_template

RESN = "GLY"

# === GENERATED IC_* named templates - do not edit below ===
# Regenerate with: /opt/envs/BioRazer/bin/python scripts/generate_internal_coord_template_named.py

IC_Gly_HelixAlpha = build_template('GLY', 'alpha-helix', 'canonical')
IC_Gly_Helix310 = build_template('GLY', '3-10-helix', 'canonical')
IC_Gly_HelixPi = build_template('GLY', 'pi-helix', 'canonical')
IC_Gly_HelixPPII = build_template('GLY', 'polyproline-II', 'canonical')
IC_Gly_Strand = build_template('GLY', 'beta-strand', 'canonical')
IC_Gly_StrandParallel = build_template('GLY', 'parallel-beta-strand', 'canonical')
IC_Gly_StrandAntiParallel = build_template('GLY', 'antiparallel-beta-strand', 'canonical')
IC_Gly_Bridge = build_template('GLY', 'beta-bridge', 'canonical')
IC_Gly_Turn = build_template('GLY', 'turn', 'canonical')
IC_Gly_Bend = build_template('GLY', 'bend', 'canonical')
IC_Gly_Coil = build_template('GLY', 'coil', 'canonical')
IC_Gly_CisPeptide = build_template('GLY', 'cis-peptide-bond', 'canonical')

__all__ = [
    'IC_Gly_Bend',
    'IC_Gly_Bridge',
    'IC_Gly_CisPeptide',
    'IC_Gly_Coil',
    'IC_Gly_Helix310',
    'IC_Gly_HelixAlpha',
    'IC_Gly_HelixPPII',
    'IC_Gly_HelixPi',
    'IC_Gly_Strand',
    'IC_Gly_StrandAntiParallel',
    'IC_Gly_StrandParallel',
    'IC_Gly_Turn',
]
