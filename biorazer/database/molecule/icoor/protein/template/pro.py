# -*- coding: utf-8 -*-
"""Ideal {secondary-structure x rotamer} {names[aa]} {var} ``InternalCoord`` templates.

Each template is the ideal heavy-atom topology for {aa} anchored at
{{N, CA, C}}, with the side chain grown to ``rotamer`` and the
``phi``/``psi``/``omega`` of ``ss`` carried on the instance.  See
:mod:`biorazer.database.molecule.icoor.protein.template` and
:mod:`._builder` for the build rule and the per-conformer caveat.
"""

from biorazer.database.molecule.icoor.protein.template._builder import build_template

# === GENERATED IC_* named templates - do not edit below ===
# Regenerate with: /opt/envs/BioRazer/bin/python scripts/generate_internal_coord_template_named.py

IC_Pro_HelixAlpha = build_template('PRO', 'alpha-helix', 'canonical')
IC_Pro_Helix310 = build_template('PRO', '3-10-helix', 'canonical')
IC_Pro_HelixPi = build_template('PRO', 'pi-helix', 'canonical')
IC_Pro_HelixPPII = build_template('PRO', 'polyproline-II', 'canonical')
IC_Pro_Strand = build_template('PRO', 'beta-strand', 'canonical')
IC_Pro_StrandParallel = build_template('PRO', 'parallel-beta-strand', 'canonical')
IC_Pro_StrandAntiParallel = build_template('PRO', 'antiparallel-beta-strand', 'canonical')
IC_Pro_Bridge = build_template('PRO', 'beta-bridge', 'canonical')
IC_Pro_Turn = build_template('PRO', 'turn', 'canonical')
IC_Pro_Bend = build_template('PRO', 'bend', 'canonical')
IC_Pro_Coil = build_template('PRO', 'coil', 'canonical')
IC_Pro_CisPeptide = build_template('PRO', 'cis-peptide-bond', 'canonical')

__all__ = [
    'IC_Pro_Bend',
    'IC_Pro_Bridge',
    'IC_Pro_CisPeptide',
    'IC_Pro_Coil',
    'IC_Pro_Helix310',
    'IC_Pro_HelixAlpha',
    'IC_Pro_HelixPPII',
    'IC_Pro_HelixPi',
    'IC_Pro_Strand',
    'IC_Pro_StrandAntiParallel',
    'IC_Pro_StrandParallel',
    'IC_Pro_Turn',
]
