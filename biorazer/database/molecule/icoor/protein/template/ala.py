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

IC_Ala_HelixAlpha = build_template('ALA', 'alpha-helix', 'canonical')
IC_Ala_Helix310 = build_template('ALA', '3-10-helix', 'canonical')
IC_Ala_HelixPi = build_template('ALA', 'pi-helix', 'canonical')
IC_Ala_HelixPPII = build_template('ALA', 'polyproline-II', 'canonical')
IC_Ala_Strand = build_template('ALA', 'beta-strand', 'canonical')
IC_Ala_StrandParallel = build_template('ALA', 'parallel-beta-strand', 'canonical')
IC_Ala_StrandAntiParallel = build_template('ALA', 'antiparallel-beta-strand', 'canonical')
IC_Ala_Bridge = build_template('ALA', 'beta-bridge', 'canonical')
IC_Ala_Turn = build_template('ALA', 'turn', 'canonical')
IC_Ala_Bend = build_template('ALA', 'bend', 'canonical')
IC_Ala_Coil = build_template('ALA', 'coil', 'canonical')
IC_Ala_CisPeptide = build_template('ALA', 'cis-peptide-bond', 'canonical')

__all__ = [
    'IC_Ala_Bend',
    'IC_Ala_Bridge',
    'IC_Ala_CisPeptide',
    'IC_Ala_Coil',
    'IC_Ala_Helix310',
    'IC_Ala_HelixAlpha',
    'IC_Ala_HelixPPII',
    'IC_Ala_HelixPi',
    'IC_Ala_Strand',
    'IC_Ala_StrandAntiParallel',
    'IC_Ala_StrandParallel',
    'IC_Ala_Turn',
]
