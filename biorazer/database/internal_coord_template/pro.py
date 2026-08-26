# -*- coding: utf-8 -*-
"""Ideal {secondary-structure x rotamer} {names[aa]} {var} ``InternalCoord`` templates.

Each template is the ideal heavy-atom topology for {aa} anchored at
{{N, CA, C}}, with the side chain grown to ``rotamer`` and the
``phi``/``psi``/``omega`` of ``ss`` carried on the instance.  See
:mod:`biorazer.database.internal_coord_template` and
:mod:`._builder` for the build rule and the per-conformer caveat.
"""

from biorazer.database.internal_coord_template._builder import make_residue_templates

RESN = "PRO"
#: {ss: {rotamer: :class:`~biorazer.structure.objects.InternalCoord`}}
PRO_TEMPLATES = make_residue_templates(RESN)

# === GENERATED IC_* named templates - do not edit below ===
# Regenerate with: /opt/envs/BioRazer/bin/python scripts/generate_internal_coord_template_named.py

IC_Pro_HelixAlpha = PRO_TEMPLATES['alpha-helix']['canonical']
IC_Pro_Helix310 = PRO_TEMPLATES['3-10-helix']['canonical']
IC_Pro_HelixPi = PRO_TEMPLATES['pi-helix']['canonical']
IC_Pro_HelixPPII = PRO_TEMPLATES['polyproline-II']['canonical']
IC_Pro_Strand = PRO_TEMPLATES['beta-strand']['canonical']
IC_Pro_StrandParallel = PRO_TEMPLATES['parallel-beta-strand']['canonical']
IC_Pro_StrandAntiParallel = PRO_TEMPLATES['antiparallel-beta-strand']['canonical']
IC_Pro_Bridge = PRO_TEMPLATES['beta-bridge']['canonical']
IC_Pro_Turn = PRO_TEMPLATES['turn']['canonical']
IC_Pro_Bend = PRO_TEMPLATES['bend']['canonical']
IC_Pro_Coil = PRO_TEMPLATES['coil']['canonical']
IC_Pro_CisPeptide = PRO_TEMPLATES['cis-peptide-bond']['canonical']

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
