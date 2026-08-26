# -*- coding: utf-8 -*-
"""Ideal {secondary-structure x rotamer} {names[aa]} {var} ``InternalCoord`` templates.

Each template is the ideal heavy-atom topology for {aa} anchored at
{{N, CA, C}}, with the side chain grown to ``rotamer`` and the
``phi``/``psi``/``omega`` of ``ss`` carried on the instance.  See
:mod:`biorazer.database.internal_coord_template` and
:mod:`._builder` for the build rule and the per-conformer caveat.
"""

from biorazer.database.internal_coord_template._builder import make_residue_templates

RESN = "GLY"
#: {ss: {rotamer: :class:`~biorazer.structure.objects.InternalCoord`}}
GLY_TEMPLATES = make_residue_templates(RESN)

# === GENERATED IC_* named templates - do not edit below ===
# Regenerate with: /opt/envs/BioRazer/bin/python scripts/generate_internal_coord_template_named.py

IC_Gly_HelixAlpha = GLY_TEMPLATES['alpha-helix']['canonical']
IC_Gly_Helix310 = GLY_TEMPLATES['3-10-helix']['canonical']
IC_Gly_HelixPi = GLY_TEMPLATES['pi-helix']['canonical']
IC_Gly_HelixPPII = GLY_TEMPLATES['polyproline-II']['canonical']
IC_Gly_Strand = GLY_TEMPLATES['beta-strand']['canonical']
IC_Gly_StrandParallel = GLY_TEMPLATES['parallel-beta-strand']['canonical']
IC_Gly_StrandAntiParallel = GLY_TEMPLATES['antiparallel-beta-strand']['canonical']
IC_Gly_Bridge = GLY_TEMPLATES['beta-bridge']['canonical']
IC_Gly_Turn = GLY_TEMPLATES['turn']['canonical']
IC_Gly_Bend = GLY_TEMPLATES['bend']['canonical']
IC_Gly_Coil = GLY_TEMPLATES['coil']['canonical']
IC_Gly_CisPeptide = GLY_TEMPLATES['cis-peptide-bond']['canonical']

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
