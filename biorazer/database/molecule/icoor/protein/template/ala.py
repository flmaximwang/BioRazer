# -*- coding: utf-8 -*-
"""Ideal {secondary-structure x rotamer} {names[aa]} {var} ``InternalCoord`` templates.

Each template is the ideal heavy-atom topology for {aa} anchored at
{{N, CA, C}}, with the side chain grown to ``rotamer`` and the
``phi``/``psi``/``omega`` of ``ss`` carried on the instance.  See
:mod:`biorazer.database.internal_coord_template` and
:mod:`._builder` for the build rule and the per-conformer caveat.
"""

from biorazer.database.molecule.icoor.protein.template._builder import make_residue_templates

RESN = "ALA"
#: {ss: {rotamer: :class:`~biorazer.structure.objects.InternalCoord`}}
ALA_TEMPLATES = make_residue_templates(RESN)

# === GENERATED IC_* named templates - do not edit below ===
# Regenerate with: /opt/envs/BioRazer/bin/python scripts/generate_internal_coord_template_named.py

IC_Ala_HelixAlpha = ALA_TEMPLATES['alpha-helix']['canonical']
IC_Ala_Helix310 = ALA_TEMPLATES['3-10-helix']['canonical']
IC_Ala_HelixPi = ALA_TEMPLATES['pi-helix']['canonical']
IC_Ala_HelixPPII = ALA_TEMPLATES['polyproline-II']['canonical']
IC_Ala_Strand = ALA_TEMPLATES['beta-strand']['canonical']
IC_Ala_StrandParallel = ALA_TEMPLATES['parallel-beta-strand']['canonical']
IC_Ala_StrandAntiParallel = ALA_TEMPLATES['antiparallel-beta-strand']['canonical']
IC_Ala_Bridge = ALA_TEMPLATES['beta-bridge']['canonical']
IC_Ala_Turn = ALA_TEMPLATES['turn']['canonical']
IC_Ala_Bend = ALA_TEMPLATES['bend']['canonical']
IC_Ala_Coil = ALA_TEMPLATES['coil']['canonical']
IC_Ala_CisPeptide = ALA_TEMPLATES['cis-peptide-bond']['canonical']

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
