# -*- coding: utf-8 -*-
"""Ideal {secondary-structure x rotamer} {names[aa]} {var} ``InternalCoord`` templates.

Each template is the ideal heavy-atom topology for {aa} anchored at
{{N, CA, C}}, with the side chain grown to ``rotamer`` and the
``phi``/``psi``/``omega`` of ``ss`` carried on the instance.  See
:mod:`biorazer.database.internal_coord_template` and
:mod:`._builder` for the build rule and the per-conformer caveat.
"""

from biorazer.database.molecule.icoor.protein.template._builder import make_residue_templates

RESN = "VAL"
#: {ss: {rotamer: :class:`~biorazer.structure.objects.InternalCoord`}}
VAL_TEMPLATES = make_residue_templates(RESN)

# === GENERATED IC_* named templates - do not edit below ===
# Regenerate with: /opt/envs/BioRazer/bin/python scripts/generate_internal_coord_template_named.py

IC_Val_HelixAlpha = VAL_TEMPLATES['alpha-helix']['canonical']
IC_Val_HelixAlpha_gminus = VAL_TEMPLATES['alpha-helix']['g-']
IC_Val_HelixAlpha_trans = VAL_TEMPLATES['alpha-helix']['t']
IC_Val_HelixAlpha_gplus = VAL_TEMPLATES['alpha-helix']['g+']
IC_Val_Helix310 = VAL_TEMPLATES['3-10-helix']['canonical']
IC_Val_Helix310_gminus = VAL_TEMPLATES['3-10-helix']['g-']
IC_Val_Helix310_trans = VAL_TEMPLATES['3-10-helix']['t']
IC_Val_Helix310_gplus = VAL_TEMPLATES['3-10-helix']['g+']
IC_Val_HelixPi = VAL_TEMPLATES['pi-helix']['canonical']
IC_Val_HelixPi_gminus = VAL_TEMPLATES['pi-helix']['g-']
IC_Val_HelixPi_trans = VAL_TEMPLATES['pi-helix']['t']
IC_Val_HelixPi_gplus = VAL_TEMPLATES['pi-helix']['g+']
IC_Val_HelixPPII = VAL_TEMPLATES['polyproline-II']['canonical']
IC_Val_HelixPPII_gminus = VAL_TEMPLATES['polyproline-II']['g-']
IC_Val_HelixPPII_trans = VAL_TEMPLATES['polyproline-II']['t']
IC_Val_HelixPPII_gplus = VAL_TEMPLATES['polyproline-II']['g+']
IC_Val_Strand = VAL_TEMPLATES['beta-strand']['canonical']
IC_Val_Strand_gminus = VAL_TEMPLATES['beta-strand']['g-']
IC_Val_Strand_trans = VAL_TEMPLATES['beta-strand']['t']
IC_Val_Strand_gplus = VAL_TEMPLATES['beta-strand']['g+']
IC_Val_StrandParallel = VAL_TEMPLATES['parallel-beta-strand']['canonical']
IC_Val_StrandParallel_gminus = VAL_TEMPLATES['parallel-beta-strand']['g-']
IC_Val_StrandParallel_trans = VAL_TEMPLATES['parallel-beta-strand']['t']
IC_Val_StrandParallel_gplus = VAL_TEMPLATES['parallel-beta-strand']['g+']
IC_Val_StrandAntiParallel = VAL_TEMPLATES['antiparallel-beta-strand']['canonical']
IC_Val_StrandAntiParallel_gminus = VAL_TEMPLATES['antiparallel-beta-strand']['g-']
IC_Val_StrandAntiParallel_trans = VAL_TEMPLATES['antiparallel-beta-strand']['t']
IC_Val_StrandAntiParallel_gplus = VAL_TEMPLATES['antiparallel-beta-strand']['g+']
IC_Val_Bridge = VAL_TEMPLATES['beta-bridge']['canonical']
IC_Val_Bridge_gminus = VAL_TEMPLATES['beta-bridge']['g-']
IC_Val_Bridge_trans = VAL_TEMPLATES['beta-bridge']['t']
IC_Val_Bridge_gplus = VAL_TEMPLATES['beta-bridge']['g+']
IC_Val_Turn = VAL_TEMPLATES['turn']['canonical']
IC_Val_Turn_gminus = VAL_TEMPLATES['turn']['g-']
IC_Val_Turn_trans = VAL_TEMPLATES['turn']['t']
IC_Val_Turn_gplus = VAL_TEMPLATES['turn']['g+']
IC_Val_Bend = VAL_TEMPLATES['bend']['canonical']
IC_Val_Bend_gminus = VAL_TEMPLATES['bend']['g-']
IC_Val_Bend_trans = VAL_TEMPLATES['bend']['t']
IC_Val_Bend_gplus = VAL_TEMPLATES['bend']['g+']
IC_Val_Coil = VAL_TEMPLATES['coil']['canonical']
IC_Val_Coil_gminus = VAL_TEMPLATES['coil']['g-']
IC_Val_Coil_trans = VAL_TEMPLATES['coil']['t']
IC_Val_Coil_gplus = VAL_TEMPLATES['coil']['g+']
IC_Val_CisPeptide = VAL_TEMPLATES['cis-peptide-bond']['canonical']
IC_Val_CisPeptide_gminus = VAL_TEMPLATES['cis-peptide-bond']['g-']
IC_Val_CisPeptide_trans = VAL_TEMPLATES['cis-peptide-bond']['t']
IC_Val_CisPeptide_gplus = VAL_TEMPLATES['cis-peptide-bond']['g+']

__all__ = [
    'IC_Val_Bend',
    'IC_Val_Bend_gminus',
    'IC_Val_Bend_gplus',
    'IC_Val_Bend_trans',
    'IC_Val_Bridge',
    'IC_Val_Bridge_gminus',
    'IC_Val_Bridge_gplus',
    'IC_Val_Bridge_trans',
    'IC_Val_CisPeptide',
    'IC_Val_CisPeptide_gminus',
    'IC_Val_CisPeptide_gplus',
    'IC_Val_CisPeptide_trans',
    'IC_Val_Coil',
    'IC_Val_Coil_gminus',
    'IC_Val_Coil_gplus',
    'IC_Val_Coil_trans',
    'IC_Val_Helix310',
    'IC_Val_Helix310_gminus',
    'IC_Val_Helix310_gplus',
    'IC_Val_Helix310_trans',
    'IC_Val_HelixAlpha',
    'IC_Val_HelixAlpha_gminus',
    'IC_Val_HelixAlpha_gplus',
    'IC_Val_HelixAlpha_trans',
    'IC_Val_HelixPPII',
    'IC_Val_HelixPPII_gminus',
    'IC_Val_HelixPPII_gplus',
    'IC_Val_HelixPPII_trans',
    'IC_Val_HelixPi',
    'IC_Val_HelixPi_gminus',
    'IC_Val_HelixPi_gplus',
    'IC_Val_HelixPi_trans',
    'IC_Val_Strand',
    'IC_Val_StrandAntiParallel',
    'IC_Val_StrandAntiParallel_gminus',
    'IC_Val_StrandAntiParallel_gplus',
    'IC_Val_StrandAntiParallel_trans',
    'IC_Val_StrandParallel',
    'IC_Val_StrandParallel_gminus',
    'IC_Val_StrandParallel_gplus',
    'IC_Val_StrandParallel_trans',
    'IC_Val_Strand_gminus',
    'IC_Val_Strand_gplus',
    'IC_Val_Strand_trans',
    'IC_Val_Turn',
    'IC_Val_Turn_gminus',
    'IC_Val_Turn_gplus',
    'IC_Val_Turn_trans',
]
