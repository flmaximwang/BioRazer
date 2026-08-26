# -*- coding: utf-8 -*-
"""Ideal {secondary-structure x rotamer} {names[aa]} {var} ``InternalCoord`` templates.

Each template is the ideal heavy-atom topology for {aa} anchored at
{{N, CA, C}}, with the side chain grown to ``rotamer`` and the
``phi``/``psi``/``omega`` of ``ss`` carried on the instance.  See
:mod:`biorazer.database.internal_coord_template` and
:mod:`._builder` for the build rule and the per-conformer caveat.
"""

from biorazer.database.internal_coord_template._builder import make_residue_templates

RESN = "CYS"
#: {ss: {rotamer: :class:`~biorazer.structure.objects.InternalCoord`}}
CYS_TEMPLATES = make_residue_templates(RESN)

# === GENERATED IC_* named templates - do not edit below ===
# Regenerate with: /opt/envs/BioRazer/bin/python scripts/generate_internal_coord_template_named.py

IC_Cys_HelixAlpha = CYS_TEMPLATES['alpha-helix']['canonical']
IC_Cys_HelixAlpha_gminus = CYS_TEMPLATES['alpha-helix']['g-']
IC_Cys_HelixAlpha_trans = CYS_TEMPLATES['alpha-helix']['t']
IC_Cys_HelixAlpha_gplus = CYS_TEMPLATES['alpha-helix']['g+']
IC_Cys_Helix310 = CYS_TEMPLATES['3-10-helix']['canonical']
IC_Cys_Helix310_gminus = CYS_TEMPLATES['3-10-helix']['g-']
IC_Cys_Helix310_trans = CYS_TEMPLATES['3-10-helix']['t']
IC_Cys_Helix310_gplus = CYS_TEMPLATES['3-10-helix']['g+']
IC_Cys_HelixPi = CYS_TEMPLATES['pi-helix']['canonical']
IC_Cys_HelixPi_gminus = CYS_TEMPLATES['pi-helix']['g-']
IC_Cys_HelixPi_trans = CYS_TEMPLATES['pi-helix']['t']
IC_Cys_HelixPi_gplus = CYS_TEMPLATES['pi-helix']['g+']
IC_Cys_HelixPPII = CYS_TEMPLATES['polyproline-II']['canonical']
IC_Cys_HelixPPII_gminus = CYS_TEMPLATES['polyproline-II']['g-']
IC_Cys_HelixPPII_trans = CYS_TEMPLATES['polyproline-II']['t']
IC_Cys_HelixPPII_gplus = CYS_TEMPLATES['polyproline-II']['g+']
IC_Cys_Strand = CYS_TEMPLATES['beta-strand']['canonical']
IC_Cys_Strand_gminus = CYS_TEMPLATES['beta-strand']['g-']
IC_Cys_Strand_trans = CYS_TEMPLATES['beta-strand']['t']
IC_Cys_Strand_gplus = CYS_TEMPLATES['beta-strand']['g+']
IC_Cys_StrandParallel = CYS_TEMPLATES['parallel-beta-strand']['canonical']
IC_Cys_StrandParallel_gminus = CYS_TEMPLATES['parallel-beta-strand']['g-']
IC_Cys_StrandParallel_trans = CYS_TEMPLATES['parallel-beta-strand']['t']
IC_Cys_StrandParallel_gplus = CYS_TEMPLATES['parallel-beta-strand']['g+']
IC_Cys_StrandAntiParallel = CYS_TEMPLATES['antiparallel-beta-strand']['canonical']
IC_Cys_StrandAntiParallel_gminus = CYS_TEMPLATES['antiparallel-beta-strand']['g-']
IC_Cys_StrandAntiParallel_trans = CYS_TEMPLATES['antiparallel-beta-strand']['t']
IC_Cys_StrandAntiParallel_gplus = CYS_TEMPLATES['antiparallel-beta-strand']['g+']
IC_Cys_Bridge = CYS_TEMPLATES['beta-bridge']['canonical']
IC_Cys_Bridge_gminus = CYS_TEMPLATES['beta-bridge']['g-']
IC_Cys_Bridge_trans = CYS_TEMPLATES['beta-bridge']['t']
IC_Cys_Bridge_gplus = CYS_TEMPLATES['beta-bridge']['g+']
IC_Cys_Turn = CYS_TEMPLATES['turn']['canonical']
IC_Cys_Turn_gminus = CYS_TEMPLATES['turn']['g-']
IC_Cys_Turn_trans = CYS_TEMPLATES['turn']['t']
IC_Cys_Turn_gplus = CYS_TEMPLATES['turn']['g+']
IC_Cys_Bend = CYS_TEMPLATES['bend']['canonical']
IC_Cys_Bend_gminus = CYS_TEMPLATES['bend']['g-']
IC_Cys_Bend_trans = CYS_TEMPLATES['bend']['t']
IC_Cys_Bend_gplus = CYS_TEMPLATES['bend']['g+']
IC_Cys_Coil = CYS_TEMPLATES['coil']['canonical']
IC_Cys_Coil_gminus = CYS_TEMPLATES['coil']['g-']
IC_Cys_Coil_trans = CYS_TEMPLATES['coil']['t']
IC_Cys_Coil_gplus = CYS_TEMPLATES['coil']['g+']
IC_Cys_CisPeptide = CYS_TEMPLATES['cis-peptide-bond']['canonical']
IC_Cys_CisPeptide_gminus = CYS_TEMPLATES['cis-peptide-bond']['g-']
IC_Cys_CisPeptide_trans = CYS_TEMPLATES['cis-peptide-bond']['t']
IC_Cys_CisPeptide_gplus = CYS_TEMPLATES['cis-peptide-bond']['g+']

__all__ = [
    'IC_Cys_Bend',
    'IC_Cys_Bend_gminus',
    'IC_Cys_Bend_gplus',
    'IC_Cys_Bend_trans',
    'IC_Cys_Bridge',
    'IC_Cys_Bridge_gminus',
    'IC_Cys_Bridge_gplus',
    'IC_Cys_Bridge_trans',
    'IC_Cys_CisPeptide',
    'IC_Cys_CisPeptide_gminus',
    'IC_Cys_CisPeptide_gplus',
    'IC_Cys_CisPeptide_trans',
    'IC_Cys_Coil',
    'IC_Cys_Coil_gminus',
    'IC_Cys_Coil_gplus',
    'IC_Cys_Coil_trans',
    'IC_Cys_Helix310',
    'IC_Cys_Helix310_gminus',
    'IC_Cys_Helix310_gplus',
    'IC_Cys_Helix310_trans',
    'IC_Cys_HelixAlpha',
    'IC_Cys_HelixAlpha_gminus',
    'IC_Cys_HelixAlpha_gplus',
    'IC_Cys_HelixAlpha_trans',
    'IC_Cys_HelixPPII',
    'IC_Cys_HelixPPII_gminus',
    'IC_Cys_HelixPPII_gplus',
    'IC_Cys_HelixPPII_trans',
    'IC_Cys_HelixPi',
    'IC_Cys_HelixPi_gminus',
    'IC_Cys_HelixPi_gplus',
    'IC_Cys_HelixPi_trans',
    'IC_Cys_Strand',
    'IC_Cys_StrandAntiParallel',
    'IC_Cys_StrandAntiParallel_gminus',
    'IC_Cys_StrandAntiParallel_gplus',
    'IC_Cys_StrandAntiParallel_trans',
    'IC_Cys_StrandParallel',
    'IC_Cys_StrandParallel_gminus',
    'IC_Cys_StrandParallel_gplus',
    'IC_Cys_StrandParallel_trans',
    'IC_Cys_Strand_gminus',
    'IC_Cys_Strand_gplus',
    'IC_Cys_Strand_trans',
    'IC_Cys_Turn',
    'IC_Cys_Turn_gminus',
    'IC_Cys_Turn_gplus',
    'IC_Cys_Turn_trans',
]
