# -*- coding: utf-8 -*-
"""Ideal {secondary-structure x rotamer} {names[aa]} {var} ``InternalCoord`` templates.

Each template is the ideal heavy-atom topology for {aa} anchored at
{{N, CA, C}}, with the side chain grown to ``rotamer`` and the
``phi``/``psi``/``omega`` of ``ss`` carried on the instance.  See
:mod:`biorazer.database.internal_coord_template` and
:mod:`._builder` for the build rule and the per-conformer caveat.
"""

from biorazer.database.molecule.icoor.protein.template._builder import make_residue_templates

RESN = "ASP"
#: {ss: {rotamer: :class:`~biorazer.structure.objects.InternalCoord`}}
ASP_TEMPLATES = make_residue_templates(RESN)

# === GENERATED IC_* named templates - do not edit below ===
# Regenerate with: /opt/envs/BioRazer/bin/python scripts/generate_internal_coord_template_named.py

IC_Asp_HelixAlpha = ASP_TEMPLATES['alpha-helix']['canonical']
IC_Asp_HelixAlpha_gminus_gminus = ASP_TEMPLATES['alpha-helix']['g-/g-']
IC_Asp_HelixAlpha_gminus_trans = ASP_TEMPLATES['alpha-helix']['g-/t']
IC_Asp_HelixAlpha_gminus_gplus = ASP_TEMPLATES['alpha-helix']['g-/g+']
IC_Asp_HelixAlpha_trans_gminus = ASP_TEMPLATES['alpha-helix']['t/g-']
IC_Asp_HelixAlpha_trans_trans = ASP_TEMPLATES['alpha-helix']['t/t']
IC_Asp_HelixAlpha_trans_gplus = ASP_TEMPLATES['alpha-helix']['t/g+']
IC_Asp_HelixAlpha_gplus_gminus = ASP_TEMPLATES['alpha-helix']['g+/g-']
IC_Asp_HelixAlpha_gplus_trans = ASP_TEMPLATES['alpha-helix']['g+/t']
IC_Asp_HelixAlpha_gplus_gplus = ASP_TEMPLATES['alpha-helix']['g+/g+']
IC_Asp_Helix310 = ASP_TEMPLATES['3-10-helix']['canonical']
IC_Asp_Helix310_gminus_gminus = ASP_TEMPLATES['3-10-helix']['g-/g-']
IC_Asp_Helix310_gminus_trans = ASP_TEMPLATES['3-10-helix']['g-/t']
IC_Asp_Helix310_gminus_gplus = ASP_TEMPLATES['3-10-helix']['g-/g+']
IC_Asp_Helix310_trans_gminus = ASP_TEMPLATES['3-10-helix']['t/g-']
IC_Asp_Helix310_trans_trans = ASP_TEMPLATES['3-10-helix']['t/t']
IC_Asp_Helix310_trans_gplus = ASP_TEMPLATES['3-10-helix']['t/g+']
IC_Asp_Helix310_gplus_gminus = ASP_TEMPLATES['3-10-helix']['g+/g-']
IC_Asp_Helix310_gplus_trans = ASP_TEMPLATES['3-10-helix']['g+/t']
IC_Asp_Helix310_gplus_gplus = ASP_TEMPLATES['3-10-helix']['g+/g+']
IC_Asp_HelixPi = ASP_TEMPLATES['pi-helix']['canonical']
IC_Asp_HelixPi_gminus_gminus = ASP_TEMPLATES['pi-helix']['g-/g-']
IC_Asp_HelixPi_gminus_trans = ASP_TEMPLATES['pi-helix']['g-/t']
IC_Asp_HelixPi_gminus_gplus = ASP_TEMPLATES['pi-helix']['g-/g+']
IC_Asp_HelixPi_trans_gminus = ASP_TEMPLATES['pi-helix']['t/g-']
IC_Asp_HelixPi_trans_trans = ASP_TEMPLATES['pi-helix']['t/t']
IC_Asp_HelixPi_trans_gplus = ASP_TEMPLATES['pi-helix']['t/g+']
IC_Asp_HelixPi_gplus_gminus = ASP_TEMPLATES['pi-helix']['g+/g-']
IC_Asp_HelixPi_gplus_trans = ASP_TEMPLATES['pi-helix']['g+/t']
IC_Asp_HelixPi_gplus_gplus = ASP_TEMPLATES['pi-helix']['g+/g+']
IC_Asp_HelixPPII = ASP_TEMPLATES['polyproline-II']['canonical']
IC_Asp_HelixPPII_gminus_gminus = ASP_TEMPLATES['polyproline-II']['g-/g-']
IC_Asp_HelixPPII_gminus_trans = ASP_TEMPLATES['polyproline-II']['g-/t']
IC_Asp_HelixPPII_gminus_gplus = ASP_TEMPLATES['polyproline-II']['g-/g+']
IC_Asp_HelixPPII_trans_gminus = ASP_TEMPLATES['polyproline-II']['t/g-']
IC_Asp_HelixPPII_trans_trans = ASP_TEMPLATES['polyproline-II']['t/t']
IC_Asp_HelixPPII_trans_gplus = ASP_TEMPLATES['polyproline-II']['t/g+']
IC_Asp_HelixPPII_gplus_gminus = ASP_TEMPLATES['polyproline-II']['g+/g-']
IC_Asp_HelixPPII_gplus_trans = ASP_TEMPLATES['polyproline-II']['g+/t']
IC_Asp_HelixPPII_gplus_gplus = ASP_TEMPLATES['polyproline-II']['g+/g+']
IC_Asp_Strand = ASP_TEMPLATES['beta-strand']['canonical']
IC_Asp_Strand_gminus_gminus = ASP_TEMPLATES['beta-strand']['g-/g-']
IC_Asp_Strand_gminus_trans = ASP_TEMPLATES['beta-strand']['g-/t']
IC_Asp_Strand_gminus_gplus = ASP_TEMPLATES['beta-strand']['g-/g+']
IC_Asp_Strand_trans_gminus = ASP_TEMPLATES['beta-strand']['t/g-']
IC_Asp_Strand_trans_trans = ASP_TEMPLATES['beta-strand']['t/t']
IC_Asp_Strand_trans_gplus = ASP_TEMPLATES['beta-strand']['t/g+']
IC_Asp_Strand_gplus_gminus = ASP_TEMPLATES['beta-strand']['g+/g-']
IC_Asp_Strand_gplus_trans = ASP_TEMPLATES['beta-strand']['g+/t']
IC_Asp_Strand_gplus_gplus = ASP_TEMPLATES['beta-strand']['g+/g+']
IC_Asp_StrandParallel = ASP_TEMPLATES['parallel-beta-strand']['canonical']
IC_Asp_StrandParallel_gminus_gminus = ASP_TEMPLATES['parallel-beta-strand']['g-/g-']
IC_Asp_StrandParallel_gminus_trans = ASP_TEMPLATES['parallel-beta-strand']['g-/t']
IC_Asp_StrandParallel_gminus_gplus = ASP_TEMPLATES['parallel-beta-strand']['g-/g+']
IC_Asp_StrandParallel_trans_gminus = ASP_TEMPLATES['parallel-beta-strand']['t/g-']
IC_Asp_StrandParallel_trans_trans = ASP_TEMPLATES['parallel-beta-strand']['t/t']
IC_Asp_StrandParallel_trans_gplus = ASP_TEMPLATES['parallel-beta-strand']['t/g+']
IC_Asp_StrandParallel_gplus_gminus = ASP_TEMPLATES['parallel-beta-strand']['g+/g-']
IC_Asp_StrandParallel_gplus_trans = ASP_TEMPLATES['parallel-beta-strand']['g+/t']
IC_Asp_StrandParallel_gplus_gplus = ASP_TEMPLATES['parallel-beta-strand']['g+/g+']
IC_Asp_StrandAntiParallel = ASP_TEMPLATES['antiparallel-beta-strand']['canonical']
IC_Asp_StrandAntiParallel_gminus_gminus = ASP_TEMPLATES['antiparallel-beta-strand']['g-/g-']
IC_Asp_StrandAntiParallel_gminus_trans = ASP_TEMPLATES['antiparallel-beta-strand']['g-/t']
IC_Asp_StrandAntiParallel_gminus_gplus = ASP_TEMPLATES['antiparallel-beta-strand']['g-/g+']
IC_Asp_StrandAntiParallel_trans_gminus = ASP_TEMPLATES['antiparallel-beta-strand']['t/g-']
IC_Asp_StrandAntiParallel_trans_trans = ASP_TEMPLATES['antiparallel-beta-strand']['t/t']
IC_Asp_StrandAntiParallel_trans_gplus = ASP_TEMPLATES['antiparallel-beta-strand']['t/g+']
IC_Asp_StrandAntiParallel_gplus_gminus = ASP_TEMPLATES['antiparallel-beta-strand']['g+/g-']
IC_Asp_StrandAntiParallel_gplus_trans = ASP_TEMPLATES['antiparallel-beta-strand']['g+/t']
IC_Asp_StrandAntiParallel_gplus_gplus = ASP_TEMPLATES['antiparallel-beta-strand']['g+/g+']
IC_Asp_Bridge = ASP_TEMPLATES['beta-bridge']['canonical']
IC_Asp_Bridge_gminus_gminus = ASP_TEMPLATES['beta-bridge']['g-/g-']
IC_Asp_Bridge_gminus_trans = ASP_TEMPLATES['beta-bridge']['g-/t']
IC_Asp_Bridge_gminus_gplus = ASP_TEMPLATES['beta-bridge']['g-/g+']
IC_Asp_Bridge_trans_gminus = ASP_TEMPLATES['beta-bridge']['t/g-']
IC_Asp_Bridge_trans_trans = ASP_TEMPLATES['beta-bridge']['t/t']
IC_Asp_Bridge_trans_gplus = ASP_TEMPLATES['beta-bridge']['t/g+']
IC_Asp_Bridge_gplus_gminus = ASP_TEMPLATES['beta-bridge']['g+/g-']
IC_Asp_Bridge_gplus_trans = ASP_TEMPLATES['beta-bridge']['g+/t']
IC_Asp_Bridge_gplus_gplus = ASP_TEMPLATES['beta-bridge']['g+/g+']
IC_Asp_Turn = ASP_TEMPLATES['turn']['canonical']
IC_Asp_Turn_gminus_gminus = ASP_TEMPLATES['turn']['g-/g-']
IC_Asp_Turn_gminus_trans = ASP_TEMPLATES['turn']['g-/t']
IC_Asp_Turn_gminus_gplus = ASP_TEMPLATES['turn']['g-/g+']
IC_Asp_Turn_trans_gminus = ASP_TEMPLATES['turn']['t/g-']
IC_Asp_Turn_trans_trans = ASP_TEMPLATES['turn']['t/t']
IC_Asp_Turn_trans_gplus = ASP_TEMPLATES['turn']['t/g+']
IC_Asp_Turn_gplus_gminus = ASP_TEMPLATES['turn']['g+/g-']
IC_Asp_Turn_gplus_trans = ASP_TEMPLATES['turn']['g+/t']
IC_Asp_Turn_gplus_gplus = ASP_TEMPLATES['turn']['g+/g+']
IC_Asp_Bend = ASP_TEMPLATES['bend']['canonical']
IC_Asp_Bend_gminus_gminus = ASP_TEMPLATES['bend']['g-/g-']
IC_Asp_Bend_gminus_trans = ASP_TEMPLATES['bend']['g-/t']
IC_Asp_Bend_gminus_gplus = ASP_TEMPLATES['bend']['g-/g+']
IC_Asp_Bend_trans_gminus = ASP_TEMPLATES['bend']['t/g-']
IC_Asp_Bend_trans_trans = ASP_TEMPLATES['bend']['t/t']
IC_Asp_Bend_trans_gplus = ASP_TEMPLATES['bend']['t/g+']
IC_Asp_Bend_gplus_gminus = ASP_TEMPLATES['bend']['g+/g-']
IC_Asp_Bend_gplus_trans = ASP_TEMPLATES['bend']['g+/t']
IC_Asp_Bend_gplus_gplus = ASP_TEMPLATES['bend']['g+/g+']
IC_Asp_Coil = ASP_TEMPLATES['coil']['canonical']
IC_Asp_Coil_gminus_gminus = ASP_TEMPLATES['coil']['g-/g-']
IC_Asp_Coil_gminus_trans = ASP_TEMPLATES['coil']['g-/t']
IC_Asp_Coil_gminus_gplus = ASP_TEMPLATES['coil']['g-/g+']
IC_Asp_Coil_trans_gminus = ASP_TEMPLATES['coil']['t/g-']
IC_Asp_Coil_trans_trans = ASP_TEMPLATES['coil']['t/t']
IC_Asp_Coil_trans_gplus = ASP_TEMPLATES['coil']['t/g+']
IC_Asp_Coil_gplus_gminus = ASP_TEMPLATES['coil']['g+/g-']
IC_Asp_Coil_gplus_trans = ASP_TEMPLATES['coil']['g+/t']
IC_Asp_Coil_gplus_gplus = ASP_TEMPLATES['coil']['g+/g+']
IC_Asp_CisPeptide = ASP_TEMPLATES['cis-peptide-bond']['canonical']
IC_Asp_CisPeptide_gminus_gminus = ASP_TEMPLATES['cis-peptide-bond']['g-/g-']
IC_Asp_CisPeptide_gminus_trans = ASP_TEMPLATES['cis-peptide-bond']['g-/t']
IC_Asp_CisPeptide_gminus_gplus = ASP_TEMPLATES['cis-peptide-bond']['g-/g+']
IC_Asp_CisPeptide_trans_gminus = ASP_TEMPLATES['cis-peptide-bond']['t/g-']
IC_Asp_CisPeptide_trans_trans = ASP_TEMPLATES['cis-peptide-bond']['t/t']
IC_Asp_CisPeptide_trans_gplus = ASP_TEMPLATES['cis-peptide-bond']['t/g+']
IC_Asp_CisPeptide_gplus_gminus = ASP_TEMPLATES['cis-peptide-bond']['g+/g-']
IC_Asp_CisPeptide_gplus_trans = ASP_TEMPLATES['cis-peptide-bond']['g+/t']
IC_Asp_CisPeptide_gplus_gplus = ASP_TEMPLATES['cis-peptide-bond']['g+/g+']

__all__ = [
    'IC_Asp_Bend',
    'IC_Asp_Bend_gminus_gminus',
    'IC_Asp_Bend_gminus_gplus',
    'IC_Asp_Bend_gminus_trans',
    'IC_Asp_Bend_gplus_gminus',
    'IC_Asp_Bend_gplus_gplus',
    'IC_Asp_Bend_gplus_trans',
    'IC_Asp_Bend_trans_gminus',
    'IC_Asp_Bend_trans_gplus',
    'IC_Asp_Bend_trans_trans',
    'IC_Asp_Bridge',
    'IC_Asp_Bridge_gminus_gminus',
    'IC_Asp_Bridge_gminus_gplus',
    'IC_Asp_Bridge_gminus_trans',
    'IC_Asp_Bridge_gplus_gminus',
    'IC_Asp_Bridge_gplus_gplus',
    'IC_Asp_Bridge_gplus_trans',
    'IC_Asp_Bridge_trans_gminus',
    'IC_Asp_Bridge_trans_gplus',
    'IC_Asp_Bridge_trans_trans',
    'IC_Asp_CisPeptide',
    'IC_Asp_CisPeptide_gminus_gminus',
    'IC_Asp_CisPeptide_gminus_gplus',
    'IC_Asp_CisPeptide_gminus_trans',
    'IC_Asp_CisPeptide_gplus_gminus',
    'IC_Asp_CisPeptide_gplus_gplus',
    'IC_Asp_CisPeptide_gplus_trans',
    'IC_Asp_CisPeptide_trans_gminus',
    'IC_Asp_CisPeptide_trans_gplus',
    'IC_Asp_CisPeptide_trans_trans',
    'IC_Asp_Coil',
    'IC_Asp_Coil_gminus_gminus',
    'IC_Asp_Coil_gminus_gplus',
    'IC_Asp_Coil_gminus_trans',
    'IC_Asp_Coil_gplus_gminus',
    'IC_Asp_Coil_gplus_gplus',
    'IC_Asp_Coil_gplus_trans',
    'IC_Asp_Coil_trans_gminus',
    'IC_Asp_Coil_trans_gplus',
    'IC_Asp_Coil_trans_trans',
    'IC_Asp_Helix310',
    'IC_Asp_Helix310_gminus_gminus',
    'IC_Asp_Helix310_gminus_gplus',
    'IC_Asp_Helix310_gminus_trans',
    'IC_Asp_Helix310_gplus_gminus',
    'IC_Asp_Helix310_gplus_gplus',
    'IC_Asp_Helix310_gplus_trans',
    'IC_Asp_Helix310_trans_gminus',
    'IC_Asp_Helix310_trans_gplus',
    'IC_Asp_Helix310_trans_trans',
    'IC_Asp_HelixAlpha',
    'IC_Asp_HelixAlpha_gminus_gminus',
    'IC_Asp_HelixAlpha_gminus_gplus',
    'IC_Asp_HelixAlpha_gminus_trans',
    'IC_Asp_HelixAlpha_gplus_gminus',
    'IC_Asp_HelixAlpha_gplus_gplus',
    'IC_Asp_HelixAlpha_gplus_trans',
    'IC_Asp_HelixAlpha_trans_gminus',
    'IC_Asp_HelixAlpha_trans_gplus',
    'IC_Asp_HelixAlpha_trans_trans',
    'IC_Asp_HelixPPII',
    'IC_Asp_HelixPPII_gminus_gminus',
    'IC_Asp_HelixPPII_gminus_gplus',
    'IC_Asp_HelixPPII_gminus_trans',
    'IC_Asp_HelixPPII_gplus_gminus',
    'IC_Asp_HelixPPII_gplus_gplus',
    'IC_Asp_HelixPPII_gplus_trans',
    'IC_Asp_HelixPPII_trans_gminus',
    'IC_Asp_HelixPPII_trans_gplus',
    'IC_Asp_HelixPPII_trans_trans',
    'IC_Asp_HelixPi',
    'IC_Asp_HelixPi_gminus_gminus',
    'IC_Asp_HelixPi_gminus_gplus',
    'IC_Asp_HelixPi_gminus_trans',
    'IC_Asp_HelixPi_gplus_gminus',
    'IC_Asp_HelixPi_gplus_gplus',
    'IC_Asp_HelixPi_gplus_trans',
    'IC_Asp_HelixPi_trans_gminus',
    'IC_Asp_HelixPi_trans_gplus',
    'IC_Asp_HelixPi_trans_trans',
    'IC_Asp_Strand',
    'IC_Asp_StrandAntiParallel',
    'IC_Asp_StrandAntiParallel_gminus_gminus',
    'IC_Asp_StrandAntiParallel_gminus_gplus',
    'IC_Asp_StrandAntiParallel_gminus_trans',
    'IC_Asp_StrandAntiParallel_gplus_gminus',
    'IC_Asp_StrandAntiParallel_gplus_gplus',
    'IC_Asp_StrandAntiParallel_gplus_trans',
    'IC_Asp_StrandAntiParallel_trans_gminus',
    'IC_Asp_StrandAntiParallel_trans_gplus',
    'IC_Asp_StrandAntiParallel_trans_trans',
    'IC_Asp_StrandParallel',
    'IC_Asp_StrandParallel_gminus_gminus',
    'IC_Asp_StrandParallel_gminus_gplus',
    'IC_Asp_StrandParallel_gminus_trans',
    'IC_Asp_StrandParallel_gplus_gminus',
    'IC_Asp_StrandParallel_gplus_gplus',
    'IC_Asp_StrandParallel_gplus_trans',
    'IC_Asp_StrandParallel_trans_gminus',
    'IC_Asp_StrandParallel_trans_gplus',
    'IC_Asp_StrandParallel_trans_trans',
    'IC_Asp_Strand_gminus_gminus',
    'IC_Asp_Strand_gminus_gplus',
    'IC_Asp_Strand_gminus_trans',
    'IC_Asp_Strand_gplus_gminus',
    'IC_Asp_Strand_gplus_gplus',
    'IC_Asp_Strand_gplus_trans',
    'IC_Asp_Strand_trans_gminus',
    'IC_Asp_Strand_trans_gplus',
    'IC_Asp_Strand_trans_trans',
    'IC_Asp_Turn',
    'IC_Asp_Turn_gminus_gminus',
    'IC_Asp_Turn_gminus_gplus',
    'IC_Asp_Turn_gminus_trans',
    'IC_Asp_Turn_gplus_gminus',
    'IC_Asp_Turn_gplus_gplus',
    'IC_Asp_Turn_gplus_trans',
    'IC_Asp_Turn_trans_gminus',
    'IC_Asp_Turn_trans_gplus',
    'IC_Asp_Turn_trans_trans',
]
