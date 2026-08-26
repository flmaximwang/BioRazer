# -*- coding: utf-8 -*-
"""Ideal {secondary-structure x rotamer} {names[aa]} {var} ``InternalCoord`` templates.

Each template is the ideal heavy-atom topology for {aa} anchored at
{{N, CA, C}}, with the side chain grown to ``rotamer`` and the
``phi``/``psi``/``omega`` of ``ss`` carried on the instance.  See
:mod:`biorazer.database.internal_coord_template` and
:mod:`._builder` for the build rule and the per-conformer caveat.
"""

from biorazer.database.internal_coord_template._builder import make_residue_templates

RESN = "TRP"
#: {ss: {rotamer: :class:`~biorazer.structure.objects.InternalCoord`}}
TRP_TEMPLATES = make_residue_templates(RESN)

# === GENERATED IC_* named templates - do not edit below ===
# Regenerate with: /opt/envs/BioRazer/bin/python scripts/generate_internal_coord_template_named.py

IC_Trp_HelixAlpha = TRP_TEMPLATES['alpha-helix']['canonical']
IC_Trp_HelixAlpha_gminus_gminus = TRP_TEMPLATES['alpha-helix']['g-/g-']
IC_Trp_HelixAlpha_gminus_trans = TRP_TEMPLATES['alpha-helix']['g-/t']
IC_Trp_HelixAlpha_gminus_gplus = TRP_TEMPLATES['alpha-helix']['g-/g+']
IC_Trp_HelixAlpha_trans_gminus = TRP_TEMPLATES['alpha-helix']['t/g-']
IC_Trp_HelixAlpha_trans_trans = TRP_TEMPLATES['alpha-helix']['t/t']
IC_Trp_HelixAlpha_trans_gplus = TRP_TEMPLATES['alpha-helix']['t/g+']
IC_Trp_HelixAlpha_gplus_gminus = TRP_TEMPLATES['alpha-helix']['g+/g-']
IC_Trp_HelixAlpha_gplus_trans = TRP_TEMPLATES['alpha-helix']['g+/t']
IC_Trp_HelixAlpha_gplus_gplus = TRP_TEMPLATES['alpha-helix']['g+/g+']
IC_Trp_Helix310 = TRP_TEMPLATES['3-10-helix']['canonical']
IC_Trp_Helix310_gminus_gminus = TRP_TEMPLATES['3-10-helix']['g-/g-']
IC_Trp_Helix310_gminus_trans = TRP_TEMPLATES['3-10-helix']['g-/t']
IC_Trp_Helix310_gminus_gplus = TRP_TEMPLATES['3-10-helix']['g-/g+']
IC_Trp_Helix310_trans_gminus = TRP_TEMPLATES['3-10-helix']['t/g-']
IC_Trp_Helix310_trans_trans = TRP_TEMPLATES['3-10-helix']['t/t']
IC_Trp_Helix310_trans_gplus = TRP_TEMPLATES['3-10-helix']['t/g+']
IC_Trp_Helix310_gplus_gminus = TRP_TEMPLATES['3-10-helix']['g+/g-']
IC_Trp_Helix310_gplus_trans = TRP_TEMPLATES['3-10-helix']['g+/t']
IC_Trp_Helix310_gplus_gplus = TRP_TEMPLATES['3-10-helix']['g+/g+']
IC_Trp_HelixPi = TRP_TEMPLATES['pi-helix']['canonical']
IC_Trp_HelixPi_gminus_gminus = TRP_TEMPLATES['pi-helix']['g-/g-']
IC_Trp_HelixPi_gminus_trans = TRP_TEMPLATES['pi-helix']['g-/t']
IC_Trp_HelixPi_gminus_gplus = TRP_TEMPLATES['pi-helix']['g-/g+']
IC_Trp_HelixPi_trans_gminus = TRP_TEMPLATES['pi-helix']['t/g-']
IC_Trp_HelixPi_trans_trans = TRP_TEMPLATES['pi-helix']['t/t']
IC_Trp_HelixPi_trans_gplus = TRP_TEMPLATES['pi-helix']['t/g+']
IC_Trp_HelixPi_gplus_gminus = TRP_TEMPLATES['pi-helix']['g+/g-']
IC_Trp_HelixPi_gplus_trans = TRP_TEMPLATES['pi-helix']['g+/t']
IC_Trp_HelixPi_gplus_gplus = TRP_TEMPLATES['pi-helix']['g+/g+']
IC_Trp_HelixPPII = TRP_TEMPLATES['polyproline-II']['canonical']
IC_Trp_HelixPPII_gminus_gminus = TRP_TEMPLATES['polyproline-II']['g-/g-']
IC_Trp_HelixPPII_gminus_trans = TRP_TEMPLATES['polyproline-II']['g-/t']
IC_Trp_HelixPPII_gminus_gplus = TRP_TEMPLATES['polyproline-II']['g-/g+']
IC_Trp_HelixPPII_trans_gminus = TRP_TEMPLATES['polyproline-II']['t/g-']
IC_Trp_HelixPPII_trans_trans = TRP_TEMPLATES['polyproline-II']['t/t']
IC_Trp_HelixPPII_trans_gplus = TRP_TEMPLATES['polyproline-II']['t/g+']
IC_Trp_HelixPPII_gplus_gminus = TRP_TEMPLATES['polyproline-II']['g+/g-']
IC_Trp_HelixPPII_gplus_trans = TRP_TEMPLATES['polyproline-II']['g+/t']
IC_Trp_HelixPPII_gplus_gplus = TRP_TEMPLATES['polyproline-II']['g+/g+']
IC_Trp_Strand = TRP_TEMPLATES['beta-strand']['canonical']
IC_Trp_Strand_gminus_gminus = TRP_TEMPLATES['beta-strand']['g-/g-']
IC_Trp_Strand_gminus_trans = TRP_TEMPLATES['beta-strand']['g-/t']
IC_Trp_Strand_gminus_gplus = TRP_TEMPLATES['beta-strand']['g-/g+']
IC_Trp_Strand_trans_gminus = TRP_TEMPLATES['beta-strand']['t/g-']
IC_Trp_Strand_trans_trans = TRP_TEMPLATES['beta-strand']['t/t']
IC_Trp_Strand_trans_gplus = TRP_TEMPLATES['beta-strand']['t/g+']
IC_Trp_Strand_gplus_gminus = TRP_TEMPLATES['beta-strand']['g+/g-']
IC_Trp_Strand_gplus_trans = TRP_TEMPLATES['beta-strand']['g+/t']
IC_Trp_Strand_gplus_gplus = TRP_TEMPLATES['beta-strand']['g+/g+']
IC_Trp_StrandParallel = TRP_TEMPLATES['parallel-beta-strand']['canonical']
IC_Trp_StrandParallel_gminus_gminus = TRP_TEMPLATES['parallel-beta-strand']['g-/g-']
IC_Trp_StrandParallel_gminus_trans = TRP_TEMPLATES['parallel-beta-strand']['g-/t']
IC_Trp_StrandParallel_gminus_gplus = TRP_TEMPLATES['parallel-beta-strand']['g-/g+']
IC_Trp_StrandParallel_trans_gminus = TRP_TEMPLATES['parallel-beta-strand']['t/g-']
IC_Trp_StrandParallel_trans_trans = TRP_TEMPLATES['parallel-beta-strand']['t/t']
IC_Trp_StrandParallel_trans_gplus = TRP_TEMPLATES['parallel-beta-strand']['t/g+']
IC_Trp_StrandParallel_gplus_gminus = TRP_TEMPLATES['parallel-beta-strand']['g+/g-']
IC_Trp_StrandParallel_gplus_trans = TRP_TEMPLATES['parallel-beta-strand']['g+/t']
IC_Trp_StrandParallel_gplus_gplus = TRP_TEMPLATES['parallel-beta-strand']['g+/g+']
IC_Trp_StrandAntiParallel = TRP_TEMPLATES['antiparallel-beta-strand']['canonical']
IC_Trp_StrandAntiParallel_gminus_gminus = TRP_TEMPLATES['antiparallel-beta-strand']['g-/g-']
IC_Trp_StrandAntiParallel_gminus_trans = TRP_TEMPLATES['antiparallel-beta-strand']['g-/t']
IC_Trp_StrandAntiParallel_gminus_gplus = TRP_TEMPLATES['antiparallel-beta-strand']['g-/g+']
IC_Trp_StrandAntiParallel_trans_gminus = TRP_TEMPLATES['antiparallel-beta-strand']['t/g-']
IC_Trp_StrandAntiParallel_trans_trans = TRP_TEMPLATES['antiparallel-beta-strand']['t/t']
IC_Trp_StrandAntiParallel_trans_gplus = TRP_TEMPLATES['antiparallel-beta-strand']['t/g+']
IC_Trp_StrandAntiParallel_gplus_gminus = TRP_TEMPLATES['antiparallel-beta-strand']['g+/g-']
IC_Trp_StrandAntiParallel_gplus_trans = TRP_TEMPLATES['antiparallel-beta-strand']['g+/t']
IC_Trp_StrandAntiParallel_gplus_gplus = TRP_TEMPLATES['antiparallel-beta-strand']['g+/g+']
IC_Trp_Bridge = TRP_TEMPLATES['beta-bridge']['canonical']
IC_Trp_Bridge_gminus_gminus = TRP_TEMPLATES['beta-bridge']['g-/g-']
IC_Trp_Bridge_gminus_trans = TRP_TEMPLATES['beta-bridge']['g-/t']
IC_Trp_Bridge_gminus_gplus = TRP_TEMPLATES['beta-bridge']['g-/g+']
IC_Trp_Bridge_trans_gminus = TRP_TEMPLATES['beta-bridge']['t/g-']
IC_Trp_Bridge_trans_trans = TRP_TEMPLATES['beta-bridge']['t/t']
IC_Trp_Bridge_trans_gplus = TRP_TEMPLATES['beta-bridge']['t/g+']
IC_Trp_Bridge_gplus_gminus = TRP_TEMPLATES['beta-bridge']['g+/g-']
IC_Trp_Bridge_gplus_trans = TRP_TEMPLATES['beta-bridge']['g+/t']
IC_Trp_Bridge_gplus_gplus = TRP_TEMPLATES['beta-bridge']['g+/g+']
IC_Trp_Turn = TRP_TEMPLATES['turn']['canonical']
IC_Trp_Turn_gminus_gminus = TRP_TEMPLATES['turn']['g-/g-']
IC_Trp_Turn_gminus_trans = TRP_TEMPLATES['turn']['g-/t']
IC_Trp_Turn_gminus_gplus = TRP_TEMPLATES['turn']['g-/g+']
IC_Trp_Turn_trans_gminus = TRP_TEMPLATES['turn']['t/g-']
IC_Trp_Turn_trans_trans = TRP_TEMPLATES['turn']['t/t']
IC_Trp_Turn_trans_gplus = TRP_TEMPLATES['turn']['t/g+']
IC_Trp_Turn_gplus_gminus = TRP_TEMPLATES['turn']['g+/g-']
IC_Trp_Turn_gplus_trans = TRP_TEMPLATES['turn']['g+/t']
IC_Trp_Turn_gplus_gplus = TRP_TEMPLATES['turn']['g+/g+']
IC_Trp_Bend = TRP_TEMPLATES['bend']['canonical']
IC_Trp_Bend_gminus_gminus = TRP_TEMPLATES['bend']['g-/g-']
IC_Trp_Bend_gminus_trans = TRP_TEMPLATES['bend']['g-/t']
IC_Trp_Bend_gminus_gplus = TRP_TEMPLATES['bend']['g-/g+']
IC_Trp_Bend_trans_gminus = TRP_TEMPLATES['bend']['t/g-']
IC_Trp_Bend_trans_trans = TRP_TEMPLATES['bend']['t/t']
IC_Trp_Bend_trans_gplus = TRP_TEMPLATES['bend']['t/g+']
IC_Trp_Bend_gplus_gminus = TRP_TEMPLATES['bend']['g+/g-']
IC_Trp_Bend_gplus_trans = TRP_TEMPLATES['bend']['g+/t']
IC_Trp_Bend_gplus_gplus = TRP_TEMPLATES['bend']['g+/g+']
IC_Trp_Coil = TRP_TEMPLATES['coil']['canonical']
IC_Trp_Coil_gminus_gminus = TRP_TEMPLATES['coil']['g-/g-']
IC_Trp_Coil_gminus_trans = TRP_TEMPLATES['coil']['g-/t']
IC_Trp_Coil_gminus_gplus = TRP_TEMPLATES['coil']['g-/g+']
IC_Trp_Coil_trans_gminus = TRP_TEMPLATES['coil']['t/g-']
IC_Trp_Coil_trans_trans = TRP_TEMPLATES['coil']['t/t']
IC_Trp_Coil_trans_gplus = TRP_TEMPLATES['coil']['t/g+']
IC_Trp_Coil_gplus_gminus = TRP_TEMPLATES['coil']['g+/g-']
IC_Trp_Coil_gplus_trans = TRP_TEMPLATES['coil']['g+/t']
IC_Trp_Coil_gplus_gplus = TRP_TEMPLATES['coil']['g+/g+']
IC_Trp_CisPeptide = TRP_TEMPLATES['cis-peptide-bond']['canonical']
IC_Trp_CisPeptide_gminus_gminus = TRP_TEMPLATES['cis-peptide-bond']['g-/g-']
IC_Trp_CisPeptide_gminus_trans = TRP_TEMPLATES['cis-peptide-bond']['g-/t']
IC_Trp_CisPeptide_gminus_gplus = TRP_TEMPLATES['cis-peptide-bond']['g-/g+']
IC_Trp_CisPeptide_trans_gminus = TRP_TEMPLATES['cis-peptide-bond']['t/g-']
IC_Trp_CisPeptide_trans_trans = TRP_TEMPLATES['cis-peptide-bond']['t/t']
IC_Trp_CisPeptide_trans_gplus = TRP_TEMPLATES['cis-peptide-bond']['t/g+']
IC_Trp_CisPeptide_gplus_gminus = TRP_TEMPLATES['cis-peptide-bond']['g+/g-']
IC_Trp_CisPeptide_gplus_trans = TRP_TEMPLATES['cis-peptide-bond']['g+/t']
IC_Trp_CisPeptide_gplus_gplus = TRP_TEMPLATES['cis-peptide-bond']['g+/g+']

__all__ = [
    'IC_Trp_Bend',
    'IC_Trp_Bend_gminus_gminus',
    'IC_Trp_Bend_gminus_gplus',
    'IC_Trp_Bend_gminus_trans',
    'IC_Trp_Bend_gplus_gminus',
    'IC_Trp_Bend_gplus_gplus',
    'IC_Trp_Bend_gplus_trans',
    'IC_Trp_Bend_trans_gminus',
    'IC_Trp_Bend_trans_gplus',
    'IC_Trp_Bend_trans_trans',
    'IC_Trp_Bridge',
    'IC_Trp_Bridge_gminus_gminus',
    'IC_Trp_Bridge_gminus_gplus',
    'IC_Trp_Bridge_gminus_trans',
    'IC_Trp_Bridge_gplus_gminus',
    'IC_Trp_Bridge_gplus_gplus',
    'IC_Trp_Bridge_gplus_trans',
    'IC_Trp_Bridge_trans_gminus',
    'IC_Trp_Bridge_trans_gplus',
    'IC_Trp_Bridge_trans_trans',
    'IC_Trp_CisPeptide',
    'IC_Trp_CisPeptide_gminus_gminus',
    'IC_Trp_CisPeptide_gminus_gplus',
    'IC_Trp_CisPeptide_gminus_trans',
    'IC_Trp_CisPeptide_gplus_gminus',
    'IC_Trp_CisPeptide_gplus_gplus',
    'IC_Trp_CisPeptide_gplus_trans',
    'IC_Trp_CisPeptide_trans_gminus',
    'IC_Trp_CisPeptide_trans_gplus',
    'IC_Trp_CisPeptide_trans_trans',
    'IC_Trp_Coil',
    'IC_Trp_Coil_gminus_gminus',
    'IC_Trp_Coil_gminus_gplus',
    'IC_Trp_Coil_gminus_trans',
    'IC_Trp_Coil_gplus_gminus',
    'IC_Trp_Coil_gplus_gplus',
    'IC_Trp_Coil_gplus_trans',
    'IC_Trp_Coil_trans_gminus',
    'IC_Trp_Coil_trans_gplus',
    'IC_Trp_Coil_trans_trans',
    'IC_Trp_Helix310',
    'IC_Trp_Helix310_gminus_gminus',
    'IC_Trp_Helix310_gminus_gplus',
    'IC_Trp_Helix310_gminus_trans',
    'IC_Trp_Helix310_gplus_gminus',
    'IC_Trp_Helix310_gplus_gplus',
    'IC_Trp_Helix310_gplus_trans',
    'IC_Trp_Helix310_trans_gminus',
    'IC_Trp_Helix310_trans_gplus',
    'IC_Trp_Helix310_trans_trans',
    'IC_Trp_HelixAlpha',
    'IC_Trp_HelixAlpha_gminus_gminus',
    'IC_Trp_HelixAlpha_gminus_gplus',
    'IC_Trp_HelixAlpha_gminus_trans',
    'IC_Trp_HelixAlpha_gplus_gminus',
    'IC_Trp_HelixAlpha_gplus_gplus',
    'IC_Trp_HelixAlpha_gplus_trans',
    'IC_Trp_HelixAlpha_trans_gminus',
    'IC_Trp_HelixAlpha_trans_gplus',
    'IC_Trp_HelixAlpha_trans_trans',
    'IC_Trp_HelixPPII',
    'IC_Trp_HelixPPII_gminus_gminus',
    'IC_Trp_HelixPPII_gminus_gplus',
    'IC_Trp_HelixPPII_gminus_trans',
    'IC_Trp_HelixPPII_gplus_gminus',
    'IC_Trp_HelixPPII_gplus_gplus',
    'IC_Trp_HelixPPII_gplus_trans',
    'IC_Trp_HelixPPII_trans_gminus',
    'IC_Trp_HelixPPII_trans_gplus',
    'IC_Trp_HelixPPII_trans_trans',
    'IC_Trp_HelixPi',
    'IC_Trp_HelixPi_gminus_gminus',
    'IC_Trp_HelixPi_gminus_gplus',
    'IC_Trp_HelixPi_gminus_trans',
    'IC_Trp_HelixPi_gplus_gminus',
    'IC_Trp_HelixPi_gplus_gplus',
    'IC_Trp_HelixPi_gplus_trans',
    'IC_Trp_HelixPi_trans_gminus',
    'IC_Trp_HelixPi_trans_gplus',
    'IC_Trp_HelixPi_trans_trans',
    'IC_Trp_Strand',
    'IC_Trp_StrandAntiParallel',
    'IC_Trp_StrandAntiParallel_gminus_gminus',
    'IC_Trp_StrandAntiParallel_gminus_gplus',
    'IC_Trp_StrandAntiParallel_gminus_trans',
    'IC_Trp_StrandAntiParallel_gplus_gminus',
    'IC_Trp_StrandAntiParallel_gplus_gplus',
    'IC_Trp_StrandAntiParallel_gplus_trans',
    'IC_Trp_StrandAntiParallel_trans_gminus',
    'IC_Trp_StrandAntiParallel_trans_gplus',
    'IC_Trp_StrandAntiParallel_trans_trans',
    'IC_Trp_StrandParallel',
    'IC_Trp_StrandParallel_gminus_gminus',
    'IC_Trp_StrandParallel_gminus_gplus',
    'IC_Trp_StrandParallel_gminus_trans',
    'IC_Trp_StrandParallel_gplus_gminus',
    'IC_Trp_StrandParallel_gplus_gplus',
    'IC_Trp_StrandParallel_gplus_trans',
    'IC_Trp_StrandParallel_trans_gminus',
    'IC_Trp_StrandParallel_trans_gplus',
    'IC_Trp_StrandParallel_trans_trans',
    'IC_Trp_Strand_gminus_gminus',
    'IC_Trp_Strand_gminus_gplus',
    'IC_Trp_Strand_gminus_trans',
    'IC_Trp_Strand_gplus_gminus',
    'IC_Trp_Strand_gplus_gplus',
    'IC_Trp_Strand_gplus_trans',
    'IC_Trp_Strand_trans_gminus',
    'IC_Trp_Strand_trans_gplus',
    'IC_Trp_Strand_trans_trans',
    'IC_Trp_Turn',
    'IC_Trp_Turn_gminus_gminus',
    'IC_Trp_Turn_gminus_gplus',
    'IC_Trp_Turn_gminus_trans',
    'IC_Trp_Turn_gplus_gminus',
    'IC_Trp_Turn_gplus_gplus',
    'IC_Trp_Turn_gplus_trans',
    'IC_Trp_Turn_trans_gminus',
    'IC_Trp_Turn_trans_gplus',
    'IC_Trp_Turn_trans_trans',
]
