# -*- coding: utf-8 -*-
"""Ideal {secondary-structure x rotamer} {names[aa]} {var} ``InternalCoord`` templates.

Each template is the ideal heavy-atom topology for {aa} anchored at
{{N, CA, C}}, with the side chain grown to ``rotamer`` and the
``phi``/``psi``/``omega`` of ``ss`` carried on the instance.  See
:mod:`biorazer.database.internal_coord_template` and
:mod:`._builder` for the build rule and the per-conformer caveat.
"""

from biorazer.database.internal_coord_template._builder import make_residue_templates

RESN = "TYR"
#: {ss: {rotamer: :class:`~biorazer.structure.objects.InternalCoord`}}
TYR_TEMPLATES = make_residue_templates(RESN)

# === GENERATED IC_* named templates - do not edit below ===
# Regenerate with: /opt/envs/BioRazer/bin/python scripts/generate_internal_coord_template_named.py

IC_Tyr_HelixAlpha = TYR_TEMPLATES['alpha-helix']['canonical']
IC_Tyr_HelixAlpha_gminus_gminus = TYR_TEMPLATES['alpha-helix']['g-/g-']
IC_Tyr_HelixAlpha_gminus_trans = TYR_TEMPLATES['alpha-helix']['g-/t']
IC_Tyr_HelixAlpha_gminus_gplus = TYR_TEMPLATES['alpha-helix']['g-/g+']
IC_Tyr_HelixAlpha_trans_gminus = TYR_TEMPLATES['alpha-helix']['t/g-']
IC_Tyr_HelixAlpha_trans_trans = TYR_TEMPLATES['alpha-helix']['t/t']
IC_Tyr_HelixAlpha_trans_gplus = TYR_TEMPLATES['alpha-helix']['t/g+']
IC_Tyr_HelixAlpha_gplus_gminus = TYR_TEMPLATES['alpha-helix']['g+/g-']
IC_Tyr_HelixAlpha_gplus_trans = TYR_TEMPLATES['alpha-helix']['g+/t']
IC_Tyr_HelixAlpha_gplus_gplus = TYR_TEMPLATES['alpha-helix']['g+/g+']
IC_Tyr_Helix310 = TYR_TEMPLATES['3-10-helix']['canonical']
IC_Tyr_Helix310_gminus_gminus = TYR_TEMPLATES['3-10-helix']['g-/g-']
IC_Tyr_Helix310_gminus_trans = TYR_TEMPLATES['3-10-helix']['g-/t']
IC_Tyr_Helix310_gminus_gplus = TYR_TEMPLATES['3-10-helix']['g-/g+']
IC_Tyr_Helix310_trans_gminus = TYR_TEMPLATES['3-10-helix']['t/g-']
IC_Tyr_Helix310_trans_trans = TYR_TEMPLATES['3-10-helix']['t/t']
IC_Tyr_Helix310_trans_gplus = TYR_TEMPLATES['3-10-helix']['t/g+']
IC_Tyr_Helix310_gplus_gminus = TYR_TEMPLATES['3-10-helix']['g+/g-']
IC_Tyr_Helix310_gplus_trans = TYR_TEMPLATES['3-10-helix']['g+/t']
IC_Tyr_Helix310_gplus_gplus = TYR_TEMPLATES['3-10-helix']['g+/g+']
IC_Tyr_HelixPi = TYR_TEMPLATES['pi-helix']['canonical']
IC_Tyr_HelixPi_gminus_gminus = TYR_TEMPLATES['pi-helix']['g-/g-']
IC_Tyr_HelixPi_gminus_trans = TYR_TEMPLATES['pi-helix']['g-/t']
IC_Tyr_HelixPi_gminus_gplus = TYR_TEMPLATES['pi-helix']['g-/g+']
IC_Tyr_HelixPi_trans_gminus = TYR_TEMPLATES['pi-helix']['t/g-']
IC_Tyr_HelixPi_trans_trans = TYR_TEMPLATES['pi-helix']['t/t']
IC_Tyr_HelixPi_trans_gplus = TYR_TEMPLATES['pi-helix']['t/g+']
IC_Tyr_HelixPi_gplus_gminus = TYR_TEMPLATES['pi-helix']['g+/g-']
IC_Tyr_HelixPi_gplus_trans = TYR_TEMPLATES['pi-helix']['g+/t']
IC_Tyr_HelixPi_gplus_gplus = TYR_TEMPLATES['pi-helix']['g+/g+']
IC_Tyr_HelixPPII = TYR_TEMPLATES['polyproline-II']['canonical']
IC_Tyr_HelixPPII_gminus_gminus = TYR_TEMPLATES['polyproline-II']['g-/g-']
IC_Tyr_HelixPPII_gminus_trans = TYR_TEMPLATES['polyproline-II']['g-/t']
IC_Tyr_HelixPPII_gminus_gplus = TYR_TEMPLATES['polyproline-II']['g-/g+']
IC_Tyr_HelixPPII_trans_gminus = TYR_TEMPLATES['polyproline-II']['t/g-']
IC_Tyr_HelixPPII_trans_trans = TYR_TEMPLATES['polyproline-II']['t/t']
IC_Tyr_HelixPPII_trans_gplus = TYR_TEMPLATES['polyproline-II']['t/g+']
IC_Tyr_HelixPPII_gplus_gminus = TYR_TEMPLATES['polyproline-II']['g+/g-']
IC_Tyr_HelixPPII_gplus_trans = TYR_TEMPLATES['polyproline-II']['g+/t']
IC_Tyr_HelixPPII_gplus_gplus = TYR_TEMPLATES['polyproline-II']['g+/g+']
IC_Tyr_Strand = TYR_TEMPLATES['beta-strand']['canonical']
IC_Tyr_Strand_gminus_gminus = TYR_TEMPLATES['beta-strand']['g-/g-']
IC_Tyr_Strand_gminus_trans = TYR_TEMPLATES['beta-strand']['g-/t']
IC_Tyr_Strand_gminus_gplus = TYR_TEMPLATES['beta-strand']['g-/g+']
IC_Tyr_Strand_trans_gminus = TYR_TEMPLATES['beta-strand']['t/g-']
IC_Tyr_Strand_trans_trans = TYR_TEMPLATES['beta-strand']['t/t']
IC_Tyr_Strand_trans_gplus = TYR_TEMPLATES['beta-strand']['t/g+']
IC_Tyr_Strand_gplus_gminus = TYR_TEMPLATES['beta-strand']['g+/g-']
IC_Tyr_Strand_gplus_trans = TYR_TEMPLATES['beta-strand']['g+/t']
IC_Tyr_Strand_gplus_gplus = TYR_TEMPLATES['beta-strand']['g+/g+']
IC_Tyr_StrandParallel = TYR_TEMPLATES['parallel-beta-strand']['canonical']
IC_Tyr_StrandParallel_gminus_gminus = TYR_TEMPLATES['parallel-beta-strand']['g-/g-']
IC_Tyr_StrandParallel_gminus_trans = TYR_TEMPLATES['parallel-beta-strand']['g-/t']
IC_Tyr_StrandParallel_gminus_gplus = TYR_TEMPLATES['parallel-beta-strand']['g-/g+']
IC_Tyr_StrandParallel_trans_gminus = TYR_TEMPLATES['parallel-beta-strand']['t/g-']
IC_Tyr_StrandParallel_trans_trans = TYR_TEMPLATES['parallel-beta-strand']['t/t']
IC_Tyr_StrandParallel_trans_gplus = TYR_TEMPLATES['parallel-beta-strand']['t/g+']
IC_Tyr_StrandParallel_gplus_gminus = TYR_TEMPLATES['parallel-beta-strand']['g+/g-']
IC_Tyr_StrandParallel_gplus_trans = TYR_TEMPLATES['parallel-beta-strand']['g+/t']
IC_Tyr_StrandParallel_gplus_gplus = TYR_TEMPLATES['parallel-beta-strand']['g+/g+']
IC_Tyr_StrandAntiParallel = TYR_TEMPLATES['antiparallel-beta-strand']['canonical']
IC_Tyr_StrandAntiParallel_gminus_gminus = TYR_TEMPLATES['antiparallel-beta-strand']['g-/g-']
IC_Tyr_StrandAntiParallel_gminus_trans = TYR_TEMPLATES['antiparallel-beta-strand']['g-/t']
IC_Tyr_StrandAntiParallel_gminus_gplus = TYR_TEMPLATES['antiparallel-beta-strand']['g-/g+']
IC_Tyr_StrandAntiParallel_trans_gminus = TYR_TEMPLATES['antiparallel-beta-strand']['t/g-']
IC_Tyr_StrandAntiParallel_trans_trans = TYR_TEMPLATES['antiparallel-beta-strand']['t/t']
IC_Tyr_StrandAntiParallel_trans_gplus = TYR_TEMPLATES['antiparallel-beta-strand']['t/g+']
IC_Tyr_StrandAntiParallel_gplus_gminus = TYR_TEMPLATES['antiparallel-beta-strand']['g+/g-']
IC_Tyr_StrandAntiParallel_gplus_trans = TYR_TEMPLATES['antiparallel-beta-strand']['g+/t']
IC_Tyr_StrandAntiParallel_gplus_gplus = TYR_TEMPLATES['antiparallel-beta-strand']['g+/g+']
IC_Tyr_Bridge = TYR_TEMPLATES['beta-bridge']['canonical']
IC_Tyr_Bridge_gminus_gminus = TYR_TEMPLATES['beta-bridge']['g-/g-']
IC_Tyr_Bridge_gminus_trans = TYR_TEMPLATES['beta-bridge']['g-/t']
IC_Tyr_Bridge_gminus_gplus = TYR_TEMPLATES['beta-bridge']['g-/g+']
IC_Tyr_Bridge_trans_gminus = TYR_TEMPLATES['beta-bridge']['t/g-']
IC_Tyr_Bridge_trans_trans = TYR_TEMPLATES['beta-bridge']['t/t']
IC_Tyr_Bridge_trans_gplus = TYR_TEMPLATES['beta-bridge']['t/g+']
IC_Tyr_Bridge_gplus_gminus = TYR_TEMPLATES['beta-bridge']['g+/g-']
IC_Tyr_Bridge_gplus_trans = TYR_TEMPLATES['beta-bridge']['g+/t']
IC_Tyr_Bridge_gplus_gplus = TYR_TEMPLATES['beta-bridge']['g+/g+']
IC_Tyr_Turn = TYR_TEMPLATES['turn']['canonical']
IC_Tyr_Turn_gminus_gminus = TYR_TEMPLATES['turn']['g-/g-']
IC_Tyr_Turn_gminus_trans = TYR_TEMPLATES['turn']['g-/t']
IC_Tyr_Turn_gminus_gplus = TYR_TEMPLATES['turn']['g-/g+']
IC_Tyr_Turn_trans_gminus = TYR_TEMPLATES['turn']['t/g-']
IC_Tyr_Turn_trans_trans = TYR_TEMPLATES['turn']['t/t']
IC_Tyr_Turn_trans_gplus = TYR_TEMPLATES['turn']['t/g+']
IC_Tyr_Turn_gplus_gminus = TYR_TEMPLATES['turn']['g+/g-']
IC_Tyr_Turn_gplus_trans = TYR_TEMPLATES['turn']['g+/t']
IC_Tyr_Turn_gplus_gplus = TYR_TEMPLATES['turn']['g+/g+']
IC_Tyr_Bend = TYR_TEMPLATES['bend']['canonical']
IC_Tyr_Bend_gminus_gminus = TYR_TEMPLATES['bend']['g-/g-']
IC_Tyr_Bend_gminus_trans = TYR_TEMPLATES['bend']['g-/t']
IC_Tyr_Bend_gminus_gplus = TYR_TEMPLATES['bend']['g-/g+']
IC_Tyr_Bend_trans_gminus = TYR_TEMPLATES['bend']['t/g-']
IC_Tyr_Bend_trans_trans = TYR_TEMPLATES['bend']['t/t']
IC_Tyr_Bend_trans_gplus = TYR_TEMPLATES['bend']['t/g+']
IC_Tyr_Bend_gplus_gminus = TYR_TEMPLATES['bend']['g+/g-']
IC_Tyr_Bend_gplus_trans = TYR_TEMPLATES['bend']['g+/t']
IC_Tyr_Bend_gplus_gplus = TYR_TEMPLATES['bend']['g+/g+']
IC_Tyr_Coil = TYR_TEMPLATES['coil']['canonical']
IC_Tyr_Coil_gminus_gminus = TYR_TEMPLATES['coil']['g-/g-']
IC_Tyr_Coil_gminus_trans = TYR_TEMPLATES['coil']['g-/t']
IC_Tyr_Coil_gminus_gplus = TYR_TEMPLATES['coil']['g-/g+']
IC_Tyr_Coil_trans_gminus = TYR_TEMPLATES['coil']['t/g-']
IC_Tyr_Coil_trans_trans = TYR_TEMPLATES['coil']['t/t']
IC_Tyr_Coil_trans_gplus = TYR_TEMPLATES['coil']['t/g+']
IC_Tyr_Coil_gplus_gminus = TYR_TEMPLATES['coil']['g+/g-']
IC_Tyr_Coil_gplus_trans = TYR_TEMPLATES['coil']['g+/t']
IC_Tyr_Coil_gplus_gplus = TYR_TEMPLATES['coil']['g+/g+']
IC_Tyr_CisPeptide = TYR_TEMPLATES['cis-peptide-bond']['canonical']
IC_Tyr_CisPeptide_gminus_gminus = TYR_TEMPLATES['cis-peptide-bond']['g-/g-']
IC_Tyr_CisPeptide_gminus_trans = TYR_TEMPLATES['cis-peptide-bond']['g-/t']
IC_Tyr_CisPeptide_gminus_gplus = TYR_TEMPLATES['cis-peptide-bond']['g-/g+']
IC_Tyr_CisPeptide_trans_gminus = TYR_TEMPLATES['cis-peptide-bond']['t/g-']
IC_Tyr_CisPeptide_trans_trans = TYR_TEMPLATES['cis-peptide-bond']['t/t']
IC_Tyr_CisPeptide_trans_gplus = TYR_TEMPLATES['cis-peptide-bond']['t/g+']
IC_Tyr_CisPeptide_gplus_gminus = TYR_TEMPLATES['cis-peptide-bond']['g+/g-']
IC_Tyr_CisPeptide_gplus_trans = TYR_TEMPLATES['cis-peptide-bond']['g+/t']
IC_Tyr_CisPeptide_gplus_gplus = TYR_TEMPLATES['cis-peptide-bond']['g+/g+']

__all__ = [
    'IC_Tyr_Bend',
    'IC_Tyr_Bend_gminus_gminus',
    'IC_Tyr_Bend_gminus_gplus',
    'IC_Tyr_Bend_gminus_trans',
    'IC_Tyr_Bend_gplus_gminus',
    'IC_Tyr_Bend_gplus_gplus',
    'IC_Tyr_Bend_gplus_trans',
    'IC_Tyr_Bend_trans_gminus',
    'IC_Tyr_Bend_trans_gplus',
    'IC_Tyr_Bend_trans_trans',
    'IC_Tyr_Bridge',
    'IC_Tyr_Bridge_gminus_gminus',
    'IC_Tyr_Bridge_gminus_gplus',
    'IC_Tyr_Bridge_gminus_trans',
    'IC_Tyr_Bridge_gplus_gminus',
    'IC_Tyr_Bridge_gplus_gplus',
    'IC_Tyr_Bridge_gplus_trans',
    'IC_Tyr_Bridge_trans_gminus',
    'IC_Tyr_Bridge_trans_gplus',
    'IC_Tyr_Bridge_trans_trans',
    'IC_Tyr_CisPeptide',
    'IC_Tyr_CisPeptide_gminus_gminus',
    'IC_Tyr_CisPeptide_gminus_gplus',
    'IC_Tyr_CisPeptide_gminus_trans',
    'IC_Tyr_CisPeptide_gplus_gminus',
    'IC_Tyr_CisPeptide_gplus_gplus',
    'IC_Tyr_CisPeptide_gplus_trans',
    'IC_Tyr_CisPeptide_trans_gminus',
    'IC_Tyr_CisPeptide_trans_gplus',
    'IC_Tyr_CisPeptide_trans_trans',
    'IC_Tyr_Coil',
    'IC_Tyr_Coil_gminus_gminus',
    'IC_Tyr_Coil_gminus_gplus',
    'IC_Tyr_Coil_gminus_trans',
    'IC_Tyr_Coil_gplus_gminus',
    'IC_Tyr_Coil_gplus_gplus',
    'IC_Tyr_Coil_gplus_trans',
    'IC_Tyr_Coil_trans_gminus',
    'IC_Tyr_Coil_trans_gplus',
    'IC_Tyr_Coil_trans_trans',
    'IC_Tyr_Helix310',
    'IC_Tyr_Helix310_gminus_gminus',
    'IC_Tyr_Helix310_gminus_gplus',
    'IC_Tyr_Helix310_gminus_trans',
    'IC_Tyr_Helix310_gplus_gminus',
    'IC_Tyr_Helix310_gplus_gplus',
    'IC_Tyr_Helix310_gplus_trans',
    'IC_Tyr_Helix310_trans_gminus',
    'IC_Tyr_Helix310_trans_gplus',
    'IC_Tyr_Helix310_trans_trans',
    'IC_Tyr_HelixAlpha',
    'IC_Tyr_HelixAlpha_gminus_gminus',
    'IC_Tyr_HelixAlpha_gminus_gplus',
    'IC_Tyr_HelixAlpha_gminus_trans',
    'IC_Tyr_HelixAlpha_gplus_gminus',
    'IC_Tyr_HelixAlpha_gplus_gplus',
    'IC_Tyr_HelixAlpha_gplus_trans',
    'IC_Tyr_HelixAlpha_trans_gminus',
    'IC_Tyr_HelixAlpha_trans_gplus',
    'IC_Tyr_HelixAlpha_trans_trans',
    'IC_Tyr_HelixPPII',
    'IC_Tyr_HelixPPII_gminus_gminus',
    'IC_Tyr_HelixPPII_gminus_gplus',
    'IC_Tyr_HelixPPII_gminus_trans',
    'IC_Tyr_HelixPPII_gplus_gminus',
    'IC_Tyr_HelixPPII_gplus_gplus',
    'IC_Tyr_HelixPPII_gplus_trans',
    'IC_Tyr_HelixPPII_trans_gminus',
    'IC_Tyr_HelixPPII_trans_gplus',
    'IC_Tyr_HelixPPII_trans_trans',
    'IC_Tyr_HelixPi',
    'IC_Tyr_HelixPi_gminus_gminus',
    'IC_Tyr_HelixPi_gminus_gplus',
    'IC_Tyr_HelixPi_gminus_trans',
    'IC_Tyr_HelixPi_gplus_gminus',
    'IC_Tyr_HelixPi_gplus_gplus',
    'IC_Tyr_HelixPi_gplus_trans',
    'IC_Tyr_HelixPi_trans_gminus',
    'IC_Tyr_HelixPi_trans_gplus',
    'IC_Tyr_HelixPi_trans_trans',
    'IC_Tyr_Strand',
    'IC_Tyr_StrandAntiParallel',
    'IC_Tyr_StrandAntiParallel_gminus_gminus',
    'IC_Tyr_StrandAntiParallel_gminus_gplus',
    'IC_Tyr_StrandAntiParallel_gminus_trans',
    'IC_Tyr_StrandAntiParallel_gplus_gminus',
    'IC_Tyr_StrandAntiParallel_gplus_gplus',
    'IC_Tyr_StrandAntiParallel_gplus_trans',
    'IC_Tyr_StrandAntiParallel_trans_gminus',
    'IC_Tyr_StrandAntiParallel_trans_gplus',
    'IC_Tyr_StrandAntiParallel_trans_trans',
    'IC_Tyr_StrandParallel',
    'IC_Tyr_StrandParallel_gminus_gminus',
    'IC_Tyr_StrandParallel_gminus_gplus',
    'IC_Tyr_StrandParallel_gminus_trans',
    'IC_Tyr_StrandParallel_gplus_gminus',
    'IC_Tyr_StrandParallel_gplus_gplus',
    'IC_Tyr_StrandParallel_gplus_trans',
    'IC_Tyr_StrandParallel_trans_gminus',
    'IC_Tyr_StrandParallel_trans_gplus',
    'IC_Tyr_StrandParallel_trans_trans',
    'IC_Tyr_Strand_gminus_gminus',
    'IC_Tyr_Strand_gminus_gplus',
    'IC_Tyr_Strand_gminus_trans',
    'IC_Tyr_Strand_gplus_gminus',
    'IC_Tyr_Strand_gplus_gplus',
    'IC_Tyr_Strand_gplus_trans',
    'IC_Tyr_Strand_trans_gminus',
    'IC_Tyr_Strand_trans_gplus',
    'IC_Tyr_Strand_trans_trans',
    'IC_Tyr_Turn',
    'IC_Tyr_Turn_gminus_gminus',
    'IC_Tyr_Turn_gminus_gplus',
    'IC_Tyr_Turn_gminus_trans',
    'IC_Tyr_Turn_gplus_gminus',
    'IC_Tyr_Turn_gplus_gplus',
    'IC_Tyr_Turn_gplus_trans',
    'IC_Tyr_Turn_trans_gminus',
    'IC_Tyr_Turn_trans_gplus',
    'IC_Tyr_Turn_trans_trans',
]
