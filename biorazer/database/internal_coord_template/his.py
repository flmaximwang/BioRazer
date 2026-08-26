# -*- coding: utf-8 -*-
"""Ideal {secondary-structure x rotamer} {names[aa]} {var} ``InternalCoord`` templates.

Each template is the ideal heavy-atom topology for {aa} anchored at
{{N, CA, C}}, with the side chain grown to ``rotamer`` and the
``phi``/``psi``/``omega`` of ``ss`` carried on the instance.  See
:mod:`biorazer.database.internal_coord_template` and
:mod:`._builder` for the build rule and the per-conformer caveat.
"""

from biorazer.database.internal_coord_template._builder import make_residue_templates

RESN = "HIS"
#: {ss: {rotamer: :class:`~biorazer.structure.objects.InternalCoord`}}
HIS_TEMPLATES = make_residue_templates(RESN)

# === GENERATED IC_* named templates - do not edit below ===
# Regenerate with: /opt/envs/BioRazer/bin/python scripts/generate_internal_coord_template_named.py

IC_His_HelixAlpha = HIS_TEMPLATES['alpha-helix']['canonical']
IC_His_HelixAlpha_gminus_gminus = HIS_TEMPLATES['alpha-helix']['g-/g-']
IC_His_HelixAlpha_gminus_trans = HIS_TEMPLATES['alpha-helix']['g-/t']
IC_His_HelixAlpha_gminus_gplus = HIS_TEMPLATES['alpha-helix']['g-/g+']
IC_His_HelixAlpha_trans_gminus = HIS_TEMPLATES['alpha-helix']['t/g-']
IC_His_HelixAlpha_trans_trans = HIS_TEMPLATES['alpha-helix']['t/t']
IC_His_HelixAlpha_trans_gplus = HIS_TEMPLATES['alpha-helix']['t/g+']
IC_His_HelixAlpha_gplus_gminus = HIS_TEMPLATES['alpha-helix']['g+/g-']
IC_His_HelixAlpha_gplus_trans = HIS_TEMPLATES['alpha-helix']['g+/t']
IC_His_HelixAlpha_gplus_gplus = HIS_TEMPLATES['alpha-helix']['g+/g+']
IC_His_Helix310 = HIS_TEMPLATES['3-10-helix']['canonical']
IC_His_Helix310_gminus_gminus = HIS_TEMPLATES['3-10-helix']['g-/g-']
IC_His_Helix310_gminus_trans = HIS_TEMPLATES['3-10-helix']['g-/t']
IC_His_Helix310_gminus_gplus = HIS_TEMPLATES['3-10-helix']['g-/g+']
IC_His_Helix310_trans_gminus = HIS_TEMPLATES['3-10-helix']['t/g-']
IC_His_Helix310_trans_trans = HIS_TEMPLATES['3-10-helix']['t/t']
IC_His_Helix310_trans_gplus = HIS_TEMPLATES['3-10-helix']['t/g+']
IC_His_Helix310_gplus_gminus = HIS_TEMPLATES['3-10-helix']['g+/g-']
IC_His_Helix310_gplus_trans = HIS_TEMPLATES['3-10-helix']['g+/t']
IC_His_Helix310_gplus_gplus = HIS_TEMPLATES['3-10-helix']['g+/g+']
IC_His_HelixPi = HIS_TEMPLATES['pi-helix']['canonical']
IC_His_HelixPi_gminus_gminus = HIS_TEMPLATES['pi-helix']['g-/g-']
IC_His_HelixPi_gminus_trans = HIS_TEMPLATES['pi-helix']['g-/t']
IC_His_HelixPi_gminus_gplus = HIS_TEMPLATES['pi-helix']['g-/g+']
IC_His_HelixPi_trans_gminus = HIS_TEMPLATES['pi-helix']['t/g-']
IC_His_HelixPi_trans_trans = HIS_TEMPLATES['pi-helix']['t/t']
IC_His_HelixPi_trans_gplus = HIS_TEMPLATES['pi-helix']['t/g+']
IC_His_HelixPi_gplus_gminus = HIS_TEMPLATES['pi-helix']['g+/g-']
IC_His_HelixPi_gplus_trans = HIS_TEMPLATES['pi-helix']['g+/t']
IC_His_HelixPi_gplus_gplus = HIS_TEMPLATES['pi-helix']['g+/g+']
IC_His_HelixPPII = HIS_TEMPLATES['polyproline-II']['canonical']
IC_His_HelixPPII_gminus_gminus = HIS_TEMPLATES['polyproline-II']['g-/g-']
IC_His_HelixPPII_gminus_trans = HIS_TEMPLATES['polyproline-II']['g-/t']
IC_His_HelixPPII_gminus_gplus = HIS_TEMPLATES['polyproline-II']['g-/g+']
IC_His_HelixPPII_trans_gminus = HIS_TEMPLATES['polyproline-II']['t/g-']
IC_His_HelixPPII_trans_trans = HIS_TEMPLATES['polyproline-II']['t/t']
IC_His_HelixPPII_trans_gplus = HIS_TEMPLATES['polyproline-II']['t/g+']
IC_His_HelixPPII_gplus_gminus = HIS_TEMPLATES['polyproline-II']['g+/g-']
IC_His_HelixPPII_gplus_trans = HIS_TEMPLATES['polyproline-II']['g+/t']
IC_His_HelixPPII_gplus_gplus = HIS_TEMPLATES['polyproline-II']['g+/g+']
IC_His_Strand = HIS_TEMPLATES['beta-strand']['canonical']
IC_His_Strand_gminus_gminus = HIS_TEMPLATES['beta-strand']['g-/g-']
IC_His_Strand_gminus_trans = HIS_TEMPLATES['beta-strand']['g-/t']
IC_His_Strand_gminus_gplus = HIS_TEMPLATES['beta-strand']['g-/g+']
IC_His_Strand_trans_gminus = HIS_TEMPLATES['beta-strand']['t/g-']
IC_His_Strand_trans_trans = HIS_TEMPLATES['beta-strand']['t/t']
IC_His_Strand_trans_gplus = HIS_TEMPLATES['beta-strand']['t/g+']
IC_His_Strand_gplus_gminus = HIS_TEMPLATES['beta-strand']['g+/g-']
IC_His_Strand_gplus_trans = HIS_TEMPLATES['beta-strand']['g+/t']
IC_His_Strand_gplus_gplus = HIS_TEMPLATES['beta-strand']['g+/g+']
IC_His_StrandParallel = HIS_TEMPLATES['parallel-beta-strand']['canonical']
IC_His_StrandParallel_gminus_gminus = HIS_TEMPLATES['parallel-beta-strand']['g-/g-']
IC_His_StrandParallel_gminus_trans = HIS_TEMPLATES['parallel-beta-strand']['g-/t']
IC_His_StrandParallel_gminus_gplus = HIS_TEMPLATES['parallel-beta-strand']['g-/g+']
IC_His_StrandParallel_trans_gminus = HIS_TEMPLATES['parallel-beta-strand']['t/g-']
IC_His_StrandParallel_trans_trans = HIS_TEMPLATES['parallel-beta-strand']['t/t']
IC_His_StrandParallel_trans_gplus = HIS_TEMPLATES['parallel-beta-strand']['t/g+']
IC_His_StrandParallel_gplus_gminus = HIS_TEMPLATES['parallel-beta-strand']['g+/g-']
IC_His_StrandParallel_gplus_trans = HIS_TEMPLATES['parallel-beta-strand']['g+/t']
IC_His_StrandParallel_gplus_gplus = HIS_TEMPLATES['parallel-beta-strand']['g+/g+']
IC_His_StrandAntiParallel = HIS_TEMPLATES['antiparallel-beta-strand']['canonical']
IC_His_StrandAntiParallel_gminus_gminus = HIS_TEMPLATES['antiparallel-beta-strand']['g-/g-']
IC_His_StrandAntiParallel_gminus_trans = HIS_TEMPLATES['antiparallel-beta-strand']['g-/t']
IC_His_StrandAntiParallel_gminus_gplus = HIS_TEMPLATES['antiparallel-beta-strand']['g-/g+']
IC_His_StrandAntiParallel_trans_gminus = HIS_TEMPLATES['antiparallel-beta-strand']['t/g-']
IC_His_StrandAntiParallel_trans_trans = HIS_TEMPLATES['antiparallel-beta-strand']['t/t']
IC_His_StrandAntiParallel_trans_gplus = HIS_TEMPLATES['antiparallel-beta-strand']['t/g+']
IC_His_StrandAntiParallel_gplus_gminus = HIS_TEMPLATES['antiparallel-beta-strand']['g+/g-']
IC_His_StrandAntiParallel_gplus_trans = HIS_TEMPLATES['antiparallel-beta-strand']['g+/t']
IC_His_StrandAntiParallel_gplus_gplus = HIS_TEMPLATES['antiparallel-beta-strand']['g+/g+']
IC_His_Bridge = HIS_TEMPLATES['beta-bridge']['canonical']
IC_His_Bridge_gminus_gminus = HIS_TEMPLATES['beta-bridge']['g-/g-']
IC_His_Bridge_gminus_trans = HIS_TEMPLATES['beta-bridge']['g-/t']
IC_His_Bridge_gminus_gplus = HIS_TEMPLATES['beta-bridge']['g-/g+']
IC_His_Bridge_trans_gminus = HIS_TEMPLATES['beta-bridge']['t/g-']
IC_His_Bridge_trans_trans = HIS_TEMPLATES['beta-bridge']['t/t']
IC_His_Bridge_trans_gplus = HIS_TEMPLATES['beta-bridge']['t/g+']
IC_His_Bridge_gplus_gminus = HIS_TEMPLATES['beta-bridge']['g+/g-']
IC_His_Bridge_gplus_trans = HIS_TEMPLATES['beta-bridge']['g+/t']
IC_His_Bridge_gplus_gplus = HIS_TEMPLATES['beta-bridge']['g+/g+']
IC_His_Turn = HIS_TEMPLATES['turn']['canonical']
IC_His_Turn_gminus_gminus = HIS_TEMPLATES['turn']['g-/g-']
IC_His_Turn_gminus_trans = HIS_TEMPLATES['turn']['g-/t']
IC_His_Turn_gminus_gplus = HIS_TEMPLATES['turn']['g-/g+']
IC_His_Turn_trans_gminus = HIS_TEMPLATES['turn']['t/g-']
IC_His_Turn_trans_trans = HIS_TEMPLATES['turn']['t/t']
IC_His_Turn_trans_gplus = HIS_TEMPLATES['turn']['t/g+']
IC_His_Turn_gplus_gminus = HIS_TEMPLATES['turn']['g+/g-']
IC_His_Turn_gplus_trans = HIS_TEMPLATES['turn']['g+/t']
IC_His_Turn_gplus_gplus = HIS_TEMPLATES['turn']['g+/g+']
IC_His_Bend = HIS_TEMPLATES['bend']['canonical']
IC_His_Bend_gminus_gminus = HIS_TEMPLATES['bend']['g-/g-']
IC_His_Bend_gminus_trans = HIS_TEMPLATES['bend']['g-/t']
IC_His_Bend_gminus_gplus = HIS_TEMPLATES['bend']['g-/g+']
IC_His_Bend_trans_gminus = HIS_TEMPLATES['bend']['t/g-']
IC_His_Bend_trans_trans = HIS_TEMPLATES['bend']['t/t']
IC_His_Bend_trans_gplus = HIS_TEMPLATES['bend']['t/g+']
IC_His_Bend_gplus_gminus = HIS_TEMPLATES['bend']['g+/g-']
IC_His_Bend_gplus_trans = HIS_TEMPLATES['bend']['g+/t']
IC_His_Bend_gplus_gplus = HIS_TEMPLATES['bend']['g+/g+']
IC_His_Coil = HIS_TEMPLATES['coil']['canonical']
IC_His_Coil_gminus_gminus = HIS_TEMPLATES['coil']['g-/g-']
IC_His_Coil_gminus_trans = HIS_TEMPLATES['coil']['g-/t']
IC_His_Coil_gminus_gplus = HIS_TEMPLATES['coil']['g-/g+']
IC_His_Coil_trans_gminus = HIS_TEMPLATES['coil']['t/g-']
IC_His_Coil_trans_trans = HIS_TEMPLATES['coil']['t/t']
IC_His_Coil_trans_gplus = HIS_TEMPLATES['coil']['t/g+']
IC_His_Coil_gplus_gminus = HIS_TEMPLATES['coil']['g+/g-']
IC_His_Coil_gplus_trans = HIS_TEMPLATES['coil']['g+/t']
IC_His_Coil_gplus_gplus = HIS_TEMPLATES['coil']['g+/g+']
IC_His_CisPeptide = HIS_TEMPLATES['cis-peptide-bond']['canonical']
IC_His_CisPeptide_gminus_gminus = HIS_TEMPLATES['cis-peptide-bond']['g-/g-']
IC_His_CisPeptide_gminus_trans = HIS_TEMPLATES['cis-peptide-bond']['g-/t']
IC_His_CisPeptide_gminus_gplus = HIS_TEMPLATES['cis-peptide-bond']['g-/g+']
IC_His_CisPeptide_trans_gminus = HIS_TEMPLATES['cis-peptide-bond']['t/g-']
IC_His_CisPeptide_trans_trans = HIS_TEMPLATES['cis-peptide-bond']['t/t']
IC_His_CisPeptide_trans_gplus = HIS_TEMPLATES['cis-peptide-bond']['t/g+']
IC_His_CisPeptide_gplus_gminus = HIS_TEMPLATES['cis-peptide-bond']['g+/g-']
IC_His_CisPeptide_gplus_trans = HIS_TEMPLATES['cis-peptide-bond']['g+/t']
IC_His_CisPeptide_gplus_gplus = HIS_TEMPLATES['cis-peptide-bond']['g+/g+']

__all__ = [
    'IC_His_Bend',
    'IC_His_Bend_gminus_gminus',
    'IC_His_Bend_gminus_gplus',
    'IC_His_Bend_gminus_trans',
    'IC_His_Bend_gplus_gminus',
    'IC_His_Bend_gplus_gplus',
    'IC_His_Bend_gplus_trans',
    'IC_His_Bend_trans_gminus',
    'IC_His_Bend_trans_gplus',
    'IC_His_Bend_trans_trans',
    'IC_His_Bridge',
    'IC_His_Bridge_gminus_gminus',
    'IC_His_Bridge_gminus_gplus',
    'IC_His_Bridge_gminus_trans',
    'IC_His_Bridge_gplus_gminus',
    'IC_His_Bridge_gplus_gplus',
    'IC_His_Bridge_gplus_trans',
    'IC_His_Bridge_trans_gminus',
    'IC_His_Bridge_trans_gplus',
    'IC_His_Bridge_trans_trans',
    'IC_His_CisPeptide',
    'IC_His_CisPeptide_gminus_gminus',
    'IC_His_CisPeptide_gminus_gplus',
    'IC_His_CisPeptide_gminus_trans',
    'IC_His_CisPeptide_gplus_gminus',
    'IC_His_CisPeptide_gplus_gplus',
    'IC_His_CisPeptide_gplus_trans',
    'IC_His_CisPeptide_trans_gminus',
    'IC_His_CisPeptide_trans_gplus',
    'IC_His_CisPeptide_trans_trans',
    'IC_His_Coil',
    'IC_His_Coil_gminus_gminus',
    'IC_His_Coil_gminus_gplus',
    'IC_His_Coil_gminus_trans',
    'IC_His_Coil_gplus_gminus',
    'IC_His_Coil_gplus_gplus',
    'IC_His_Coil_gplus_trans',
    'IC_His_Coil_trans_gminus',
    'IC_His_Coil_trans_gplus',
    'IC_His_Coil_trans_trans',
    'IC_His_Helix310',
    'IC_His_Helix310_gminus_gminus',
    'IC_His_Helix310_gminus_gplus',
    'IC_His_Helix310_gminus_trans',
    'IC_His_Helix310_gplus_gminus',
    'IC_His_Helix310_gplus_gplus',
    'IC_His_Helix310_gplus_trans',
    'IC_His_Helix310_trans_gminus',
    'IC_His_Helix310_trans_gplus',
    'IC_His_Helix310_trans_trans',
    'IC_His_HelixAlpha',
    'IC_His_HelixAlpha_gminus_gminus',
    'IC_His_HelixAlpha_gminus_gplus',
    'IC_His_HelixAlpha_gminus_trans',
    'IC_His_HelixAlpha_gplus_gminus',
    'IC_His_HelixAlpha_gplus_gplus',
    'IC_His_HelixAlpha_gplus_trans',
    'IC_His_HelixAlpha_trans_gminus',
    'IC_His_HelixAlpha_trans_gplus',
    'IC_His_HelixAlpha_trans_trans',
    'IC_His_HelixPPII',
    'IC_His_HelixPPII_gminus_gminus',
    'IC_His_HelixPPII_gminus_gplus',
    'IC_His_HelixPPII_gminus_trans',
    'IC_His_HelixPPII_gplus_gminus',
    'IC_His_HelixPPII_gplus_gplus',
    'IC_His_HelixPPII_gplus_trans',
    'IC_His_HelixPPII_trans_gminus',
    'IC_His_HelixPPII_trans_gplus',
    'IC_His_HelixPPII_trans_trans',
    'IC_His_HelixPi',
    'IC_His_HelixPi_gminus_gminus',
    'IC_His_HelixPi_gminus_gplus',
    'IC_His_HelixPi_gminus_trans',
    'IC_His_HelixPi_gplus_gminus',
    'IC_His_HelixPi_gplus_gplus',
    'IC_His_HelixPi_gplus_trans',
    'IC_His_HelixPi_trans_gminus',
    'IC_His_HelixPi_trans_gplus',
    'IC_His_HelixPi_trans_trans',
    'IC_His_Strand',
    'IC_His_StrandAntiParallel',
    'IC_His_StrandAntiParallel_gminus_gminus',
    'IC_His_StrandAntiParallel_gminus_gplus',
    'IC_His_StrandAntiParallel_gminus_trans',
    'IC_His_StrandAntiParallel_gplus_gminus',
    'IC_His_StrandAntiParallel_gplus_gplus',
    'IC_His_StrandAntiParallel_gplus_trans',
    'IC_His_StrandAntiParallel_trans_gminus',
    'IC_His_StrandAntiParallel_trans_gplus',
    'IC_His_StrandAntiParallel_trans_trans',
    'IC_His_StrandParallel',
    'IC_His_StrandParallel_gminus_gminus',
    'IC_His_StrandParallel_gminus_gplus',
    'IC_His_StrandParallel_gminus_trans',
    'IC_His_StrandParallel_gplus_gminus',
    'IC_His_StrandParallel_gplus_gplus',
    'IC_His_StrandParallel_gplus_trans',
    'IC_His_StrandParallel_trans_gminus',
    'IC_His_StrandParallel_trans_gplus',
    'IC_His_StrandParallel_trans_trans',
    'IC_His_Strand_gminus_gminus',
    'IC_His_Strand_gminus_gplus',
    'IC_His_Strand_gminus_trans',
    'IC_His_Strand_gplus_gminus',
    'IC_His_Strand_gplus_gplus',
    'IC_His_Strand_gplus_trans',
    'IC_His_Strand_trans_gminus',
    'IC_His_Strand_trans_gplus',
    'IC_His_Strand_trans_trans',
    'IC_His_Turn',
    'IC_His_Turn_gminus_gminus',
    'IC_His_Turn_gminus_gplus',
    'IC_His_Turn_gminus_trans',
    'IC_His_Turn_gplus_gminus',
    'IC_His_Turn_gplus_gplus',
    'IC_His_Turn_gplus_trans',
    'IC_His_Turn_trans_gminus',
    'IC_His_Turn_trans_gplus',
    'IC_His_Turn_trans_trans',
]
