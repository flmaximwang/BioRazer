# -*- coding: utf-8 -*-
"""Ideal {secondary-structure x rotamer} {names[aa]} {var} ``InternalCoord`` templates.

Each template is the ideal heavy-atom topology for {aa} anchored at
{{N, CA, C}}, with the side chain grown to ``rotamer`` and the
``phi``/``psi``/``omega`` of ``ss`` carried on the instance.  See
:mod:`biorazer.database.internal_coord_template` and
:mod:`._builder` for the build rule and the per-conformer caveat.
"""

from biorazer.database.molecule.icoor.protein.template._builder import make_residue_templates

RESN = "LYS"
#: {ss: {rotamer: :class:`~biorazer.structure.objects.InternalCoord`}}
LYS_TEMPLATES = make_residue_templates(RESN)

# === GENERATED IC_* named templates - do not edit below ===
# Regenerate with: /opt/envs/BioRazer/bin/python scripts/generate_internal_coord_template_named.py

IC_Lys_HelixAlpha = LYS_TEMPLATES['alpha-helix']['canonical']
IC_Lys_HelixAlpha_gminus_gminus = LYS_TEMPLATES['alpha-helix']['g-/g-']
IC_Lys_HelixAlpha_gminus_trans = LYS_TEMPLATES['alpha-helix']['g-/t']
IC_Lys_HelixAlpha_gminus_gplus = LYS_TEMPLATES['alpha-helix']['g-/g+']
IC_Lys_HelixAlpha_trans_gminus = LYS_TEMPLATES['alpha-helix']['t/g-']
IC_Lys_HelixAlpha_trans_trans = LYS_TEMPLATES['alpha-helix']['t/t']
IC_Lys_HelixAlpha_trans_gplus = LYS_TEMPLATES['alpha-helix']['t/g+']
IC_Lys_HelixAlpha_gplus_gminus = LYS_TEMPLATES['alpha-helix']['g+/g-']
IC_Lys_HelixAlpha_gplus_trans = LYS_TEMPLATES['alpha-helix']['g+/t']
IC_Lys_HelixAlpha_gplus_gplus = LYS_TEMPLATES['alpha-helix']['g+/g+']
IC_Lys_Helix310 = LYS_TEMPLATES['3-10-helix']['canonical']
IC_Lys_Helix310_gminus_gminus = LYS_TEMPLATES['3-10-helix']['g-/g-']
IC_Lys_Helix310_gminus_trans = LYS_TEMPLATES['3-10-helix']['g-/t']
IC_Lys_Helix310_gminus_gplus = LYS_TEMPLATES['3-10-helix']['g-/g+']
IC_Lys_Helix310_trans_gminus = LYS_TEMPLATES['3-10-helix']['t/g-']
IC_Lys_Helix310_trans_trans = LYS_TEMPLATES['3-10-helix']['t/t']
IC_Lys_Helix310_trans_gplus = LYS_TEMPLATES['3-10-helix']['t/g+']
IC_Lys_Helix310_gplus_gminus = LYS_TEMPLATES['3-10-helix']['g+/g-']
IC_Lys_Helix310_gplus_trans = LYS_TEMPLATES['3-10-helix']['g+/t']
IC_Lys_Helix310_gplus_gplus = LYS_TEMPLATES['3-10-helix']['g+/g+']
IC_Lys_HelixPi = LYS_TEMPLATES['pi-helix']['canonical']
IC_Lys_HelixPi_gminus_gminus = LYS_TEMPLATES['pi-helix']['g-/g-']
IC_Lys_HelixPi_gminus_trans = LYS_TEMPLATES['pi-helix']['g-/t']
IC_Lys_HelixPi_gminus_gplus = LYS_TEMPLATES['pi-helix']['g-/g+']
IC_Lys_HelixPi_trans_gminus = LYS_TEMPLATES['pi-helix']['t/g-']
IC_Lys_HelixPi_trans_trans = LYS_TEMPLATES['pi-helix']['t/t']
IC_Lys_HelixPi_trans_gplus = LYS_TEMPLATES['pi-helix']['t/g+']
IC_Lys_HelixPi_gplus_gminus = LYS_TEMPLATES['pi-helix']['g+/g-']
IC_Lys_HelixPi_gplus_trans = LYS_TEMPLATES['pi-helix']['g+/t']
IC_Lys_HelixPi_gplus_gplus = LYS_TEMPLATES['pi-helix']['g+/g+']
IC_Lys_HelixPPII = LYS_TEMPLATES['polyproline-II']['canonical']
IC_Lys_HelixPPII_gminus_gminus = LYS_TEMPLATES['polyproline-II']['g-/g-']
IC_Lys_HelixPPII_gminus_trans = LYS_TEMPLATES['polyproline-II']['g-/t']
IC_Lys_HelixPPII_gminus_gplus = LYS_TEMPLATES['polyproline-II']['g-/g+']
IC_Lys_HelixPPII_trans_gminus = LYS_TEMPLATES['polyproline-II']['t/g-']
IC_Lys_HelixPPII_trans_trans = LYS_TEMPLATES['polyproline-II']['t/t']
IC_Lys_HelixPPII_trans_gplus = LYS_TEMPLATES['polyproline-II']['t/g+']
IC_Lys_HelixPPII_gplus_gminus = LYS_TEMPLATES['polyproline-II']['g+/g-']
IC_Lys_HelixPPII_gplus_trans = LYS_TEMPLATES['polyproline-II']['g+/t']
IC_Lys_HelixPPII_gplus_gplus = LYS_TEMPLATES['polyproline-II']['g+/g+']
IC_Lys_Strand = LYS_TEMPLATES['beta-strand']['canonical']
IC_Lys_Strand_gminus_gminus = LYS_TEMPLATES['beta-strand']['g-/g-']
IC_Lys_Strand_gminus_trans = LYS_TEMPLATES['beta-strand']['g-/t']
IC_Lys_Strand_gminus_gplus = LYS_TEMPLATES['beta-strand']['g-/g+']
IC_Lys_Strand_trans_gminus = LYS_TEMPLATES['beta-strand']['t/g-']
IC_Lys_Strand_trans_trans = LYS_TEMPLATES['beta-strand']['t/t']
IC_Lys_Strand_trans_gplus = LYS_TEMPLATES['beta-strand']['t/g+']
IC_Lys_Strand_gplus_gminus = LYS_TEMPLATES['beta-strand']['g+/g-']
IC_Lys_Strand_gplus_trans = LYS_TEMPLATES['beta-strand']['g+/t']
IC_Lys_Strand_gplus_gplus = LYS_TEMPLATES['beta-strand']['g+/g+']
IC_Lys_StrandParallel = LYS_TEMPLATES['parallel-beta-strand']['canonical']
IC_Lys_StrandParallel_gminus_gminus = LYS_TEMPLATES['parallel-beta-strand']['g-/g-']
IC_Lys_StrandParallel_gminus_trans = LYS_TEMPLATES['parallel-beta-strand']['g-/t']
IC_Lys_StrandParallel_gminus_gplus = LYS_TEMPLATES['parallel-beta-strand']['g-/g+']
IC_Lys_StrandParallel_trans_gminus = LYS_TEMPLATES['parallel-beta-strand']['t/g-']
IC_Lys_StrandParallel_trans_trans = LYS_TEMPLATES['parallel-beta-strand']['t/t']
IC_Lys_StrandParallel_trans_gplus = LYS_TEMPLATES['parallel-beta-strand']['t/g+']
IC_Lys_StrandParallel_gplus_gminus = LYS_TEMPLATES['parallel-beta-strand']['g+/g-']
IC_Lys_StrandParallel_gplus_trans = LYS_TEMPLATES['parallel-beta-strand']['g+/t']
IC_Lys_StrandParallel_gplus_gplus = LYS_TEMPLATES['parallel-beta-strand']['g+/g+']
IC_Lys_StrandAntiParallel = LYS_TEMPLATES['antiparallel-beta-strand']['canonical']
IC_Lys_StrandAntiParallel_gminus_gminus = LYS_TEMPLATES['antiparallel-beta-strand']['g-/g-']
IC_Lys_StrandAntiParallel_gminus_trans = LYS_TEMPLATES['antiparallel-beta-strand']['g-/t']
IC_Lys_StrandAntiParallel_gminus_gplus = LYS_TEMPLATES['antiparallel-beta-strand']['g-/g+']
IC_Lys_StrandAntiParallel_trans_gminus = LYS_TEMPLATES['antiparallel-beta-strand']['t/g-']
IC_Lys_StrandAntiParallel_trans_trans = LYS_TEMPLATES['antiparallel-beta-strand']['t/t']
IC_Lys_StrandAntiParallel_trans_gplus = LYS_TEMPLATES['antiparallel-beta-strand']['t/g+']
IC_Lys_StrandAntiParallel_gplus_gminus = LYS_TEMPLATES['antiparallel-beta-strand']['g+/g-']
IC_Lys_StrandAntiParallel_gplus_trans = LYS_TEMPLATES['antiparallel-beta-strand']['g+/t']
IC_Lys_StrandAntiParallel_gplus_gplus = LYS_TEMPLATES['antiparallel-beta-strand']['g+/g+']
IC_Lys_Bridge = LYS_TEMPLATES['beta-bridge']['canonical']
IC_Lys_Bridge_gminus_gminus = LYS_TEMPLATES['beta-bridge']['g-/g-']
IC_Lys_Bridge_gminus_trans = LYS_TEMPLATES['beta-bridge']['g-/t']
IC_Lys_Bridge_gminus_gplus = LYS_TEMPLATES['beta-bridge']['g-/g+']
IC_Lys_Bridge_trans_gminus = LYS_TEMPLATES['beta-bridge']['t/g-']
IC_Lys_Bridge_trans_trans = LYS_TEMPLATES['beta-bridge']['t/t']
IC_Lys_Bridge_trans_gplus = LYS_TEMPLATES['beta-bridge']['t/g+']
IC_Lys_Bridge_gplus_gminus = LYS_TEMPLATES['beta-bridge']['g+/g-']
IC_Lys_Bridge_gplus_trans = LYS_TEMPLATES['beta-bridge']['g+/t']
IC_Lys_Bridge_gplus_gplus = LYS_TEMPLATES['beta-bridge']['g+/g+']
IC_Lys_Turn = LYS_TEMPLATES['turn']['canonical']
IC_Lys_Turn_gminus_gminus = LYS_TEMPLATES['turn']['g-/g-']
IC_Lys_Turn_gminus_trans = LYS_TEMPLATES['turn']['g-/t']
IC_Lys_Turn_gminus_gplus = LYS_TEMPLATES['turn']['g-/g+']
IC_Lys_Turn_trans_gminus = LYS_TEMPLATES['turn']['t/g-']
IC_Lys_Turn_trans_trans = LYS_TEMPLATES['turn']['t/t']
IC_Lys_Turn_trans_gplus = LYS_TEMPLATES['turn']['t/g+']
IC_Lys_Turn_gplus_gminus = LYS_TEMPLATES['turn']['g+/g-']
IC_Lys_Turn_gplus_trans = LYS_TEMPLATES['turn']['g+/t']
IC_Lys_Turn_gplus_gplus = LYS_TEMPLATES['turn']['g+/g+']
IC_Lys_Bend = LYS_TEMPLATES['bend']['canonical']
IC_Lys_Bend_gminus_gminus = LYS_TEMPLATES['bend']['g-/g-']
IC_Lys_Bend_gminus_trans = LYS_TEMPLATES['bend']['g-/t']
IC_Lys_Bend_gminus_gplus = LYS_TEMPLATES['bend']['g-/g+']
IC_Lys_Bend_trans_gminus = LYS_TEMPLATES['bend']['t/g-']
IC_Lys_Bend_trans_trans = LYS_TEMPLATES['bend']['t/t']
IC_Lys_Bend_trans_gplus = LYS_TEMPLATES['bend']['t/g+']
IC_Lys_Bend_gplus_gminus = LYS_TEMPLATES['bend']['g+/g-']
IC_Lys_Bend_gplus_trans = LYS_TEMPLATES['bend']['g+/t']
IC_Lys_Bend_gplus_gplus = LYS_TEMPLATES['bend']['g+/g+']
IC_Lys_Coil = LYS_TEMPLATES['coil']['canonical']
IC_Lys_Coil_gminus_gminus = LYS_TEMPLATES['coil']['g-/g-']
IC_Lys_Coil_gminus_trans = LYS_TEMPLATES['coil']['g-/t']
IC_Lys_Coil_gminus_gplus = LYS_TEMPLATES['coil']['g-/g+']
IC_Lys_Coil_trans_gminus = LYS_TEMPLATES['coil']['t/g-']
IC_Lys_Coil_trans_trans = LYS_TEMPLATES['coil']['t/t']
IC_Lys_Coil_trans_gplus = LYS_TEMPLATES['coil']['t/g+']
IC_Lys_Coil_gplus_gminus = LYS_TEMPLATES['coil']['g+/g-']
IC_Lys_Coil_gplus_trans = LYS_TEMPLATES['coil']['g+/t']
IC_Lys_Coil_gplus_gplus = LYS_TEMPLATES['coil']['g+/g+']
IC_Lys_CisPeptide = LYS_TEMPLATES['cis-peptide-bond']['canonical']
IC_Lys_CisPeptide_gminus_gminus = LYS_TEMPLATES['cis-peptide-bond']['g-/g-']
IC_Lys_CisPeptide_gminus_trans = LYS_TEMPLATES['cis-peptide-bond']['g-/t']
IC_Lys_CisPeptide_gminus_gplus = LYS_TEMPLATES['cis-peptide-bond']['g-/g+']
IC_Lys_CisPeptide_trans_gminus = LYS_TEMPLATES['cis-peptide-bond']['t/g-']
IC_Lys_CisPeptide_trans_trans = LYS_TEMPLATES['cis-peptide-bond']['t/t']
IC_Lys_CisPeptide_trans_gplus = LYS_TEMPLATES['cis-peptide-bond']['t/g+']
IC_Lys_CisPeptide_gplus_gminus = LYS_TEMPLATES['cis-peptide-bond']['g+/g-']
IC_Lys_CisPeptide_gplus_trans = LYS_TEMPLATES['cis-peptide-bond']['g+/t']
IC_Lys_CisPeptide_gplus_gplus = LYS_TEMPLATES['cis-peptide-bond']['g+/g+']

__all__ = [
    'IC_Lys_Bend',
    'IC_Lys_Bend_gminus_gminus',
    'IC_Lys_Bend_gminus_gplus',
    'IC_Lys_Bend_gminus_trans',
    'IC_Lys_Bend_gplus_gminus',
    'IC_Lys_Bend_gplus_gplus',
    'IC_Lys_Bend_gplus_trans',
    'IC_Lys_Bend_trans_gminus',
    'IC_Lys_Bend_trans_gplus',
    'IC_Lys_Bend_trans_trans',
    'IC_Lys_Bridge',
    'IC_Lys_Bridge_gminus_gminus',
    'IC_Lys_Bridge_gminus_gplus',
    'IC_Lys_Bridge_gminus_trans',
    'IC_Lys_Bridge_gplus_gminus',
    'IC_Lys_Bridge_gplus_gplus',
    'IC_Lys_Bridge_gplus_trans',
    'IC_Lys_Bridge_trans_gminus',
    'IC_Lys_Bridge_trans_gplus',
    'IC_Lys_Bridge_trans_trans',
    'IC_Lys_CisPeptide',
    'IC_Lys_CisPeptide_gminus_gminus',
    'IC_Lys_CisPeptide_gminus_gplus',
    'IC_Lys_CisPeptide_gminus_trans',
    'IC_Lys_CisPeptide_gplus_gminus',
    'IC_Lys_CisPeptide_gplus_gplus',
    'IC_Lys_CisPeptide_gplus_trans',
    'IC_Lys_CisPeptide_trans_gminus',
    'IC_Lys_CisPeptide_trans_gplus',
    'IC_Lys_CisPeptide_trans_trans',
    'IC_Lys_Coil',
    'IC_Lys_Coil_gminus_gminus',
    'IC_Lys_Coil_gminus_gplus',
    'IC_Lys_Coil_gminus_trans',
    'IC_Lys_Coil_gplus_gminus',
    'IC_Lys_Coil_gplus_gplus',
    'IC_Lys_Coil_gplus_trans',
    'IC_Lys_Coil_trans_gminus',
    'IC_Lys_Coil_trans_gplus',
    'IC_Lys_Coil_trans_trans',
    'IC_Lys_Helix310',
    'IC_Lys_Helix310_gminus_gminus',
    'IC_Lys_Helix310_gminus_gplus',
    'IC_Lys_Helix310_gminus_trans',
    'IC_Lys_Helix310_gplus_gminus',
    'IC_Lys_Helix310_gplus_gplus',
    'IC_Lys_Helix310_gplus_trans',
    'IC_Lys_Helix310_trans_gminus',
    'IC_Lys_Helix310_trans_gplus',
    'IC_Lys_Helix310_trans_trans',
    'IC_Lys_HelixAlpha',
    'IC_Lys_HelixAlpha_gminus_gminus',
    'IC_Lys_HelixAlpha_gminus_gplus',
    'IC_Lys_HelixAlpha_gminus_trans',
    'IC_Lys_HelixAlpha_gplus_gminus',
    'IC_Lys_HelixAlpha_gplus_gplus',
    'IC_Lys_HelixAlpha_gplus_trans',
    'IC_Lys_HelixAlpha_trans_gminus',
    'IC_Lys_HelixAlpha_trans_gplus',
    'IC_Lys_HelixAlpha_trans_trans',
    'IC_Lys_HelixPPII',
    'IC_Lys_HelixPPII_gminus_gminus',
    'IC_Lys_HelixPPII_gminus_gplus',
    'IC_Lys_HelixPPII_gminus_trans',
    'IC_Lys_HelixPPII_gplus_gminus',
    'IC_Lys_HelixPPII_gplus_gplus',
    'IC_Lys_HelixPPII_gplus_trans',
    'IC_Lys_HelixPPII_trans_gminus',
    'IC_Lys_HelixPPII_trans_gplus',
    'IC_Lys_HelixPPII_trans_trans',
    'IC_Lys_HelixPi',
    'IC_Lys_HelixPi_gminus_gminus',
    'IC_Lys_HelixPi_gminus_gplus',
    'IC_Lys_HelixPi_gminus_trans',
    'IC_Lys_HelixPi_gplus_gminus',
    'IC_Lys_HelixPi_gplus_gplus',
    'IC_Lys_HelixPi_gplus_trans',
    'IC_Lys_HelixPi_trans_gminus',
    'IC_Lys_HelixPi_trans_gplus',
    'IC_Lys_HelixPi_trans_trans',
    'IC_Lys_Strand',
    'IC_Lys_StrandAntiParallel',
    'IC_Lys_StrandAntiParallel_gminus_gminus',
    'IC_Lys_StrandAntiParallel_gminus_gplus',
    'IC_Lys_StrandAntiParallel_gminus_trans',
    'IC_Lys_StrandAntiParallel_gplus_gminus',
    'IC_Lys_StrandAntiParallel_gplus_gplus',
    'IC_Lys_StrandAntiParallel_gplus_trans',
    'IC_Lys_StrandAntiParallel_trans_gminus',
    'IC_Lys_StrandAntiParallel_trans_gplus',
    'IC_Lys_StrandAntiParallel_trans_trans',
    'IC_Lys_StrandParallel',
    'IC_Lys_StrandParallel_gminus_gminus',
    'IC_Lys_StrandParallel_gminus_gplus',
    'IC_Lys_StrandParallel_gminus_trans',
    'IC_Lys_StrandParallel_gplus_gminus',
    'IC_Lys_StrandParallel_gplus_gplus',
    'IC_Lys_StrandParallel_gplus_trans',
    'IC_Lys_StrandParallel_trans_gminus',
    'IC_Lys_StrandParallel_trans_gplus',
    'IC_Lys_StrandParallel_trans_trans',
    'IC_Lys_Strand_gminus_gminus',
    'IC_Lys_Strand_gminus_gplus',
    'IC_Lys_Strand_gminus_trans',
    'IC_Lys_Strand_gplus_gminus',
    'IC_Lys_Strand_gplus_gplus',
    'IC_Lys_Strand_gplus_trans',
    'IC_Lys_Strand_trans_gminus',
    'IC_Lys_Strand_trans_gplus',
    'IC_Lys_Strand_trans_trans',
    'IC_Lys_Turn',
    'IC_Lys_Turn_gminus_gminus',
    'IC_Lys_Turn_gminus_gplus',
    'IC_Lys_Turn_gminus_trans',
    'IC_Lys_Turn_gplus_gminus',
    'IC_Lys_Turn_gplus_gplus',
    'IC_Lys_Turn_gplus_trans',
    'IC_Lys_Turn_trans_gminus',
    'IC_Lys_Turn_trans_gplus',
    'IC_Lys_Turn_trans_trans',
]
