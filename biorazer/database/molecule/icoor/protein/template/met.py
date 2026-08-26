# -*- coding: utf-8 -*-
"""Ideal {secondary-structure x rotamer} {names[aa]} {var} ``InternalCoord`` templates.

Each template is the ideal heavy-atom topology for {aa} anchored at
{{N, CA, C}}, with the side chain grown to ``rotamer`` and the
``phi``/``psi``/``omega`` of ``ss`` carried on the instance.  See
:mod:`biorazer.database.internal_coord_template` and
:mod:`._builder` for the build rule and the per-conformer caveat.
"""

from biorazer.database.molecule.icoor.protein.template._builder import make_residue_templates

RESN = "MET"
#: {ss: {rotamer: :class:`~biorazer.structure.objects.InternalCoord`}}
MET_TEMPLATES = make_residue_templates(RESN)

# === GENERATED IC_* named templates - do not edit below ===
# Regenerate with: /opt/envs/BioRazer/bin/python scripts/generate_internal_coord_template_named.py

IC_Met_HelixAlpha = MET_TEMPLATES['alpha-helix']['canonical']
IC_Met_HelixAlpha_gminus_gminus = MET_TEMPLATES['alpha-helix']['g-/g-']
IC_Met_HelixAlpha_gminus_trans = MET_TEMPLATES['alpha-helix']['g-/t']
IC_Met_HelixAlpha_gminus_gplus = MET_TEMPLATES['alpha-helix']['g-/g+']
IC_Met_HelixAlpha_trans_gminus = MET_TEMPLATES['alpha-helix']['t/g-']
IC_Met_HelixAlpha_trans_trans = MET_TEMPLATES['alpha-helix']['t/t']
IC_Met_HelixAlpha_trans_gplus = MET_TEMPLATES['alpha-helix']['t/g+']
IC_Met_HelixAlpha_gplus_gminus = MET_TEMPLATES['alpha-helix']['g+/g-']
IC_Met_HelixAlpha_gplus_trans = MET_TEMPLATES['alpha-helix']['g+/t']
IC_Met_HelixAlpha_gplus_gplus = MET_TEMPLATES['alpha-helix']['g+/g+']
IC_Met_Helix310 = MET_TEMPLATES['3-10-helix']['canonical']
IC_Met_Helix310_gminus_gminus = MET_TEMPLATES['3-10-helix']['g-/g-']
IC_Met_Helix310_gminus_trans = MET_TEMPLATES['3-10-helix']['g-/t']
IC_Met_Helix310_gminus_gplus = MET_TEMPLATES['3-10-helix']['g-/g+']
IC_Met_Helix310_trans_gminus = MET_TEMPLATES['3-10-helix']['t/g-']
IC_Met_Helix310_trans_trans = MET_TEMPLATES['3-10-helix']['t/t']
IC_Met_Helix310_trans_gplus = MET_TEMPLATES['3-10-helix']['t/g+']
IC_Met_Helix310_gplus_gminus = MET_TEMPLATES['3-10-helix']['g+/g-']
IC_Met_Helix310_gplus_trans = MET_TEMPLATES['3-10-helix']['g+/t']
IC_Met_Helix310_gplus_gplus = MET_TEMPLATES['3-10-helix']['g+/g+']
IC_Met_HelixPi = MET_TEMPLATES['pi-helix']['canonical']
IC_Met_HelixPi_gminus_gminus = MET_TEMPLATES['pi-helix']['g-/g-']
IC_Met_HelixPi_gminus_trans = MET_TEMPLATES['pi-helix']['g-/t']
IC_Met_HelixPi_gminus_gplus = MET_TEMPLATES['pi-helix']['g-/g+']
IC_Met_HelixPi_trans_gminus = MET_TEMPLATES['pi-helix']['t/g-']
IC_Met_HelixPi_trans_trans = MET_TEMPLATES['pi-helix']['t/t']
IC_Met_HelixPi_trans_gplus = MET_TEMPLATES['pi-helix']['t/g+']
IC_Met_HelixPi_gplus_gminus = MET_TEMPLATES['pi-helix']['g+/g-']
IC_Met_HelixPi_gplus_trans = MET_TEMPLATES['pi-helix']['g+/t']
IC_Met_HelixPi_gplus_gplus = MET_TEMPLATES['pi-helix']['g+/g+']
IC_Met_HelixPPII = MET_TEMPLATES['polyproline-II']['canonical']
IC_Met_HelixPPII_gminus_gminus = MET_TEMPLATES['polyproline-II']['g-/g-']
IC_Met_HelixPPII_gminus_trans = MET_TEMPLATES['polyproline-II']['g-/t']
IC_Met_HelixPPII_gminus_gplus = MET_TEMPLATES['polyproline-II']['g-/g+']
IC_Met_HelixPPII_trans_gminus = MET_TEMPLATES['polyproline-II']['t/g-']
IC_Met_HelixPPII_trans_trans = MET_TEMPLATES['polyproline-II']['t/t']
IC_Met_HelixPPII_trans_gplus = MET_TEMPLATES['polyproline-II']['t/g+']
IC_Met_HelixPPII_gplus_gminus = MET_TEMPLATES['polyproline-II']['g+/g-']
IC_Met_HelixPPII_gplus_trans = MET_TEMPLATES['polyproline-II']['g+/t']
IC_Met_HelixPPII_gplus_gplus = MET_TEMPLATES['polyproline-II']['g+/g+']
IC_Met_Strand = MET_TEMPLATES['beta-strand']['canonical']
IC_Met_Strand_gminus_gminus = MET_TEMPLATES['beta-strand']['g-/g-']
IC_Met_Strand_gminus_trans = MET_TEMPLATES['beta-strand']['g-/t']
IC_Met_Strand_gminus_gplus = MET_TEMPLATES['beta-strand']['g-/g+']
IC_Met_Strand_trans_gminus = MET_TEMPLATES['beta-strand']['t/g-']
IC_Met_Strand_trans_trans = MET_TEMPLATES['beta-strand']['t/t']
IC_Met_Strand_trans_gplus = MET_TEMPLATES['beta-strand']['t/g+']
IC_Met_Strand_gplus_gminus = MET_TEMPLATES['beta-strand']['g+/g-']
IC_Met_Strand_gplus_trans = MET_TEMPLATES['beta-strand']['g+/t']
IC_Met_Strand_gplus_gplus = MET_TEMPLATES['beta-strand']['g+/g+']
IC_Met_StrandParallel = MET_TEMPLATES['parallel-beta-strand']['canonical']
IC_Met_StrandParallel_gminus_gminus = MET_TEMPLATES['parallel-beta-strand']['g-/g-']
IC_Met_StrandParallel_gminus_trans = MET_TEMPLATES['parallel-beta-strand']['g-/t']
IC_Met_StrandParallel_gminus_gplus = MET_TEMPLATES['parallel-beta-strand']['g-/g+']
IC_Met_StrandParallel_trans_gminus = MET_TEMPLATES['parallel-beta-strand']['t/g-']
IC_Met_StrandParallel_trans_trans = MET_TEMPLATES['parallel-beta-strand']['t/t']
IC_Met_StrandParallel_trans_gplus = MET_TEMPLATES['parallel-beta-strand']['t/g+']
IC_Met_StrandParallel_gplus_gminus = MET_TEMPLATES['parallel-beta-strand']['g+/g-']
IC_Met_StrandParallel_gplus_trans = MET_TEMPLATES['parallel-beta-strand']['g+/t']
IC_Met_StrandParallel_gplus_gplus = MET_TEMPLATES['parallel-beta-strand']['g+/g+']
IC_Met_StrandAntiParallel = MET_TEMPLATES['antiparallel-beta-strand']['canonical']
IC_Met_StrandAntiParallel_gminus_gminus = MET_TEMPLATES['antiparallel-beta-strand']['g-/g-']
IC_Met_StrandAntiParallel_gminus_trans = MET_TEMPLATES['antiparallel-beta-strand']['g-/t']
IC_Met_StrandAntiParallel_gminus_gplus = MET_TEMPLATES['antiparallel-beta-strand']['g-/g+']
IC_Met_StrandAntiParallel_trans_gminus = MET_TEMPLATES['antiparallel-beta-strand']['t/g-']
IC_Met_StrandAntiParallel_trans_trans = MET_TEMPLATES['antiparallel-beta-strand']['t/t']
IC_Met_StrandAntiParallel_trans_gplus = MET_TEMPLATES['antiparallel-beta-strand']['t/g+']
IC_Met_StrandAntiParallel_gplus_gminus = MET_TEMPLATES['antiparallel-beta-strand']['g+/g-']
IC_Met_StrandAntiParallel_gplus_trans = MET_TEMPLATES['antiparallel-beta-strand']['g+/t']
IC_Met_StrandAntiParallel_gplus_gplus = MET_TEMPLATES['antiparallel-beta-strand']['g+/g+']
IC_Met_Bridge = MET_TEMPLATES['beta-bridge']['canonical']
IC_Met_Bridge_gminus_gminus = MET_TEMPLATES['beta-bridge']['g-/g-']
IC_Met_Bridge_gminus_trans = MET_TEMPLATES['beta-bridge']['g-/t']
IC_Met_Bridge_gminus_gplus = MET_TEMPLATES['beta-bridge']['g-/g+']
IC_Met_Bridge_trans_gminus = MET_TEMPLATES['beta-bridge']['t/g-']
IC_Met_Bridge_trans_trans = MET_TEMPLATES['beta-bridge']['t/t']
IC_Met_Bridge_trans_gplus = MET_TEMPLATES['beta-bridge']['t/g+']
IC_Met_Bridge_gplus_gminus = MET_TEMPLATES['beta-bridge']['g+/g-']
IC_Met_Bridge_gplus_trans = MET_TEMPLATES['beta-bridge']['g+/t']
IC_Met_Bridge_gplus_gplus = MET_TEMPLATES['beta-bridge']['g+/g+']
IC_Met_Turn = MET_TEMPLATES['turn']['canonical']
IC_Met_Turn_gminus_gminus = MET_TEMPLATES['turn']['g-/g-']
IC_Met_Turn_gminus_trans = MET_TEMPLATES['turn']['g-/t']
IC_Met_Turn_gminus_gplus = MET_TEMPLATES['turn']['g-/g+']
IC_Met_Turn_trans_gminus = MET_TEMPLATES['turn']['t/g-']
IC_Met_Turn_trans_trans = MET_TEMPLATES['turn']['t/t']
IC_Met_Turn_trans_gplus = MET_TEMPLATES['turn']['t/g+']
IC_Met_Turn_gplus_gminus = MET_TEMPLATES['turn']['g+/g-']
IC_Met_Turn_gplus_trans = MET_TEMPLATES['turn']['g+/t']
IC_Met_Turn_gplus_gplus = MET_TEMPLATES['turn']['g+/g+']
IC_Met_Bend = MET_TEMPLATES['bend']['canonical']
IC_Met_Bend_gminus_gminus = MET_TEMPLATES['bend']['g-/g-']
IC_Met_Bend_gminus_trans = MET_TEMPLATES['bend']['g-/t']
IC_Met_Bend_gminus_gplus = MET_TEMPLATES['bend']['g-/g+']
IC_Met_Bend_trans_gminus = MET_TEMPLATES['bend']['t/g-']
IC_Met_Bend_trans_trans = MET_TEMPLATES['bend']['t/t']
IC_Met_Bend_trans_gplus = MET_TEMPLATES['bend']['t/g+']
IC_Met_Bend_gplus_gminus = MET_TEMPLATES['bend']['g+/g-']
IC_Met_Bend_gplus_trans = MET_TEMPLATES['bend']['g+/t']
IC_Met_Bend_gplus_gplus = MET_TEMPLATES['bend']['g+/g+']
IC_Met_Coil = MET_TEMPLATES['coil']['canonical']
IC_Met_Coil_gminus_gminus = MET_TEMPLATES['coil']['g-/g-']
IC_Met_Coil_gminus_trans = MET_TEMPLATES['coil']['g-/t']
IC_Met_Coil_gminus_gplus = MET_TEMPLATES['coil']['g-/g+']
IC_Met_Coil_trans_gminus = MET_TEMPLATES['coil']['t/g-']
IC_Met_Coil_trans_trans = MET_TEMPLATES['coil']['t/t']
IC_Met_Coil_trans_gplus = MET_TEMPLATES['coil']['t/g+']
IC_Met_Coil_gplus_gminus = MET_TEMPLATES['coil']['g+/g-']
IC_Met_Coil_gplus_trans = MET_TEMPLATES['coil']['g+/t']
IC_Met_Coil_gplus_gplus = MET_TEMPLATES['coil']['g+/g+']
IC_Met_CisPeptide = MET_TEMPLATES['cis-peptide-bond']['canonical']
IC_Met_CisPeptide_gminus_gminus = MET_TEMPLATES['cis-peptide-bond']['g-/g-']
IC_Met_CisPeptide_gminus_trans = MET_TEMPLATES['cis-peptide-bond']['g-/t']
IC_Met_CisPeptide_gminus_gplus = MET_TEMPLATES['cis-peptide-bond']['g-/g+']
IC_Met_CisPeptide_trans_gminus = MET_TEMPLATES['cis-peptide-bond']['t/g-']
IC_Met_CisPeptide_trans_trans = MET_TEMPLATES['cis-peptide-bond']['t/t']
IC_Met_CisPeptide_trans_gplus = MET_TEMPLATES['cis-peptide-bond']['t/g+']
IC_Met_CisPeptide_gplus_gminus = MET_TEMPLATES['cis-peptide-bond']['g+/g-']
IC_Met_CisPeptide_gplus_trans = MET_TEMPLATES['cis-peptide-bond']['g+/t']
IC_Met_CisPeptide_gplus_gplus = MET_TEMPLATES['cis-peptide-bond']['g+/g+']

__all__ = [
    'IC_Met_Bend',
    'IC_Met_Bend_gminus_gminus',
    'IC_Met_Bend_gminus_gplus',
    'IC_Met_Bend_gminus_trans',
    'IC_Met_Bend_gplus_gminus',
    'IC_Met_Bend_gplus_gplus',
    'IC_Met_Bend_gplus_trans',
    'IC_Met_Bend_trans_gminus',
    'IC_Met_Bend_trans_gplus',
    'IC_Met_Bend_trans_trans',
    'IC_Met_Bridge',
    'IC_Met_Bridge_gminus_gminus',
    'IC_Met_Bridge_gminus_gplus',
    'IC_Met_Bridge_gminus_trans',
    'IC_Met_Bridge_gplus_gminus',
    'IC_Met_Bridge_gplus_gplus',
    'IC_Met_Bridge_gplus_trans',
    'IC_Met_Bridge_trans_gminus',
    'IC_Met_Bridge_trans_gplus',
    'IC_Met_Bridge_trans_trans',
    'IC_Met_CisPeptide',
    'IC_Met_CisPeptide_gminus_gminus',
    'IC_Met_CisPeptide_gminus_gplus',
    'IC_Met_CisPeptide_gminus_trans',
    'IC_Met_CisPeptide_gplus_gminus',
    'IC_Met_CisPeptide_gplus_gplus',
    'IC_Met_CisPeptide_gplus_trans',
    'IC_Met_CisPeptide_trans_gminus',
    'IC_Met_CisPeptide_trans_gplus',
    'IC_Met_CisPeptide_trans_trans',
    'IC_Met_Coil',
    'IC_Met_Coil_gminus_gminus',
    'IC_Met_Coil_gminus_gplus',
    'IC_Met_Coil_gminus_trans',
    'IC_Met_Coil_gplus_gminus',
    'IC_Met_Coil_gplus_gplus',
    'IC_Met_Coil_gplus_trans',
    'IC_Met_Coil_trans_gminus',
    'IC_Met_Coil_trans_gplus',
    'IC_Met_Coil_trans_trans',
    'IC_Met_Helix310',
    'IC_Met_Helix310_gminus_gminus',
    'IC_Met_Helix310_gminus_gplus',
    'IC_Met_Helix310_gminus_trans',
    'IC_Met_Helix310_gplus_gminus',
    'IC_Met_Helix310_gplus_gplus',
    'IC_Met_Helix310_gplus_trans',
    'IC_Met_Helix310_trans_gminus',
    'IC_Met_Helix310_trans_gplus',
    'IC_Met_Helix310_trans_trans',
    'IC_Met_HelixAlpha',
    'IC_Met_HelixAlpha_gminus_gminus',
    'IC_Met_HelixAlpha_gminus_gplus',
    'IC_Met_HelixAlpha_gminus_trans',
    'IC_Met_HelixAlpha_gplus_gminus',
    'IC_Met_HelixAlpha_gplus_gplus',
    'IC_Met_HelixAlpha_gplus_trans',
    'IC_Met_HelixAlpha_trans_gminus',
    'IC_Met_HelixAlpha_trans_gplus',
    'IC_Met_HelixAlpha_trans_trans',
    'IC_Met_HelixPPII',
    'IC_Met_HelixPPII_gminus_gminus',
    'IC_Met_HelixPPII_gminus_gplus',
    'IC_Met_HelixPPII_gminus_trans',
    'IC_Met_HelixPPII_gplus_gminus',
    'IC_Met_HelixPPII_gplus_gplus',
    'IC_Met_HelixPPII_gplus_trans',
    'IC_Met_HelixPPII_trans_gminus',
    'IC_Met_HelixPPII_trans_gplus',
    'IC_Met_HelixPPII_trans_trans',
    'IC_Met_HelixPi',
    'IC_Met_HelixPi_gminus_gminus',
    'IC_Met_HelixPi_gminus_gplus',
    'IC_Met_HelixPi_gminus_trans',
    'IC_Met_HelixPi_gplus_gminus',
    'IC_Met_HelixPi_gplus_gplus',
    'IC_Met_HelixPi_gplus_trans',
    'IC_Met_HelixPi_trans_gminus',
    'IC_Met_HelixPi_trans_gplus',
    'IC_Met_HelixPi_trans_trans',
    'IC_Met_Strand',
    'IC_Met_StrandAntiParallel',
    'IC_Met_StrandAntiParallel_gminus_gminus',
    'IC_Met_StrandAntiParallel_gminus_gplus',
    'IC_Met_StrandAntiParallel_gminus_trans',
    'IC_Met_StrandAntiParallel_gplus_gminus',
    'IC_Met_StrandAntiParallel_gplus_gplus',
    'IC_Met_StrandAntiParallel_gplus_trans',
    'IC_Met_StrandAntiParallel_trans_gminus',
    'IC_Met_StrandAntiParallel_trans_gplus',
    'IC_Met_StrandAntiParallel_trans_trans',
    'IC_Met_StrandParallel',
    'IC_Met_StrandParallel_gminus_gminus',
    'IC_Met_StrandParallel_gminus_gplus',
    'IC_Met_StrandParallel_gminus_trans',
    'IC_Met_StrandParallel_gplus_gminus',
    'IC_Met_StrandParallel_gplus_gplus',
    'IC_Met_StrandParallel_gplus_trans',
    'IC_Met_StrandParallel_trans_gminus',
    'IC_Met_StrandParallel_trans_gplus',
    'IC_Met_StrandParallel_trans_trans',
    'IC_Met_Strand_gminus_gminus',
    'IC_Met_Strand_gminus_gplus',
    'IC_Met_Strand_gminus_trans',
    'IC_Met_Strand_gplus_gminus',
    'IC_Met_Strand_gplus_gplus',
    'IC_Met_Strand_gplus_trans',
    'IC_Met_Strand_trans_gminus',
    'IC_Met_Strand_trans_gplus',
    'IC_Met_Strand_trans_trans',
    'IC_Met_Turn',
    'IC_Met_Turn_gminus_gminus',
    'IC_Met_Turn_gminus_gplus',
    'IC_Met_Turn_gminus_trans',
    'IC_Met_Turn_gplus_gminus',
    'IC_Met_Turn_gplus_gplus',
    'IC_Met_Turn_gplus_trans',
    'IC_Met_Turn_trans_gminus',
    'IC_Met_Turn_trans_gplus',
    'IC_Met_Turn_trans_trans',
]
