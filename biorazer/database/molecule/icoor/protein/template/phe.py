# -*- coding: utf-8 -*-
"""Ideal {secondary-structure x rotamer} {names[aa]} {var} ``InternalCoord`` templates.

Each template is the ideal heavy-atom topology for {aa} anchored at
{{N, CA, C}}, with the side chain grown to ``rotamer`` and the
``phi``/``psi``/``omega`` of ``ss`` carried on the instance.  See
:mod:`biorazer.database.internal_coord_template` and
:mod:`._builder` for the build rule and the per-conformer caveat.
"""

from biorazer.database.molecule.icoor.protein.template._builder import make_residue_templates

RESN = "PHE"
#: {ss: {rotamer: :class:`~biorazer.structure.objects.InternalCoord`}}
PHE_TEMPLATES = make_residue_templates(RESN)

# === GENERATED IC_* named templates - do not edit below ===
# Regenerate with: /opt/envs/BioRazer/bin/python scripts/generate_internal_coord_template_named.py

IC_Phe_HelixAlpha = PHE_TEMPLATES['alpha-helix']['canonical']
IC_Phe_HelixAlpha_gminus_gminus = PHE_TEMPLATES['alpha-helix']['g-/g-']
IC_Phe_HelixAlpha_gminus_trans = PHE_TEMPLATES['alpha-helix']['g-/t']
IC_Phe_HelixAlpha_gminus_gplus = PHE_TEMPLATES['alpha-helix']['g-/g+']
IC_Phe_HelixAlpha_trans_gminus = PHE_TEMPLATES['alpha-helix']['t/g-']
IC_Phe_HelixAlpha_trans_trans = PHE_TEMPLATES['alpha-helix']['t/t']
IC_Phe_HelixAlpha_trans_gplus = PHE_TEMPLATES['alpha-helix']['t/g+']
IC_Phe_HelixAlpha_gplus_gminus = PHE_TEMPLATES['alpha-helix']['g+/g-']
IC_Phe_HelixAlpha_gplus_trans = PHE_TEMPLATES['alpha-helix']['g+/t']
IC_Phe_HelixAlpha_gplus_gplus = PHE_TEMPLATES['alpha-helix']['g+/g+']
IC_Phe_Helix310 = PHE_TEMPLATES['3-10-helix']['canonical']
IC_Phe_Helix310_gminus_gminus = PHE_TEMPLATES['3-10-helix']['g-/g-']
IC_Phe_Helix310_gminus_trans = PHE_TEMPLATES['3-10-helix']['g-/t']
IC_Phe_Helix310_gminus_gplus = PHE_TEMPLATES['3-10-helix']['g-/g+']
IC_Phe_Helix310_trans_gminus = PHE_TEMPLATES['3-10-helix']['t/g-']
IC_Phe_Helix310_trans_trans = PHE_TEMPLATES['3-10-helix']['t/t']
IC_Phe_Helix310_trans_gplus = PHE_TEMPLATES['3-10-helix']['t/g+']
IC_Phe_Helix310_gplus_gminus = PHE_TEMPLATES['3-10-helix']['g+/g-']
IC_Phe_Helix310_gplus_trans = PHE_TEMPLATES['3-10-helix']['g+/t']
IC_Phe_Helix310_gplus_gplus = PHE_TEMPLATES['3-10-helix']['g+/g+']
IC_Phe_HelixPi = PHE_TEMPLATES['pi-helix']['canonical']
IC_Phe_HelixPi_gminus_gminus = PHE_TEMPLATES['pi-helix']['g-/g-']
IC_Phe_HelixPi_gminus_trans = PHE_TEMPLATES['pi-helix']['g-/t']
IC_Phe_HelixPi_gminus_gplus = PHE_TEMPLATES['pi-helix']['g-/g+']
IC_Phe_HelixPi_trans_gminus = PHE_TEMPLATES['pi-helix']['t/g-']
IC_Phe_HelixPi_trans_trans = PHE_TEMPLATES['pi-helix']['t/t']
IC_Phe_HelixPi_trans_gplus = PHE_TEMPLATES['pi-helix']['t/g+']
IC_Phe_HelixPi_gplus_gminus = PHE_TEMPLATES['pi-helix']['g+/g-']
IC_Phe_HelixPi_gplus_trans = PHE_TEMPLATES['pi-helix']['g+/t']
IC_Phe_HelixPi_gplus_gplus = PHE_TEMPLATES['pi-helix']['g+/g+']
IC_Phe_HelixPPII = PHE_TEMPLATES['polyproline-II']['canonical']
IC_Phe_HelixPPII_gminus_gminus = PHE_TEMPLATES['polyproline-II']['g-/g-']
IC_Phe_HelixPPII_gminus_trans = PHE_TEMPLATES['polyproline-II']['g-/t']
IC_Phe_HelixPPII_gminus_gplus = PHE_TEMPLATES['polyproline-II']['g-/g+']
IC_Phe_HelixPPII_trans_gminus = PHE_TEMPLATES['polyproline-II']['t/g-']
IC_Phe_HelixPPII_trans_trans = PHE_TEMPLATES['polyproline-II']['t/t']
IC_Phe_HelixPPII_trans_gplus = PHE_TEMPLATES['polyproline-II']['t/g+']
IC_Phe_HelixPPII_gplus_gminus = PHE_TEMPLATES['polyproline-II']['g+/g-']
IC_Phe_HelixPPII_gplus_trans = PHE_TEMPLATES['polyproline-II']['g+/t']
IC_Phe_HelixPPII_gplus_gplus = PHE_TEMPLATES['polyproline-II']['g+/g+']
IC_Phe_Strand = PHE_TEMPLATES['beta-strand']['canonical']
IC_Phe_Strand_gminus_gminus = PHE_TEMPLATES['beta-strand']['g-/g-']
IC_Phe_Strand_gminus_trans = PHE_TEMPLATES['beta-strand']['g-/t']
IC_Phe_Strand_gminus_gplus = PHE_TEMPLATES['beta-strand']['g-/g+']
IC_Phe_Strand_trans_gminus = PHE_TEMPLATES['beta-strand']['t/g-']
IC_Phe_Strand_trans_trans = PHE_TEMPLATES['beta-strand']['t/t']
IC_Phe_Strand_trans_gplus = PHE_TEMPLATES['beta-strand']['t/g+']
IC_Phe_Strand_gplus_gminus = PHE_TEMPLATES['beta-strand']['g+/g-']
IC_Phe_Strand_gplus_trans = PHE_TEMPLATES['beta-strand']['g+/t']
IC_Phe_Strand_gplus_gplus = PHE_TEMPLATES['beta-strand']['g+/g+']
IC_Phe_StrandParallel = PHE_TEMPLATES['parallel-beta-strand']['canonical']
IC_Phe_StrandParallel_gminus_gminus = PHE_TEMPLATES['parallel-beta-strand']['g-/g-']
IC_Phe_StrandParallel_gminus_trans = PHE_TEMPLATES['parallel-beta-strand']['g-/t']
IC_Phe_StrandParallel_gminus_gplus = PHE_TEMPLATES['parallel-beta-strand']['g-/g+']
IC_Phe_StrandParallel_trans_gminus = PHE_TEMPLATES['parallel-beta-strand']['t/g-']
IC_Phe_StrandParallel_trans_trans = PHE_TEMPLATES['parallel-beta-strand']['t/t']
IC_Phe_StrandParallel_trans_gplus = PHE_TEMPLATES['parallel-beta-strand']['t/g+']
IC_Phe_StrandParallel_gplus_gminus = PHE_TEMPLATES['parallel-beta-strand']['g+/g-']
IC_Phe_StrandParallel_gplus_trans = PHE_TEMPLATES['parallel-beta-strand']['g+/t']
IC_Phe_StrandParallel_gplus_gplus = PHE_TEMPLATES['parallel-beta-strand']['g+/g+']
IC_Phe_StrandAntiParallel = PHE_TEMPLATES['antiparallel-beta-strand']['canonical']
IC_Phe_StrandAntiParallel_gminus_gminus = PHE_TEMPLATES['antiparallel-beta-strand']['g-/g-']
IC_Phe_StrandAntiParallel_gminus_trans = PHE_TEMPLATES['antiparallel-beta-strand']['g-/t']
IC_Phe_StrandAntiParallel_gminus_gplus = PHE_TEMPLATES['antiparallel-beta-strand']['g-/g+']
IC_Phe_StrandAntiParallel_trans_gminus = PHE_TEMPLATES['antiparallel-beta-strand']['t/g-']
IC_Phe_StrandAntiParallel_trans_trans = PHE_TEMPLATES['antiparallel-beta-strand']['t/t']
IC_Phe_StrandAntiParallel_trans_gplus = PHE_TEMPLATES['antiparallel-beta-strand']['t/g+']
IC_Phe_StrandAntiParallel_gplus_gminus = PHE_TEMPLATES['antiparallel-beta-strand']['g+/g-']
IC_Phe_StrandAntiParallel_gplus_trans = PHE_TEMPLATES['antiparallel-beta-strand']['g+/t']
IC_Phe_StrandAntiParallel_gplus_gplus = PHE_TEMPLATES['antiparallel-beta-strand']['g+/g+']
IC_Phe_Bridge = PHE_TEMPLATES['beta-bridge']['canonical']
IC_Phe_Bridge_gminus_gminus = PHE_TEMPLATES['beta-bridge']['g-/g-']
IC_Phe_Bridge_gminus_trans = PHE_TEMPLATES['beta-bridge']['g-/t']
IC_Phe_Bridge_gminus_gplus = PHE_TEMPLATES['beta-bridge']['g-/g+']
IC_Phe_Bridge_trans_gminus = PHE_TEMPLATES['beta-bridge']['t/g-']
IC_Phe_Bridge_trans_trans = PHE_TEMPLATES['beta-bridge']['t/t']
IC_Phe_Bridge_trans_gplus = PHE_TEMPLATES['beta-bridge']['t/g+']
IC_Phe_Bridge_gplus_gminus = PHE_TEMPLATES['beta-bridge']['g+/g-']
IC_Phe_Bridge_gplus_trans = PHE_TEMPLATES['beta-bridge']['g+/t']
IC_Phe_Bridge_gplus_gplus = PHE_TEMPLATES['beta-bridge']['g+/g+']
IC_Phe_Turn = PHE_TEMPLATES['turn']['canonical']
IC_Phe_Turn_gminus_gminus = PHE_TEMPLATES['turn']['g-/g-']
IC_Phe_Turn_gminus_trans = PHE_TEMPLATES['turn']['g-/t']
IC_Phe_Turn_gminus_gplus = PHE_TEMPLATES['turn']['g-/g+']
IC_Phe_Turn_trans_gminus = PHE_TEMPLATES['turn']['t/g-']
IC_Phe_Turn_trans_trans = PHE_TEMPLATES['turn']['t/t']
IC_Phe_Turn_trans_gplus = PHE_TEMPLATES['turn']['t/g+']
IC_Phe_Turn_gplus_gminus = PHE_TEMPLATES['turn']['g+/g-']
IC_Phe_Turn_gplus_trans = PHE_TEMPLATES['turn']['g+/t']
IC_Phe_Turn_gplus_gplus = PHE_TEMPLATES['turn']['g+/g+']
IC_Phe_Bend = PHE_TEMPLATES['bend']['canonical']
IC_Phe_Bend_gminus_gminus = PHE_TEMPLATES['bend']['g-/g-']
IC_Phe_Bend_gminus_trans = PHE_TEMPLATES['bend']['g-/t']
IC_Phe_Bend_gminus_gplus = PHE_TEMPLATES['bend']['g-/g+']
IC_Phe_Bend_trans_gminus = PHE_TEMPLATES['bend']['t/g-']
IC_Phe_Bend_trans_trans = PHE_TEMPLATES['bend']['t/t']
IC_Phe_Bend_trans_gplus = PHE_TEMPLATES['bend']['t/g+']
IC_Phe_Bend_gplus_gminus = PHE_TEMPLATES['bend']['g+/g-']
IC_Phe_Bend_gplus_trans = PHE_TEMPLATES['bend']['g+/t']
IC_Phe_Bend_gplus_gplus = PHE_TEMPLATES['bend']['g+/g+']
IC_Phe_Coil = PHE_TEMPLATES['coil']['canonical']
IC_Phe_Coil_gminus_gminus = PHE_TEMPLATES['coil']['g-/g-']
IC_Phe_Coil_gminus_trans = PHE_TEMPLATES['coil']['g-/t']
IC_Phe_Coil_gminus_gplus = PHE_TEMPLATES['coil']['g-/g+']
IC_Phe_Coil_trans_gminus = PHE_TEMPLATES['coil']['t/g-']
IC_Phe_Coil_trans_trans = PHE_TEMPLATES['coil']['t/t']
IC_Phe_Coil_trans_gplus = PHE_TEMPLATES['coil']['t/g+']
IC_Phe_Coil_gplus_gminus = PHE_TEMPLATES['coil']['g+/g-']
IC_Phe_Coil_gplus_trans = PHE_TEMPLATES['coil']['g+/t']
IC_Phe_Coil_gplus_gplus = PHE_TEMPLATES['coil']['g+/g+']
IC_Phe_CisPeptide = PHE_TEMPLATES['cis-peptide-bond']['canonical']
IC_Phe_CisPeptide_gminus_gminus = PHE_TEMPLATES['cis-peptide-bond']['g-/g-']
IC_Phe_CisPeptide_gminus_trans = PHE_TEMPLATES['cis-peptide-bond']['g-/t']
IC_Phe_CisPeptide_gminus_gplus = PHE_TEMPLATES['cis-peptide-bond']['g-/g+']
IC_Phe_CisPeptide_trans_gminus = PHE_TEMPLATES['cis-peptide-bond']['t/g-']
IC_Phe_CisPeptide_trans_trans = PHE_TEMPLATES['cis-peptide-bond']['t/t']
IC_Phe_CisPeptide_trans_gplus = PHE_TEMPLATES['cis-peptide-bond']['t/g+']
IC_Phe_CisPeptide_gplus_gminus = PHE_TEMPLATES['cis-peptide-bond']['g+/g-']
IC_Phe_CisPeptide_gplus_trans = PHE_TEMPLATES['cis-peptide-bond']['g+/t']
IC_Phe_CisPeptide_gplus_gplus = PHE_TEMPLATES['cis-peptide-bond']['g+/g+']

__all__ = [
    'IC_Phe_Bend',
    'IC_Phe_Bend_gminus_gminus',
    'IC_Phe_Bend_gminus_gplus',
    'IC_Phe_Bend_gminus_trans',
    'IC_Phe_Bend_gplus_gminus',
    'IC_Phe_Bend_gplus_gplus',
    'IC_Phe_Bend_gplus_trans',
    'IC_Phe_Bend_trans_gminus',
    'IC_Phe_Bend_trans_gplus',
    'IC_Phe_Bend_trans_trans',
    'IC_Phe_Bridge',
    'IC_Phe_Bridge_gminus_gminus',
    'IC_Phe_Bridge_gminus_gplus',
    'IC_Phe_Bridge_gminus_trans',
    'IC_Phe_Bridge_gplus_gminus',
    'IC_Phe_Bridge_gplus_gplus',
    'IC_Phe_Bridge_gplus_trans',
    'IC_Phe_Bridge_trans_gminus',
    'IC_Phe_Bridge_trans_gplus',
    'IC_Phe_Bridge_trans_trans',
    'IC_Phe_CisPeptide',
    'IC_Phe_CisPeptide_gminus_gminus',
    'IC_Phe_CisPeptide_gminus_gplus',
    'IC_Phe_CisPeptide_gminus_trans',
    'IC_Phe_CisPeptide_gplus_gminus',
    'IC_Phe_CisPeptide_gplus_gplus',
    'IC_Phe_CisPeptide_gplus_trans',
    'IC_Phe_CisPeptide_trans_gminus',
    'IC_Phe_CisPeptide_trans_gplus',
    'IC_Phe_CisPeptide_trans_trans',
    'IC_Phe_Coil',
    'IC_Phe_Coil_gminus_gminus',
    'IC_Phe_Coil_gminus_gplus',
    'IC_Phe_Coil_gminus_trans',
    'IC_Phe_Coil_gplus_gminus',
    'IC_Phe_Coil_gplus_gplus',
    'IC_Phe_Coil_gplus_trans',
    'IC_Phe_Coil_trans_gminus',
    'IC_Phe_Coil_trans_gplus',
    'IC_Phe_Coil_trans_trans',
    'IC_Phe_Helix310',
    'IC_Phe_Helix310_gminus_gminus',
    'IC_Phe_Helix310_gminus_gplus',
    'IC_Phe_Helix310_gminus_trans',
    'IC_Phe_Helix310_gplus_gminus',
    'IC_Phe_Helix310_gplus_gplus',
    'IC_Phe_Helix310_gplus_trans',
    'IC_Phe_Helix310_trans_gminus',
    'IC_Phe_Helix310_trans_gplus',
    'IC_Phe_Helix310_trans_trans',
    'IC_Phe_HelixAlpha',
    'IC_Phe_HelixAlpha_gminus_gminus',
    'IC_Phe_HelixAlpha_gminus_gplus',
    'IC_Phe_HelixAlpha_gminus_trans',
    'IC_Phe_HelixAlpha_gplus_gminus',
    'IC_Phe_HelixAlpha_gplus_gplus',
    'IC_Phe_HelixAlpha_gplus_trans',
    'IC_Phe_HelixAlpha_trans_gminus',
    'IC_Phe_HelixAlpha_trans_gplus',
    'IC_Phe_HelixAlpha_trans_trans',
    'IC_Phe_HelixPPII',
    'IC_Phe_HelixPPII_gminus_gminus',
    'IC_Phe_HelixPPII_gminus_gplus',
    'IC_Phe_HelixPPII_gminus_trans',
    'IC_Phe_HelixPPII_gplus_gminus',
    'IC_Phe_HelixPPII_gplus_gplus',
    'IC_Phe_HelixPPII_gplus_trans',
    'IC_Phe_HelixPPII_trans_gminus',
    'IC_Phe_HelixPPII_trans_gplus',
    'IC_Phe_HelixPPII_trans_trans',
    'IC_Phe_HelixPi',
    'IC_Phe_HelixPi_gminus_gminus',
    'IC_Phe_HelixPi_gminus_gplus',
    'IC_Phe_HelixPi_gminus_trans',
    'IC_Phe_HelixPi_gplus_gminus',
    'IC_Phe_HelixPi_gplus_gplus',
    'IC_Phe_HelixPi_gplus_trans',
    'IC_Phe_HelixPi_trans_gminus',
    'IC_Phe_HelixPi_trans_gplus',
    'IC_Phe_HelixPi_trans_trans',
    'IC_Phe_Strand',
    'IC_Phe_StrandAntiParallel',
    'IC_Phe_StrandAntiParallel_gminus_gminus',
    'IC_Phe_StrandAntiParallel_gminus_gplus',
    'IC_Phe_StrandAntiParallel_gminus_trans',
    'IC_Phe_StrandAntiParallel_gplus_gminus',
    'IC_Phe_StrandAntiParallel_gplus_gplus',
    'IC_Phe_StrandAntiParallel_gplus_trans',
    'IC_Phe_StrandAntiParallel_trans_gminus',
    'IC_Phe_StrandAntiParallel_trans_gplus',
    'IC_Phe_StrandAntiParallel_trans_trans',
    'IC_Phe_StrandParallel',
    'IC_Phe_StrandParallel_gminus_gminus',
    'IC_Phe_StrandParallel_gminus_gplus',
    'IC_Phe_StrandParallel_gminus_trans',
    'IC_Phe_StrandParallel_gplus_gminus',
    'IC_Phe_StrandParallel_gplus_gplus',
    'IC_Phe_StrandParallel_gplus_trans',
    'IC_Phe_StrandParallel_trans_gminus',
    'IC_Phe_StrandParallel_trans_gplus',
    'IC_Phe_StrandParallel_trans_trans',
    'IC_Phe_Strand_gminus_gminus',
    'IC_Phe_Strand_gminus_gplus',
    'IC_Phe_Strand_gminus_trans',
    'IC_Phe_Strand_gplus_gminus',
    'IC_Phe_Strand_gplus_gplus',
    'IC_Phe_Strand_gplus_trans',
    'IC_Phe_Strand_trans_gminus',
    'IC_Phe_Strand_trans_gplus',
    'IC_Phe_Strand_trans_trans',
    'IC_Phe_Turn',
    'IC_Phe_Turn_gminus_gminus',
    'IC_Phe_Turn_gminus_gplus',
    'IC_Phe_Turn_gminus_trans',
    'IC_Phe_Turn_gplus_gminus',
    'IC_Phe_Turn_gplus_gplus',
    'IC_Phe_Turn_gplus_trans',
    'IC_Phe_Turn_trans_gminus',
    'IC_Phe_Turn_trans_gplus',
    'IC_Phe_Turn_trans_trans',
]
