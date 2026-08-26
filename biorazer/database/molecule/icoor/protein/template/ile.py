# -*- coding: utf-8 -*-
"""Ideal {secondary-structure x rotamer} {names[aa]} {var} ``InternalCoord`` templates.

Each template is the ideal heavy-atom topology for {aa} anchored at
{{N, CA, C}}, with the side chain grown to ``rotamer`` and the
``phi``/``psi``/``omega`` of ``ss`` carried on the instance.  See
:mod:`biorazer.database.internal_coord_template` and
:mod:`._builder` for the build rule and the per-conformer caveat.
"""

from biorazer.database.molecule.icoor.protein.template._builder import make_residue_templates

RESN = "ILE"
#: {ss: {rotamer: :class:`~biorazer.structure.objects.InternalCoord`}}
ILE_TEMPLATES = make_residue_templates(RESN)

# === GENERATED IC_* named templates - do not edit below ===
# Regenerate with: /opt/envs/BioRazer/bin/python scripts/generate_internal_coord_template_named.py

IC_Ile_HelixAlpha = ILE_TEMPLATES['alpha-helix']['canonical']
IC_Ile_HelixAlpha_gminus_gminus = ILE_TEMPLATES['alpha-helix']['g-/g-']
IC_Ile_HelixAlpha_gminus_trans = ILE_TEMPLATES['alpha-helix']['g-/t']
IC_Ile_HelixAlpha_gminus_gplus = ILE_TEMPLATES['alpha-helix']['g-/g+']
IC_Ile_HelixAlpha_trans_gminus = ILE_TEMPLATES['alpha-helix']['t/g-']
IC_Ile_HelixAlpha_trans_trans = ILE_TEMPLATES['alpha-helix']['t/t']
IC_Ile_HelixAlpha_trans_gplus = ILE_TEMPLATES['alpha-helix']['t/g+']
IC_Ile_HelixAlpha_gplus_gminus = ILE_TEMPLATES['alpha-helix']['g+/g-']
IC_Ile_HelixAlpha_gplus_trans = ILE_TEMPLATES['alpha-helix']['g+/t']
IC_Ile_HelixAlpha_gplus_gplus = ILE_TEMPLATES['alpha-helix']['g+/g+']
IC_Ile_Helix310 = ILE_TEMPLATES['3-10-helix']['canonical']
IC_Ile_Helix310_gminus_gminus = ILE_TEMPLATES['3-10-helix']['g-/g-']
IC_Ile_Helix310_gminus_trans = ILE_TEMPLATES['3-10-helix']['g-/t']
IC_Ile_Helix310_gminus_gplus = ILE_TEMPLATES['3-10-helix']['g-/g+']
IC_Ile_Helix310_trans_gminus = ILE_TEMPLATES['3-10-helix']['t/g-']
IC_Ile_Helix310_trans_trans = ILE_TEMPLATES['3-10-helix']['t/t']
IC_Ile_Helix310_trans_gplus = ILE_TEMPLATES['3-10-helix']['t/g+']
IC_Ile_Helix310_gplus_gminus = ILE_TEMPLATES['3-10-helix']['g+/g-']
IC_Ile_Helix310_gplus_trans = ILE_TEMPLATES['3-10-helix']['g+/t']
IC_Ile_Helix310_gplus_gplus = ILE_TEMPLATES['3-10-helix']['g+/g+']
IC_Ile_HelixPi = ILE_TEMPLATES['pi-helix']['canonical']
IC_Ile_HelixPi_gminus_gminus = ILE_TEMPLATES['pi-helix']['g-/g-']
IC_Ile_HelixPi_gminus_trans = ILE_TEMPLATES['pi-helix']['g-/t']
IC_Ile_HelixPi_gminus_gplus = ILE_TEMPLATES['pi-helix']['g-/g+']
IC_Ile_HelixPi_trans_gminus = ILE_TEMPLATES['pi-helix']['t/g-']
IC_Ile_HelixPi_trans_trans = ILE_TEMPLATES['pi-helix']['t/t']
IC_Ile_HelixPi_trans_gplus = ILE_TEMPLATES['pi-helix']['t/g+']
IC_Ile_HelixPi_gplus_gminus = ILE_TEMPLATES['pi-helix']['g+/g-']
IC_Ile_HelixPi_gplus_trans = ILE_TEMPLATES['pi-helix']['g+/t']
IC_Ile_HelixPi_gplus_gplus = ILE_TEMPLATES['pi-helix']['g+/g+']
IC_Ile_HelixPPII = ILE_TEMPLATES['polyproline-II']['canonical']
IC_Ile_HelixPPII_gminus_gminus = ILE_TEMPLATES['polyproline-II']['g-/g-']
IC_Ile_HelixPPII_gminus_trans = ILE_TEMPLATES['polyproline-II']['g-/t']
IC_Ile_HelixPPII_gminus_gplus = ILE_TEMPLATES['polyproline-II']['g-/g+']
IC_Ile_HelixPPII_trans_gminus = ILE_TEMPLATES['polyproline-II']['t/g-']
IC_Ile_HelixPPII_trans_trans = ILE_TEMPLATES['polyproline-II']['t/t']
IC_Ile_HelixPPII_trans_gplus = ILE_TEMPLATES['polyproline-II']['t/g+']
IC_Ile_HelixPPII_gplus_gminus = ILE_TEMPLATES['polyproline-II']['g+/g-']
IC_Ile_HelixPPII_gplus_trans = ILE_TEMPLATES['polyproline-II']['g+/t']
IC_Ile_HelixPPII_gplus_gplus = ILE_TEMPLATES['polyproline-II']['g+/g+']
IC_Ile_Strand = ILE_TEMPLATES['beta-strand']['canonical']
IC_Ile_Strand_gminus_gminus = ILE_TEMPLATES['beta-strand']['g-/g-']
IC_Ile_Strand_gminus_trans = ILE_TEMPLATES['beta-strand']['g-/t']
IC_Ile_Strand_gminus_gplus = ILE_TEMPLATES['beta-strand']['g-/g+']
IC_Ile_Strand_trans_gminus = ILE_TEMPLATES['beta-strand']['t/g-']
IC_Ile_Strand_trans_trans = ILE_TEMPLATES['beta-strand']['t/t']
IC_Ile_Strand_trans_gplus = ILE_TEMPLATES['beta-strand']['t/g+']
IC_Ile_Strand_gplus_gminus = ILE_TEMPLATES['beta-strand']['g+/g-']
IC_Ile_Strand_gplus_trans = ILE_TEMPLATES['beta-strand']['g+/t']
IC_Ile_Strand_gplus_gplus = ILE_TEMPLATES['beta-strand']['g+/g+']
IC_Ile_StrandParallel = ILE_TEMPLATES['parallel-beta-strand']['canonical']
IC_Ile_StrandParallel_gminus_gminus = ILE_TEMPLATES['parallel-beta-strand']['g-/g-']
IC_Ile_StrandParallel_gminus_trans = ILE_TEMPLATES['parallel-beta-strand']['g-/t']
IC_Ile_StrandParallel_gminus_gplus = ILE_TEMPLATES['parallel-beta-strand']['g-/g+']
IC_Ile_StrandParallel_trans_gminus = ILE_TEMPLATES['parallel-beta-strand']['t/g-']
IC_Ile_StrandParallel_trans_trans = ILE_TEMPLATES['parallel-beta-strand']['t/t']
IC_Ile_StrandParallel_trans_gplus = ILE_TEMPLATES['parallel-beta-strand']['t/g+']
IC_Ile_StrandParallel_gplus_gminus = ILE_TEMPLATES['parallel-beta-strand']['g+/g-']
IC_Ile_StrandParallel_gplus_trans = ILE_TEMPLATES['parallel-beta-strand']['g+/t']
IC_Ile_StrandParallel_gplus_gplus = ILE_TEMPLATES['parallel-beta-strand']['g+/g+']
IC_Ile_StrandAntiParallel = ILE_TEMPLATES['antiparallel-beta-strand']['canonical']
IC_Ile_StrandAntiParallel_gminus_gminus = ILE_TEMPLATES['antiparallel-beta-strand']['g-/g-']
IC_Ile_StrandAntiParallel_gminus_trans = ILE_TEMPLATES['antiparallel-beta-strand']['g-/t']
IC_Ile_StrandAntiParallel_gminus_gplus = ILE_TEMPLATES['antiparallel-beta-strand']['g-/g+']
IC_Ile_StrandAntiParallel_trans_gminus = ILE_TEMPLATES['antiparallel-beta-strand']['t/g-']
IC_Ile_StrandAntiParallel_trans_trans = ILE_TEMPLATES['antiparallel-beta-strand']['t/t']
IC_Ile_StrandAntiParallel_trans_gplus = ILE_TEMPLATES['antiparallel-beta-strand']['t/g+']
IC_Ile_StrandAntiParallel_gplus_gminus = ILE_TEMPLATES['antiparallel-beta-strand']['g+/g-']
IC_Ile_StrandAntiParallel_gplus_trans = ILE_TEMPLATES['antiparallel-beta-strand']['g+/t']
IC_Ile_StrandAntiParallel_gplus_gplus = ILE_TEMPLATES['antiparallel-beta-strand']['g+/g+']
IC_Ile_Bridge = ILE_TEMPLATES['beta-bridge']['canonical']
IC_Ile_Bridge_gminus_gminus = ILE_TEMPLATES['beta-bridge']['g-/g-']
IC_Ile_Bridge_gminus_trans = ILE_TEMPLATES['beta-bridge']['g-/t']
IC_Ile_Bridge_gminus_gplus = ILE_TEMPLATES['beta-bridge']['g-/g+']
IC_Ile_Bridge_trans_gminus = ILE_TEMPLATES['beta-bridge']['t/g-']
IC_Ile_Bridge_trans_trans = ILE_TEMPLATES['beta-bridge']['t/t']
IC_Ile_Bridge_trans_gplus = ILE_TEMPLATES['beta-bridge']['t/g+']
IC_Ile_Bridge_gplus_gminus = ILE_TEMPLATES['beta-bridge']['g+/g-']
IC_Ile_Bridge_gplus_trans = ILE_TEMPLATES['beta-bridge']['g+/t']
IC_Ile_Bridge_gplus_gplus = ILE_TEMPLATES['beta-bridge']['g+/g+']
IC_Ile_Turn = ILE_TEMPLATES['turn']['canonical']
IC_Ile_Turn_gminus_gminus = ILE_TEMPLATES['turn']['g-/g-']
IC_Ile_Turn_gminus_trans = ILE_TEMPLATES['turn']['g-/t']
IC_Ile_Turn_gminus_gplus = ILE_TEMPLATES['turn']['g-/g+']
IC_Ile_Turn_trans_gminus = ILE_TEMPLATES['turn']['t/g-']
IC_Ile_Turn_trans_trans = ILE_TEMPLATES['turn']['t/t']
IC_Ile_Turn_trans_gplus = ILE_TEMPLATES['turn']['t/g+']
IC_Ile_Turn_gplus_gminus = ILE_TEMPLATES['turn']['g+/g-']
IC_Ile_Turn_gplus_trans = ILE_TEMPLATES['turn']['g+/t']
IC_Ile_Turn_gplus_gplus = ILE_TEMPLATES['turn']['g+/g+']
IC_Ile_Bend = ILE_TEMPLATES['bend']['canonical']
IC_Ile_Bend_gminus_gminus = ILE_TEMPLATES['bend']['g-/g-']
IC_Ile_Bend_gminus_trans = ILE_TEMPLATES['bend']['g-/t']
IC_Ile_Bend_gminus_gplus = ILE_TEMPLATES['bend']['g-/g+']
IC_Ile_Bend_trans_gminus = ILE_TEMPLATES['bend']['t/g-']
IC_Ile_Bend_trans_trans = ILE_TEMPLATES['bend']['t/t']
IC_Ile_Bend_trans_gplus = ILE_TEMPLATES['bend']['t/g+']
IC_Ile_Bend_gplus_gminus = ILE_TEMPLATES['bend']['g+/g-']
IC_Ile_Bend_gplus_trans = ILE_TEMPLATES['bend']['g+/t']
IC_Ile_Bend_gplus_gplus = ILE_TEMPLATES['bend']['g+/g+']
IC_Ile_Coil = ILE_TEMPLATES['coil']['canonical']
IC_Ile_Coil_gminus_gminus = ILE_TEMPLATES['coil']['g-/g-']
IC_Ile_Coil_gminus_trans = ILE_TEMPLATES['coil']['g-/t']
IC_Ile_Coil_gminus_gplus = ILE_TEMPLATES['coil']['g-/g+']
IC_Ile_Coil_trans_gminus = ILE_TEMPLATES['coil']['t/g-']
IC_Ile_Coil_trans_trans = ILE_TEMPLATES['coil']['t/t']
IC_Ile_Coil_trans_gplus = ILE_TEMPLATES['coil']['t/g+']
IC_Ile_Coil_gplus_gminus = ILE_TEMPLATES['coil']['g+/g-']
IC_Ile_Coil_gplus_trans = ILE_TEMPLATES['coil']['g+/t']
IC_Ile_Coil_gplus_gplus = ILE_TEMPLATES['coil']['g+/g+']
IC_Ile_CisPeptide = ILE_TEMPLATES['cis-peptide-bond']['canonical']
IC_Ile_CisPeptide_gminus_gminus = ILE_TEMPLATES['cis-peptide-bond']['g-/g-']
IC_Ile_CisPeptide_gminus_trans = ILE_TEMPLATES['cis-peptide-bond']['g-/t']
IC_Ile_CisPeptide_gminus_gplus = ILE_TEMPLATES['cis-peptide-bond']['g-/g+']
IC_Ile_CisPeptide_trans_gminus = ILE_TEMPLATES['cis-peptide-bond']['t/g-']
IC_Ile_CisPeptide_trans_trans = ILE_TEMPLATES['cis-peptide-bond']['t/t']
IC_Ile_CisPeptide_trans_gplus = ILE_TEMPLATES['cis-peptide-bond']['t/g+']
IC_Ile_CisPeptide_gplus_gminus = ILE_TEMPLATES['cis-peptide-bond']['g+/g-']
IC_Ile_CisPeptide_gplus_trans = ILE_TEMPLATES['cis-peptide-bond']['g+/t']
IC_Ile_CisPeptide_gplus_gplus = ILE_TEMPLATES['cis-peptide-bond']['g+/g+']

__all__ = [
    'IC_Ile_Bend',
    'IC_Ile_Bend_gminus_gminus',
    'IC_Ile_Bend_gminus_gplus',
    'IC_Ile_Bend_gminus_trans',
    'IC_Ile_Bend_gplus_gminus',
    'IC_Ile_Bend_gplus_gplus',
    'IC_Ile_Bend_gplus_trans',
    'IC_Ile_Bend_trans_gminus',
    'IC_Ile_Bend_trans_gplus',
    'IC_Ile_Bend_trans_trans',
    'IC_Ile_Bridge',
    'IC_Ile_Bridge_gminus_gminus',
    'IC_Ile_Bridge_gminus_gplus',
    'IC_Ile_Bridge_gminus_trans',
    'IC_Ile_Bridge_gplus_gminus',
    'IC_Ile_Bridge_gplus_gplus',
    'IC_Ile_Bridge_gplus_trans',
    'IC_Ile_Bridge_trans_gminus',
    'IC_Ile_Bridge_trans_gplus',
    'IC_Ile_Bridge_trans_trans',
    'IC_Ile_CisPeptide',
    'IC_Ile_CisPeptide_gminus_gminus',
    'IC_Ile_CisPeptide_gminus_gplus',
    'IC_Ile_CisPeptide_gminus_trans',
    'IC_Ile_CisPeptide_gplus_gminus',
    'IC_Ile_CisPeptide_gplus_gplus',
    'IC_Ile_CisPeptide_gplus_trans',
    'IC_Ile_CisPeptide_trans_gminus',
    'IC_Ile_CisPeptide_trans_gplus',
    'IC_Ile_CisPeptide_trans_trans',
    'IC_Ile_Coil',
    'IC_Ile_Coil_gminus_gminus',
    'IC_Ile_Coil_gminus_gplus',
    'IC_Ile_Coil_gminus_trans',
    'IC_Ile_Coil_gplus_gminus',
    'IC_Ile_Coil_gplus_gplus',
    'IC_Ile_Coil_gplus_trans',
    'IC_Ile_Coil_trans_gminus',
    'IC_Ile_Coil_trans_gplus',
    'IC_Ile_Coil_trans_trans',
    'IC_Ile_Helix310',
    'IC_Ile_Helix310_gminus_gminus',
    'IC_Ile_Helix310_gminus_gplus',
    'IC_Ile_Helix310_gminus_trans',
    'IC_Ile_Helix310_gplus_gminus',
    'IC_Ile_Helix310_gplus_gplus',
    'IC_Ile_Helix310_gplus_trans',
    'IC_Ile_Helix310_trans_gminus',
    'IC_Ile_Helix310_trans_gplus',
    'IC_Ile_Helix310_trans_trans',
    'IC_Ile_HelixAlpha',
    'IC_Ile_HelixAlpha_gminus_gminus',
    'IC_Ile_HelixAlpha_gminus_gplus',
    'IC_Ile_HelixAlpha_gminus_trans',
    'IC_Ile_HelixAlpha_gplus_gminus',
    'IC_Ile_HelixAlpha_gplus_gplus',
    'IC_Ile_HelixAlpha_gplus_trans',
    'IC_Ile_HelixAlpha_trans_gminus',
    'IC_Ile_HelixAlpha_trans_gplus',
    'IC_Ile_HelixAlpha_trans_trans',
    'IC_Ile_HelixPPII',
    'IC_Ile_HelixPPII_gminus_gminus',
    'IC_Ile_HelixPPII_gminus_gplus',
    'IC_Ile_HelixPPII_gminus_trans',
    'IC_Ile_HelixPPII_gplus_gminus',
    'IC_Ile_HelixPPII_gplus_gplus',
    'IC_Ile_HelixPPII_gplus_trans',
    'IC_Ile_HelixPPII_trans_gminus',
    'IC_Ile_HelixPPII_trans_gplus',
    'IC_Ile_HelixPPII_trans_trans',
    'IC_Ile_HelixPi',
    'IC_Ile_HelixPi_gminus_gminus',
    'IC_Ile_HelixPi_gminus_gplus',
    'IC_Ile_HelixPi_gminus_trans',
    'IC_Ile_HelixPi_gplus_gminus',
    'IC_Ile_HelixPi_gplus_gplus',
    'IC_Ile_HelixPi_gplus_trans',
    'IC_Ile_HelixPi_trans_gminus',
    'IC_Ile_HelixPi_trans_gplus',
    'IC_Ile_HelixPi_trans_trans',
    'IC_Ile_Strand',
    'IC_Ile_StrandAntiParallel',
    'IC_Ile_StrandAntiParallel_gminus_gminus',
    'IC_Ile_StrandAntiParallel_gminus_gplus',
    'IC_Ile_StrandAntiParallel_gminus_trans',
    'IC_Ile_StrandAntiParallel_gplus_gminus',
    'IC_Ile_StrandAntiParallel_gplus_gplus',
    'IC_Ile_StrandAntiParallel_gplus_trans',
    'IC_Ile_StrandAntiParallel_trans_gminus',
    'IC_Ile_StrandAntiParallel_trans_gplus',
    'IC_Ile_StrandAntiParallel_trans_trans',
    'IC_Ile_StrandParallel',
    'IC_Ile_StrandParallel_gminus_gminus',
    'IC_Ile_StrandParallel_gminus_gplus',
    'IC_Ile_StrandParallel_gminus_trans',
    'IC_Ile_StrandParallel_gplus_gminus',
    'IC_Ile_StrandParallel_gplus_gplus',
    'IC_Ile_StrandParallel_gplus_trans',
    'IC_Ile_StrandParallel_trans_gminus',
    'IC_Ile_StrandParallel_trans_gplus',
    'IC_Ile_StrandParallel_trans_trans',
    'IC_Ile_Strand_gminus_gminus',
    'IC_Ile_Strand_gminus_gplus',
    'IC_Ile_Strand_gminus_trans',
    'IC_Ile_Strand_gplus_gminus',
    'IC_Ile_Strand_gplus_gplus',
    'IC_Ile_Strand_gplus_trans',
    'IC_Ile_Strand_trans_gminus',
    'IC_Ile_Strand_trans_gplus',
    'IC_Ile_Strand_trans_trans',
    'IC_Ile_Turn',
    'IC_Ile_Turn_gminus_gminus',
    'IC_Ile_Turn_gminus_gplus',
    'IC_Ile_Turn_gminus_trans',
    'IC_Ile_Turn_gplus_gminus',
    'IC_Ile_Turn_gplus_gplus',
    'IC_Ile_Turn_gplus_trans',
    'IC_Ile_Turn_trans_gminus',
    'IC_Ile_Turn_trans_gplus',
    'IC_Ile_Turn_trans_trans',
]
