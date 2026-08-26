# -*- coding: utf-8 -*-
"""Ideal {secondary-structure x rotamer} {names[aa]} {var} ``InternalCoord`` templates.

Each template is the ideal heavy-atom topology for {aa} anchored at
{{N, CA, C}}, with the side chain grown to ``rotamer`` and the
``phi``/``psi``/``omega`` of ``ss`` carried on the instance.  See
:mod:`biorazer.database.molecule.icoor.protein.template` and
:mod:`._builder` for the build rule and the per-conformer caveat.
"""

from biorazer.database.molecule.icoor.protein.template._builder import make_residue_templates

RESN = "GLU"
#: {ss: {rotamer: :class:`~biorazer.structure.objects.InternalCoord`}}
GLU_TEMPLATES = make_residue_templates(RESN)

# === GENERATED IC_* named templates - do not edit below ===
# Regenerate with: /opt/envs/BioRazer/bin/python scripts/generate_internal_coord_template_named.py

IC_Glu_HelixAlpha = GLU_TEMPLATES['alpha-helix']['canonical']
IC_Glu_HelixAlpha_gminus_gminus = GLU_TEMPLATES['alpha-helix']['g-/g-']
IC_Glu_HelixAlpha_gminus_trans = GLU_TEMPLATES['alpha-helix']['g-/t']
IC_Glu_HelixAlpha_gminus_gplus = GLU_TEMPLATES['alpha-helix']['g-/g+']
IC_Glu_HelixAlpha_trans_gminus = GLU_TEMPLATES['alpha-helix']['t/g-']
IC_Glu_HelixAlpha_trans_trans = GLU_TEMPLATES['alpha-helix']['t/t']
IC_Glu_HelixAlpha_trans_gplus = GLU_TEMPLATES['alpha-helix']['t/g+']
IC_Glu_HelixAlpha_gplus_gminus = GLU_TEMPLATES['alpha-helix']['g+/g-']
IC_Glu_HelixAlpha_gplus_trans = GLU_TEMPLATES['alpha-helix']['g+/t']
IC_Glu_HelixAlpha_gplus_gplus = GLU_TEMPLATES['alpha-helix']['g+/g+']
IC_Glu_Helix310 = GLU_TEMPLATES['3-10-helix']['canonical']
IC_Glu_Helix310_gminus_gminus = GLU_TEMPLATES['3-10-helix']['g-/g-']
IC_Glu_Helix310_gminus_trans = GLU_TEMPLATES['3-10-helix']['g-/t']
IC_Glu_Helix310_gminus_gplus = GLU_TEMPLATES['3-10-helix']['g-/g+']
IC_Glu_Helix310_trans_gminus = GLU_TEMPLATES['3-10-helix']['t/g-']
IC_Glu_Helix310_trans_trans = GLU_TEMPLATES['3-10-helix']['t/t']
IC_Glu_Helix310_trans_gplus = GLU_TEMPLATES['3-10-helix']['t/g+']
IC_Glu_Helix310_gplus_gminus = GLU_TEMPLATES['3-10-helix']['g+/g-']
IC_Glu_Helix310_gplus_trans = GLU_TEMPLATES['3-10-helix']['g+/t']
IC_Glu_Helix310_gplus_gplus = GLU_TEMPLATES['3-10-helix']['g+/g+']
IC_Glu_HelixPi = GLU_TEMPLATES['pi-helix']['canonical']
IC_Glu_HelixPi_gminus_gminus = GLU_TEMPLATES['pi-helix']['g-/g-']
IC_Glu_HelixPi_gminus_trans = GLU_TEMPLATES['pi-helix']['g-/t']
IC_Glu_HelixPi_gminus_gplus = GLU_TEMPLATES['pi-helix']['g-/g+']
IC_Glu_HelixPi_trans_gminus = GLU_TEMPLATES['pi-helix']['t/g-']
IC_Glu_HelixPi_trans_trans = GLU_TEMPLATES['pi-helix']['t/t']
IC_Glu_HelixPi_trans_gplus = GLU_TEMPLATES['pi-helix']['t/g+']
IC_Glu_HelixPi_gplus_gminus = GLU_TEMPLATES['pi-helix']['g+/g-']
IC_Glu_HelixPi_gplus_trans = GLU_TEMPLATES['pi-helix']['g+/t']
IC_Glu_HelixPi_gplus_gplus = GLU_TEMPLATES['pi-helix']['g+/g+']
IC_Glu_HelixPPII = GLU_TEMPLATES['polyproline-II']['canonical']
IC_Glu_HelixPPII_gminus_gminus = GLU_TEMPLATES['polyproline-II']['g-/g-']
IC_Glu_HelixPPII_gminus_trans = GLU_TEMPLATES['polyproline-II']['g-/t']
IC_Glu_HelixPPII_gminus_gplus = GLU_TEMPLATES['polyproline-II']['g-/g+']
IC_Glu_HelixPPII_trans_gminus = GLU_TEMPLATES['polyproline-II']['t/g-']
IC_Glu_HelixPPII_trans_trans = GLU_TEMPLATES['polyproline-II']['t/t']
IC_Glu_HelixPPII_trans_gplus = GLU_TEMPLATES['polyproline-II']['t/g+']
IC_Glu_HelixPPII_gplus_gminus = GLU_TEMPLATES['polyproline-II']['g+/g-']
IC_Glu_HelixPPII_gplus_trans = GLU_TEMPLATES['polyproline-II']['g+/t']
IC_Glu_HelixPPII_gplus_gplus = GLU_TEMPLATES['polyproline-II']['g+/g+']
IC_Glu_Strand = GLU_TEMPLATES['beta-strand']['canonical']
IC_Glu_Strand_gminus_gminus = GLU_TEMPLATES['beta-strand']['g-/g-']
IC_Glu_Strand_gminus_trans = GLU_TEMPLATES['beta-strand']['g-/t']
IC_Glu_Strand_gminus_gplus = GLU_TEMPLATES['beta-strand']['g-/g+']
IC_Glu_Strand_trans_gminus = GLU_TEMPLATES['beta-strand']['t/g-']
IC_Glu_Strand_trans_trans = GLU_TEMPLATES['beta-strand']['t/t']
IC_Glu_Strand_trans_gplus = GLU_TEMPLATES['beta-strand']['t/g+']
IC_Glu_Strand_gplus_gminus = GLU_TEMPLATES['beta-strand']['g+/g-']
IC_Glu_Strand_gplus_trans = GLU_TEMPLATES['beta-strand']['g+/t']
IC_Glu_Strand_gplus_gplus = GLU_TEMPLATES['beta-strand']['g+/g+']
IC_Glu_StrandParallel = GLU_TEMPLATES['parallel-beta-strand']['canonical']
IC_Glu_StrandParallel_gminus_gminus = GLU_TEMPLATES['parallel-beta-strand']['g-/g-']
IC_Glu_StrandParallel_gminus_trans = GLU_TEMPLATES['parallel-beta-strand']['g-/t']
IC_Glu_StrandParallel_gminus_gplus = GLU_TEMPLATES['parallel-beta-strand']['g-/g+']
IC_Glu_StrandParallel_trans_gminus = GLU_TEMPLATES['parallel-beta-strand']['t/g-']
IC_Glu_StrandParallel_trans_trans = GLU_TEMPLATES['parallel-beta-strand']['t/t']
IC_Glu_StrandParallel_trans_gplus = GLU_TEMPLATES['parallel-beta-strand']['t/g+']
IC_Glu_StrandParallel_gplus_gminus = GLU_TEMPLATES['parallel-beta-strand']['g+/g-']
IC_Glu_StrandParallel_gplus_trans = GLU_TEMPLATES['parallel-beta-strand']['g+/t']
IC_Glu_StrandParallel_gplus_gplus = GLU_TEMPLATES['parallel-beta-strand']['g+/g+']
IC_Glu_StrandAntiParallel = GLU_TEMPLATES['antiparallel-beta-strand']['canonical']
IC_Glu_StrandAntiParallel_gminus_gminus = GLU_TEMPLATES['antiparallel-beta-strand']['g-/g-']
IC_Glu_StrandAntiParallel_gminus_trans = GLU_TEMPLATES['antiparallel-beta-strand']['g-/t']
IC_Glu_StrandAntiParallel_gminus_gplus = GLU_TEMPLATES['antiparallel-beta-strand']['g-/g+']
IC_Glu_StrandAntiParallel_trans_gminus = GLU_TEMPLATES['antiparallel-beta-strand']['t/g-']
IC_Glu_StrandAntiParallel_trans_trans = GLU_TEMPLATES['antiparallel-beta-strand']['t/t']
IC_Glu_StrandAntiParallel_trans_gplus = GLU_TEMPLATES['antiparallel-beta-strand']['t/g+']
IC_Glu_StrandAntiParallel_gplus_gminus = GLU_TEMPLATES['antiparallel-beta-strand']['g+/g-']
IC_Glu_StrandAntiParallel_gplus_trans = GLU_TEMPLATES['antiparallel-beta-strand']['g+/t']
IC_Glu_StrandAntiParallel_gplus_gplus = GLU_TEMPLATES['antiparallel-beta-strand']['g+/g+']
IC_Glu_Bridge = GLU_TEMPLATES['beta-bridge']['canonical']
IC_Glu_Bridge_gminus_gminus = GLU_TEMPLATES['beta-bridge']['g-/g-']
IC_Glu_Bridge_gminus_trans = GLU_TEMPLATES['beta-bridge']['g-/t']
IC_Glu_Bridge_gminus_gplus = GLU_TEMPLATES['beta-bridge']['g-/g+']
IC_Glu_Bridge_trans_gminus = GLU_TEMPLATES['beta-bridge']['t/g-']
IC_Glu_Bridge_trans_trans = GLU_TEMPLATES['beta-bridge']['t/t']
IC_Glu_Bridge_trans_gplus = GLU_TEMPLATES['beta-bridge']['t/g+']
IC_Glu_Bridge_gplus_gminus = GLU_TEMPLATES['beta-bridge']['g+/g-']
IC_Glu_Bridge_gplus_trans = GLU_TEMPLATES['beta-bridge']['g+/t']
IC_Glu_Bridge_gplus_gplus = GLU_TEMPLATES['beta-bridge']['g+/g+']
IC_Glu_Turn = GLU_TEMPLATES['turn']['canonical']
IC_Glu_Turn_gminus_gminus = GLU_TEMPLATES['turn']['g-/g-']
IC_Glu_Turn_gminus_trans = GLU_TEMPLATES['turn']['g-/t']
IC_Glu_Turn_gminus_gplus = GLU_TEMPLATES['turn']['g-/g+']
IC_Glu_Turn_trans_gminus = GLU_TEMPLATES['turn']['t/g-']
IC_Glu_Turn_trans_trans = GLU_TEMPLATES['turn']['t/t']
IC_Glu_Turn_trans_gplus = GLU_TEMPLATES['turn']['t/g+']
IC_Glu_Turn_gplus_gminus = GLU_TEMPLATES['turn']['g+/g-']
IC_Glu_Turn_gplus_trans = GLU_TEMPLATES['turn']['g+/t']
IC_Glu_Turn_gplus_gplus = GLU_TEMPLATES['turn']['g+/g+']
IC_Glu_Bend = GLU_TEMPLATES['bend']['canonical']
IC_Glu_Bend_gminus_gminus = GLU_TEMPLATES['bend']['g-/g-']
IC_Glu_Bend_gminus_trans = GLU_TEMPLATES['bend']['g-/t']
IC_Glu_Bend_gminus_gplus = GLU_TEMPLATES['bend']['g-/g+']
IC_Glu_Bend_trans_gminus = GLU_TEMPLATES['bend']['t/g-']
IC_Glu_Bend_trans_trans = GLU_TEMPLATES['bend']['t/t']
IC_Glu_Bend_trans_gplus = GLU_TEMPLATES['bend']['t/g+']
IC_Glu_Bend_gplus_gminus = GLU_TEMPLATES['bend']['g+/g-']
IC_Glu_Bend_gplus_trans = GLU_TEMPLATES['bend']['g+/t']
IC_Glu_Bend_gplus_gplus = GLU_TEMPLATES['bend']['g+/g+']
IC_Glu_Coil = GLU_TEMPLATES['coil']['canonical']
IC_Glu_Coil_gminus_gminus = GLU_TEMPLATES['coil']['g-/g-']
IC_Glu_Coil_gminus_trans = GLU_TEMPLATES['coil']['g-/t']
IC_Glu_Coil_gminus_gplus = GLU_TEMPLATES['coil']['g-/g+']
IC_Glu_Coil_trans_gminus = GLU_TEMPLATES['coil']['t/g-']
IC_Glu_Coil_trans_trans = GLU_TEMPLATES['coil']['t/t']
IC_Glu_Coil_trans_gplus = GLU_TEMPLATES['coil']['t/g+']
IC_Glu_Coil_gplus_gminus = GLU_TEMPLATES['coil']['g+/g-']
IC_Glu_Coil_gplus_trans = GLU_TEMPLATES['coil']['g+/t']
IC_Glu_Coil_gplus_gplus = GLU_TEMPLATES['coil']['g+/g+']
IC_Glu_CisPeptide = GLU_TEMPLATES['cis-peptide-bond']['canonical']
IC_Glu_CisPeptide_gminus_gminus = GLU_TEMPLATES['cis-peptide-bond']['g-/g-']
IC_Glu_CisPeptide_gminus_trans = GLU_TEMPLATES['cis-peptide-bond']['g-/t']
IC_Glu_CisPeptide_gminus_gplus = GLU_TEMPLATES['cis-peptide-bond']['g-/g+']
IC_Glu_CisPeptide_trans_gminus = GLU_TEMPLATES['cis-peptide-bond']['t/g-']
IC_Glu_CisPeptide_trans_trans = GLU_TEMPLATES['cis-peptide-bond']['t/t']
IC_Glu_CisPeptide_trans_gplus = GLU_TEMPLATES['cis-peptide-bond']['t/g+']
IC_Glu_CisPeptide_gplus_gminus = GLU_TEMPLATES['cis-peptide-bond']['g+/g-']
IC_Glu_CisPeptide_gplus_trans = GLU_TEMPLATES['cis-peptide-bond']['g+/t']
IC_Glu_CisPeptide_gplus_gplus = GLU_TEMPLATES['cis-peptide-bond']['g+/g+']

__all__ = [
    'IC_Glu_Bend',
    'IC_Glu_Bend_gminus_gminus',
    'IC_Glu_Bend_gminus_gplus',
    'IC_Glu_Bend_gminus_trans',
    'IC_Glu_Bend_gplus_gminus',
    'IC_Glu_Bend_gplus_gplus',
    'IC_Glu_Bend_gplus_trans',
    'IC_Glu_Bend_trans_gminus',
    'IC_Glu_Bend_trans_gplus',
    'IC_Glu_Bend_trans_trans',
    'IC_Glu_Bridge',
    'IC_Glu_Bridge_gminus_gminus',
    'IC_Glu_Bridge_gminus_gplus',
    'IC_Glu_Bridge_gminus_trans',
    'IC_Glu_Bridge_gplus_gminus',
    'IC_Glu_Bridge_gplus_gplus',
    'IC_Glu_Bridge_gplus_trans',
    'IC_Glu_Bridge_trans_gminus',
    'IC_Glu_Bridge_trans_gplus',
    'IC_Glu_Bridge_trans_trans',
    'IC_Glu_CisPeptide',
    'IC_Glu_CisPeptide_gminus_gminus',
    'IC_Glu_CisPeptide_gminus_gplus',
    'IC_Glu_CisPeptide_gminus_trans',
    'IC_Glu_CisPeptide_gplus_gminus',
    'IC_Glu_CisPeptide_gplus_gplus',
    'IC_Glu_CisPeptide_gplus_trans',
    'IC_Glu_CisPeptide_trans_gminus',
    'IC_Glu_CisPeptide_trans_gplus',
    'IC_Glu_CisPeptide_trans_trans',
    'IC_Glu_Coil',
    'IC_Glu_Coil_gminus_gminus',
    'IC_Glu_Coil_gminus_gplus',
    'IC_Glu_Coil_gminus_trans',
    'IC_Glu_Coil_gplus_gminus',
    'IC_Glu_Coil_gplus_gplus',
    'IC_Glu_Coil_gplus_trans',
    'IC_Glu_Coil_trans_gminus',
    'IC_Glu_Coil_trans_gplus',
    'IC_Glu_Coil_trans_trans',
    'IC_Glu_Helix310',
    'IC_Glu_Helix310_gminus_gminus',
    'IC_Glu_Helix310_gminus_gplus',
    'IC_Glu_Helix310_gminus_trans',
    'IC_Glu_Helix310_gplus_gminus',
    'IC_Glu_Helix310_gplus_gplus',
    'IC_Glu_Helix310_gplus_trans',
    'IC_Glu_Helix310_trans_gminus',
    'IC_Glu_Helix310_trans_gplus',
    'IC_Glu_Helix310_trans_trans',
    'IC_Glu_HelixAlpha',
    'IC_Glu_HelixAlpha_gminus_gminus',
    'IC_Glu_HelixAlpha_gminus_gplus',
    'IC_Glu_HelixAlpha_gminus_trans',
    'IC_Glu_HelixAlpha_gplus_gminus',
    'IC_Glu_HelixAlpha_gplus_gplus',
    'IC_Glu_HelixAlpha_gplus_trans',
    'IC_Glu_HelixAlpha_trans_gminus',
    'IC_Glu_HelixAlpha_trans_gplus',
    'IC_Glu_HelixAlpha_trans_trans',
    'IC_Glu_HelixPPII',
    'IC_Glu_HelixPPII_gminus_gminus',
    'IC_Glu_HelixPPII_gminus_gplus',
    'IC_Glu_HelixPPII_gminus_trans',
    'IC_Glu_HelixPPII_gplus_gminus',
    'IC_Glu_HelixPPII_gplus_gplus',
    'IC_Glu_HelixPPII_gplus_trans',
    'IC_Glu_HelixPPII_trans_gminus',
    'IC_Glu_HelixPPII_trans_gplus',
    'IC_Glu_HelixPPII_trans_trans',
    'IC_Glu_HelixPi',
    'IC_Glu_HelixPi_gminus_gminus',
    'IC_Glu_HelixPi_gminus_gplus',
    'IC_Glu_HelixPi_gminus_trans',
    'IC_Glu_HelixPi_gplus_gminus',
    'IC_Glu_HelixPi_gplus_gplus',
    'IC_Glu_HelixPi_gplus_trans',
    'IC_Glu_HelixPi_trans_gminus',
    'IC_Glu_HelixPi_trans_gplus',
    'IC_Glu_HelixPi_trans_trans',
    'IC_Glu_Strand',
    'IC_Glu_StrandAntiParallel',
    'IC_Glu_StrandAntiParallel_gminus_gminus',
    'IC_Glu_StrandAntiParallel_gminus_gplus',
    'IC_Glu_StrandAntiParallel_gminus_trans',
    'IC_Glu_StrandAntiParallel_gplus_gminus',
    'IC_Glu_StrandAntiParallel_gplus_gplus',
    'IC_Glu_StrandAntiParallel_gplus_trans',
    'IC_Glu_StrandAntiParallel_trans_gminus',
    'IC_Glu_StrandAntiParallel_trans_gplus',
    'IC_Glu_StrandAntiParallel_trans_trans',
    'IC_Glu_StrandParallel',
    'IC_Glu_StrandParallel_gminus_gminus',
    'IC_Glu_StrandParallel_gminus_gplus',
    'IC_Glu_StrandParallel_gminus_trans',
    'IC_Glu_StrandParallel_gplus_gminus',
    'IC_Glu_StrandParallel_gplus_gplus',
    'IC_Glu_StrandParallel_gplus_trans',
    'IC_Glu_StrandParallel_trans_gminus',
    'IC_Glu_StrandParallel_trans_gplus',
    'IC_Glu_StrandParallel_trans_trans',
    'IC_Glu_Strand_gminus_gminus',
    'IC_Glu_Strand_gminus_gplus',
    'IC_Glu_Strand_gminus_trans',
    'IC_Glu_Strand_gplus_gminus',
    'IC_Glu_Strand_gplus_gplus',
    'IC_Glu_Strand_gplus_trans',
    'IC_Glu_Strand_trans_gminus',
    'IC_Glu_Strand_trans_gplus',
    'IC_Glu_Strand_trans_trans',
    'IC_Glu_Turn',
    'IC_Glu_Turn_gminus_gminus',
    'IC_Glu_Turn_gminus_gplus',
    'IC_Glu_Turn_gminus_trans',
    'IC_Glu_Turn_gplus_gminus',
    'IC_Glu_Turn_gplus_gplus',
    'IC_Glu_Turn_gplus_trans',
    'IC_Glu_Turn_trans_gminus',
    'IC_Glu_Turn_trans_gplus',
    'IC_Glu_Turn_trans_trans',
]
