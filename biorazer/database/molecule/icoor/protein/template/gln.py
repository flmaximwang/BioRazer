# -*- coding: utf-8 -*-
"""Ideal {secondary-structure x rotamer} {names[aa]} {var} ``InternalCoord`` templates.

Each template is the ideal heavy-atom topology for {aa} anchored at
{{N, CA, C}}, with the side chain grown to ``rotamer`` and the
``phi``/``psi``/``omega`` of ``ss`` carried on the instance.  See
:mod:`biorazer.database.internal_coord_template` and
:mod:`._builder` for the build rule and the per-conformer caveat.
"""

from biorazer.database.molecule.icoor.protein.template._builder import make_residue_templates

RESN = "GLN"
#: {ss: {rotamer: :class:`~biorazer.structure.objects.InternalCoord`}}
GLN_TEMPLATES = make_residue_templates(RESN)

# === GENERATED IC_* named templates - do not edit below ===
# Regenerate with: /opt/envs/BioRazer/bin/python scripts/generate_internal_coord_template_named.py

IC_Gln_HelixAlpha = GLN_TEMPLATES['alpha-helix']['canonical']
IC_Gln_HelixAlpha_gminus_gminus = GLN_TEMPLATES['alpha-helix']['g-/g-']
IC_Gln_HelixAlpha_gminus_trans = GLN_TEMPLATES['alpha-helix']['g-/t']
IC_Gln_HelixAlpha_gminus_gplus = GLN_TEMPLATES['alpha-helix']['g-/g+']
IC_Gln_HelixAlpha_trans_gminus = GLN_TEMPLATES['alpha-helix']['t/g-']
IC_Gln_HelixAlpha_trans_trans = GLN_TEMPLATES['alpha-helix']['t/t']
IC_Gln_HelixAlpha_trans_gplus = GLN_TEMPLATES['alpha-helix']['t/g+']
IC_Gln_HelixAlpha_gplus_gminus = GLN_TEMPLATES['alpha-helix']['g+/g-']
IC_Gln_HelixAlpha_gplus_trans = GLN_TEMPLATES['alpha-helix']['g+/t']
IC_Gln_HelixAlpha_gplus_gplus = GLN_TEMPLATES['alpha-helix']['g+/g+']
IC_Gln_Helix310 = GLN_TEMPLATES['3-10-helix']['canonical']
IC_Gln_Helix310_gminus_gminus = GLN_TEMPLATES['3-10-helix']['g-/g-']
IC_Gln_Helix310_gminus_trans = GLN_TEMPLATES['3-10-helix']['g-/t']
IC_Gln_Helix310_gminus_gplus = GLN_TEMPLATES['3-10-helix']['g-/g+']
IC_Gln_Helix310_trans_gminus = GLN_TEMPLATES['3-10-helix']['t/g-']
IC_Gln_Helix310_trans_trans = GLN_TEMPLATES['3-10-helix']['t/t']
IC_Gln_Helix310_trans_gplus = GLN_TEMPLATES['3-10-helix']['t/g+']
IC_Gln_Helix310_gplus_gminus = GLN_TEMPLATES['3-10-helix']['g+/g-']
IC_Gln_Helix310_gplus_trans = GLN_TEMPLATES['3-10-helix']['g+/t']
IC_Gln_Helix310_gplus_gplus = GLN_TEMPLATES['3-10-helix']['g+/g+']
IC_Gln_HelixPi = GLN_TEMPLATES['pi-helix']['canonical']
IC_Gln_HelixPi_gminus_gminus = GLN_TEMPLATES['pi-helix']['g-/g-']
IC_Gln_HelixPi_gminus_trans = GLN_TEMPLATES['pi-helix']['g-/t']
IC_Gln_HelixPi_gminus_gplus = GLN_TEMPLATES['pi-helix']['g-/g+']
IC_Gln_HelixPi_trans_gminus = GLN_TEMPLATES['pi-helix']['t/g-']
IC_Gln_HelixPi_trans_trans = GLN_TEMPLATES['pi-helix']['t/t']
IC_Gln_HelixPi_trans_gplus = GLN_TEMPLATES['pi-helix']['t/g+']
IC_Gln_HelixPi_gplus_gminus = GLN_TEMPLATES['pi-helix']['g+/g-']
IC_Gln_HelixPi_gplus_trans = GLN_TEMPLATES['pi-helix']['g+/t']
IC_Gln_HelixPi_gplus_gplus = GLN_TEMPLATES['pi-helix']['g+/g+']
IC_Gln_HelixPPII = GLN_TEMPLATES['polyproline-II']['canonical']
IC_Gln_HelixPPII_gminus_gminus = GLN_TEMPLATES['polyproline-II']['g-/g-']
IC_Gln_HelixPPII_gminus_trans = GLN_TEMPLATES['polyproline-II']['g-/t']
IC_Gln_HelixPPII_gminus_gplus = GLN_TEMPLATES['polyproline-II']['g-/g+']
IC_Gln_HelixPPII_trans_gminus = GLN_TEMPLATES['polyproline-II']['t/g-']
IC_Gln_HelixPPII_trans_trans = GLN_TEMPLATES['polyproline-II']['t/t']
IC_Gln_HelixPPII_trans_gplus = GLN_TEMPLATES['polyproline-II']['t/g+']
IC_Gln_HelixPPII_gplus_gminus = GLN_TEMPLATES['polyproline-II']['g+/g-']
IC_Gln_HelixPPII_gplus_trans = GLN_TEMPLATES['polyproline-II']['g+/t']
IC_Gln_HelixPPII_gplus_gplus = GLN_TEMPLATES['polyproline-II']['g+/g+']
IC_Gln_Strand = GLN_TEMPLATES['beta-strand']['canonical']
IC_Gln_Strand_gminus_gminus = GLN_TEMPLATES['beta-strand']['g-/g-']
IC_Gln_Strand_gminus_trans = GLN_TEMPLATES['beta-strand']['g-/t']
IC_Gln_Strand_gminus_gplus = GLN_TEMPLATES['beta-strand']['g-/g+']
IC_Gln_Strand_trans_gminus = GLN_TEMPLATES['beta-strand']['t/g-']
IC_Gln_Strand_trans_trans = GLN_TEMPLATES['beta-strand']['t/t']
IC_Gln_Strand_trans_gplus = GLN_TEMPLATES['beta-strand']['t/g+']
IC_Gln_Strand_gplus_gminus = GLN_TEMPLATES['beta-strand']['g+/g-']
IC_Gln_Strand_gplus_trans = GLN_TEMPLATES['beta-strand']['g+/t']
IC_Gln_Strand_gplus_gplus = GLN_TEMPLATES['beta-strand']['g+/g+']
IC_Gln_StrandParallel = GLN_TEMPLATES['parallel-beta-strand']['canonical']
IC_Gln_StrandParallel_gminus_gminus = GLN_TEMPLATES['parallel-beta-strand']['g-/g-']
IC_Gln_StrandParallel_gminus_trans = GLN_TEMPLATES['parallel-beta-strand']['g-/t']
IC_Gln_StrandParallel_gminus_gplus = GLN_TEMPLATES['parallel-beta-strand']['g-/g+']
IC_Gln_StrandParallel_trans_gminus = GLN_TEMPLATES['parallel-beta-strand']['t/g-']
IC_Gln_StrandParallel_trans_trans = GLN_TEMPLATES['parallel-beta-strand']['t/t']
IC_Gln_StrandParallel_trans_gplus = GLN_TEMPLATES['parallel-beta-strand']['t/g+']
IC_Gln_StrandParallel_gplus_gminus = GLN_TEMPLATES['parallel-beta-strand']['g+/g-']
IC_Gln_StrandParallel_gplus_trans = GLN_TEMPLATES['parallel-beta-strand']['g+/t']
IC_Gln_StrandParallel_gplus_gplus = GLN_TEMPLATES['parallel-beta-strand']['g+/g+']
IC_Gln_StrandAntiParallel = GLN_TEMPLATES['antiparallel-beta-strand']['canonical']
IC_Gln_StrandAntiParallel_gminus_gminus = GLN_TEMPLATES['antiparallel-beta-strand']['g-/g-']
IC_Gln_StrandAntiParallel_gminus_trans = GLN_TEMPLATES['antiparallel-beta-strand']['g-/t']
IC_Gln_StrandAntiParallel_gminus_gplus = GLN_TEMPLATES['antiparallel-beta-strand']['g-/g+']
IC_Gln_StrandAntiParallel_trans_gminus = GLN_TEMPLATES['antiparallel-beta-strand']['t/g-']
IC_Gln_StrandAntiParallel_trans_trans = GLN_TEMPLATES['antiparallel-beta-strand']['t/t']
IC_Gln_StrandAntiParallel_trans_gplus = GLN_TEMPLATES['antiparallel-beta-strand']['t/g+']
IC_Gln_StrandAntiParallel_gplus_gminus = GLN_TEMPLATES['antiparallel-beta-strand']['g+/g-']
IC_Gln_StrandAntiParallel_gplus_trans = GLN_TEMPLATES['antiparallel-beta-strand']['g+/t']
IC_Gln_StrandAntiParallel_gplus_gplus = GLN_TEMPLATES['antiparallel-beta-strand']['g+/g+']
IC_Gln_Bridge = GLN_TEMPLATES['beta-bridge']['canonical']
IC_Gln_Bridge_gminus_gminus = GLN_TEMPLATES['beta-bridge']['g-/g-']
IC_Gln_Bridge_gminus_trans = GLN_TEMPLATES['beta-bridge']['g-/t']
IC_Gln_Bridge_gminus_gplus = GLN_TEMPLATES['beta-bridge']['g-/g+']
IC_Gln_Bridge_trans_gminus = GLN_TEMPLATES['beta-bridge']['t/g-']
IC_Gln_Bridge_trans_trans = GLN_TEMPLATES['beta-bridge']['t/t']
IC_Gln_Bridge_trans_gplus = GLN_TEMPLATES['beta-bridge']['t/g+']
IC_Gln_Bridge_gplus_gminus = GLN_TEMPLATES['beta-bridge']['g+/g-']
IC_Gln_Bridge_gplus_trans = GLN_TEMPLATES['beta-bridge']['g+/t']
IC_Gln_Bridge_gplus_gplus = GLN_TEMPLATES['beta-bridge']['g+/g+']
IC_Gln_Turn = GLN_TEMPLATES['turn']['canonical']
IC_Gln_Turn_gminus_gminus = GLN_TEMPLATES['turn']['g-/g-']
IC_Gln_Turn_gminus_trans = GLN_TEMPLATES['turn']['g-/t']
IC_Gln_Turn_gminus_gplus = GLN_TEMPLATES['turn']['g-/g+']
IC_Gln_Turn_trans_gminus = GLN_TEMPLATES['turn']['t/g-']
IC_Gln_Turn_trans_trans = GLN_TEMPLATES['turn']['t/t']
IC_Gln_Turn_trans_gplus = GLN_TEMPLATES['turn']['t/g+']
IC_Gln_Turn_gplus_gminus = GLN_TEMPLATES['turn']['g+/g-']
IC_Gln_Turn_gplus_trans = GLN_TEMPLATES['turn']['g+/t']
IC_Gln_Turn_gplus_gplus = GLN_TEMPLATES['turn']['g+/g+']
IC_Gln_Bend = GLN_TEMPLATES['bend']['canonical']
IC_Gln_Bend_gminus_gminus = GLN_TEMPLATES['bend']['g-/g-']
IC_Gln_Bend_gminus_trans = GLN_TEMPLATES['bend']['g-/t']
IC_Gln_Bend_gminus_gplus = GLN_TEMPLATES['bend']['g-/g+']
IC_Gln_Bend_trans_gminus = GLN_TEMPLATES['bend']['t/g-']
IC_Gln_Bend_trans_trans = GLN_TEMPLATES['bend']['t/t']
IC_Gln_Bend_trans_gplus = GLN_TEMPLATES['bend']['t/g+']
IC_Gln_Bend_gplus_gminus = GLN_TEMPLATES['bend']['g+/g-']
IC_Gln_Bend_gplus_trans = GLN_TEMPLATES['bend']['g+/t']
IC_Gln_Bend_gplus_gplus = GLN_TEMPLATES['bend']['g+/g+']
IC_Gln_Coil = GLN_TEMPLATES['coil']['canonical']
IC_Gln_Coil_gminus_gminus = GLN_TEMPLATES['coil']['g-/g-']
IC_Gln_Coil_gminus_trans = GLN_TEMPLATES['coil']['g-/t']
IC_Gln_Coil_gminus_gplus = GLN_TEMPLATES['coil']['g-/g+']
IC_Gln_Coil_trans_gminus = GLN_TEMPLATES['coil']['t/g-']
IC_Gln_Coil_trans_trans = GLN_TEMPLATES['coil']['t/t']
IC_Gln_Coil_trans_gplus = GLN_TEMPLATES['coil']['t/g+']
IC_Gln_Coil_gplus_gminus = GLN_TEMPLATES['coil']['g+/g-']
IC_Gln_Coil_gplus_trans = GLN_TEMPLATES['coil']['g+/t']
IC_Gln_Coil_gplus_gplus = GLN_TEMPLATES['coil']['g+/g+']
IC_Gln_CisPeptide = GLN_TEMPLATES['cis-peptide-bond']['canonical']
IC_Gln_CisPeptide_gminus_gminus = GLN_TEMPLATES['cis-peptide-bond']['g-/g-']
IC_Gln_CisPeptide_gminus_trans = GLN_TEMPLATES['cis-peptide-bond']['g-/t']
IC_Gln_CisPeptide_gminus_gplus = GLN_TEMPLATES['cis-peptide-bond']['g-/g+']
IC_Gln_CisPeptide_trans_gminus = GLN_TEMPLATES['cis-peptide-bond']['t/g-']
IC_Gln_CisPeptide_trans_trans = GLN_TEMPLATES['cis-peptide-bond']['t/t']
IC_Gln_CisPeptide_trans_gplus = GLN_TEMPLATES['cis-peptide-bond']['t/g+']
IC_Gln_CisPeptide_gplus_gminus = GLN_TEMPLATES['cis-peptide-bond']['g+/g-']
IC_Gln_CisPeptide_gplus_trans = GLN_TEMPLATES['cis-peptide-bond']['g+/t']
IC_Gln_CisPeptide_gplus_gplus = GLN_TEMPLATES['cis-peptide-bond']['g+/g+']

__all__ = [
    'IC_Gln_Bend',
    'IC_Gln_Bend_gminus_gminus',
    'IC_Gln_Bend_gminus_gplus',
    'IC_Gln_Bend_gminus_trans',
    'IC_Gln_Bend_gplus_gminus',
    'IC_Gln_Bend_gplus_gplus',
    'IC_Gln_Bend_gplus_trans',
    'IC_Gln_Bend_trans_gminus',
    'IC_Gln_Bend_trans_gplus',
    'IC_Gln_Bend_trans_trans',
    'IC_Gln_Bridge',
    'IC_Gln_Bridge_gminus_gminus',
    'IC_Gln_Bridge_gminus_gplus',
    'IC_Gln_Bridge_gminus_trans',
    'IC_Gln_Bridge_gplus_gminus',
    'IC_Gln_Bridge_gplus_gplus',
    'IC_Gln_Bridge_gplus_trans',
    'IC_Gln_Bridge_trans_gminus',
    'IC_Gln_Bridge_trans_gplus',
    'IC_Gln_Bridge_trans_trans',
    'IC_Gln_CisPeptide',
    'IC_Gln_CisPeptide_gminus_gminus',
    'IC_Gln_CisPeptide_gminus_gplus',
    'IC_Gln_CisPeptide_gminus_trans',
    'IC_Gln_CisPeptide_gplus_gminus',
    'IC_Gln_CisPeptide_gplus_gplus',
    'IC_Gln_CisPeptide_gplus_trans',
    'IC_Gln_CisPeptide_trans_gminus',
    'IC_Gln_CisPeptide_trans_gplus',
    'IC_Gln_CisPeptide_trans_trans',
    'IC_Gln_Coil',
    'IC_Gln_Coil_gminus_gminus',
    'IC_Gln_Coil_gminus_gplus',
    'IC_Gln_Coil_gminus_trans',
    'IC_Gln_Coil_gplus_gminus',
    'IC_Gln_Coil_gplus_gplus',
    'IC_Gln_Coil_gplus_trans',
    'IC_Gln_Coil_trans_gminus',
    'IC_Gln_Coil_trans_gplus',
    'IC_Gln_Coil_trans_trans',
    'IC_Gln_Helix310',
    'IC_Gln_Helix310_gminus_gminus',
    'IC_Gln_Helix310_gminus_gplus',
    'IC_Gln_Helix310_gminus_trans',
    'IC_Gln_Helix310_gplus_gminus',
    'IC_Gln_Helix310_gplus_gplus',
    'IC_Gln_Helix310_gplus_trans',
    'IC_Gln_Helix310_trans_gminus',
    'IC_Gln_Helix310_trans_gplus',
    'IC_Gln_Helix310_trans_trans',
    'IC_Gln_HelixAlpha',
    'IC_Gln_HelixAlpha_gminus_gminus',
    'IC_Gln_HelixAlpha_gminus_gplus',
    'IC_Gln_HelixAlpha_gminus_trans',
    'IC_Gln_HelixAlpha_gplus_gminus',
    'IC_Gln_HelixAlpha_gplus_gplus',
    'IC_Gln_HelixAlpha_gplus_trans',
    'IC_Gln_HelixAlpha_trans_gminus',
    'IC_Gln_HelixAlpha_trans_gplus',
    'IC_Gln_HelixAlpha_trans_trans',
    'IC_Gln_HelixPPII',
    'IC_Gln_HelixPPII_gminus_gminus',
    'IC_Gln_HelixPPII_gminus_gplus',
    'IC_Gln_HelixPPII_gminus_trans',
    'IC_Gln_HelixPPII_gplus_gminus',
    'IC_Gln_HelixPPII_gplus_gplus',
    'IC_Gln_HelixPPII_gplus_trans',
    'IC_Gln_HelixPPII_trans_gminus',
    'IC_Gln_HelixPPII_trans_gplus',
    'IC_Gln_HelixPPII_trans_trans',
    'IC_Gln_HelixPi',
    'IC_Gln_HelixPi_gminus_gminus',
    'IC_Gln_HelixPi_gminus_gplus',
    'IC_Gln_HelixPi_gminus_trans',
    'IC_Gln_HelixPi_gplus_gminus',
    'IC_Gln_HelixPi_gplus_gplus',
    'IC_Gln_HelixPi_gplus_trans',
    'IC_Gln_HelixPi_trans_gminus',
    'IC_Gln_HelixPi_trans_gplus',
    'IC_Gln_HelixPi_trans_trans',
    'IC_Gln_Strand',
    'IC_Gln_StrandAntiParallel',
    'IC_Gln_StrandAntiParallel_gminus_gminus',
    'IC_Gln_StrandAntiParallel_gminus_gplus',
    'IC_Gln_StrandAntiParallel_gminus_trans',
    'IC_Gln_StrandAntiParallel_gplus_gminus',
    'IC_Gln_StrandAntiParallel_gplus_gplus',
    'IC_Gln_StrandAntiParallel_gplus_trans',
    'IC_Gln_StrandAntiParallel_trans_gminus',
    'IC_Gln_StrandAntiParallel_trans_gplus',
    'IC_Gln_StrandAntiParallel_trans_trans',
    'IC_Gln_StrandParallel',
    'IC_Gln_StrandParallel_gminus_gminus',
    'IC_Gln_StrandParallel_gminus_gplus',
    'IC_Gln_StrandParallel_gminus_trans',
    'IC_Gln_StrandParallel_gplus_gminus',
    'IC_Gln_StrandParallel_gplus_gplus',
    'IC_Gln_StrandParallel_gplus_trans',
    'IC_Gln_StrandParallel_trans_gminus',
    'IC_Gln_StrandParallel_trans_gplus',
    'IC_Gln_StrandParallel_trans_trans',
    'IC_Gln_Strand_gminus_gminus',
    'IC_Gln_Strand_gminus_gplus',
    'IC_Gln_Strand_gminus_trans',
    'IC_Gln_Strand_gplus_gminus',
    'IC_Gln_Strand_gplus_gplus',
    'IC_Gln_Strand_gplus_trans',
    'IC_Gln_Strand_trans_gminus',
    'IC_Gln_Strand_trans_gplus',
    'IC_Gln_Strand_trans_trans',
    'IC_Gln_Turn',
    'IC_Gln_Turn_gminus_gminus',
    'IC_Gln_Turn_gminus_gplus',
    'IC_Gln_Turn_gminus_trans',
    'IC_Gln_Turn_gplus_gminus',
    'IC_Gln_Turn_gplus_gplus',
    'IC_Gln_Turn_gplus_trans',
    'IC_Gln_Turn_trans_gminus',
    'IC_Gln_Turn_trans_gplus',
    'IC_Gln_Turn_trans_trans',
]
