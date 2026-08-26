# -*- coding: utf-8 -*-
"""Ideal {secondary-structure x rotamer} {names[aa]} {var} ``InternalCoord`` templates.

Each template is the ideal heavy-atom topology for {aa} anchored at
{{N, CA, C}}, with the side chain grown to ``rotamer`` and the
``phi``/``psi``/``omega`` of ``ss`` carried on the instance.  See
:mod:`biorazer.database.molecule.icoor.protein.template` and
:mod:`._builder` for the build rule and the per-conformer caveat.
"""

from biorazer.database.molecule.icoor.protein.template._builder import make_residue_templates

RESN = "ARG"
#: {ss: {rotamer: :class:`~biorazer.structure.objects.InternalCoord`}}
ARG_TEMPLATES = make_residue_templates(RESN)

# === GENERATED IC_* named templates - do not edit below ===
# Regenerate with: /opt/envs/BioRazer/bin/python scripts/generate_internal_coord_template_named.py

IC_Arg_HelixAlpha = ARG_TEMPLATES['alpha-helix']['canonical']
IC_Arg_HelixAlpha_gminus_gminus = ARG_TEMPLATES['alpha-helix']['g-/g-']
IC_Arg_HelixAlpha_gminus_trans = ARG_TEMPLATES['alpha-helix']['g-/t']
IC_Arg_HelixAlpha_gminus_gplus = ARG_TEMPLATES['alpha-helix']['g-/g+']
IC_Arg_HelixAlpha_trans_gminus = ARG_TEMPLATES['alpha-helix']['t/g-']
IC_Arg_HelixAlpha_trans_trans = ARG_TEMPLATES['alpha-helix']['t/t']
IC_Arg_HelixAlpha_trans_gplus = ARG_TEMPLATES['alpha-helix']['t/g+']
IC_Arg_HelixAlpha_gplus_gminus = ARG_TEMPLATES['alpha-helix']['g+/g-']
IC_Arg_HelixAlpha_gplus_trans = ARG_TEMPLATES['alpha-helix']['g+/t']
IC_Arg_HelixAlpha_gplus_gplus = ARG_TEMPLATES['alpha-helix']['g+/g+']
IC_Arg_Helix310 = ARG_TEMPLATES['3-10-helix']['canonical']
IC_Arg_Helix310_gminus_gminus = ARG_TEMPLATES['3-10-helix']['g-/g-']
IC_Arg_Helix310_gminus_trans = ARG_TEMPLATES['3-10-helix']['g-/t']
IC_Arg_Helix310_gminus_gplus = ARG_TEMPLATES['3-10-helix']['g-/g+']
IC_Arg_Helix310_trans_gminus = ARG_TEMPLATES['3-10-helix']['t/g-']
IC_Arg_Helix310_trans_trans = ARG_TEMPLATES['3-10-helix']['t/t']
IC_Arg_Helix310_trans_gplus = ARG_TEMPLATES['3-10-helix']['t/g+']
IC_Arg_Helix310_gplus_gminus = ARG_TEMPLATES['3-10-helix']['g+/g-']
IC_Arg_Helix310_gplus_trans = ARG_TEMPLATES['3-10-helix']['g+/t']
IC_Arg_Helix310_gplus_gplus = ARG_TEMPLATES['3-10-helix']['g+/g+']
IC_Arg_HelixPi = ARG_TEMPLATES['pi-helix']['canonical']
IC_Arg_HelixPi_gminus_gminus = ARG_TEMPLATES['pi-helix']['g-/g-']
IC_Arg_HelixPi_gminus_trans = ARG_TEMPLATES['pi-helix']['g-/t']
IC_Arg_HelixPi_gminus_gplus = ARG_TEMPLATES['pi-helix']['g-/g+']
IC_Arg_HelixPi_trans_gminus = ARG_TEMPLATES['pi-helix']['t/g-']
IC_Arg_HelixPi_trans_trans = ARG_TEMPLATES['pi-helix']['t/t']
IC_Arg_HelixPi_trans_gplus = ARG_TEMPLATES['pi-helix']['t/g+']
IC_Arg_HelixPi_gplus_gminus = ARG_TEMPLATES['pi-helix']['g+/g-']
IC_Arg_HelixPi_gplus_trans = ARG_TEMPLATES['pi-helix']['g+/t']
IC_Arg_HelixPi_gplus_gplus = ARG_TEMPLATES['pi-helix']['g+/g+']
IC_Arg_HelixPPII = ARG_TEMPLATES['polyproline-II']['canonical']
IC_Arg_HelixPPII_gminus_gminus = ARG_TEMPLATES['polyproline-II']['g-/g-']
IC_Arg_HelixPPII_gminus_trans = ARG_TEMPLATES['polyproline-II']['g-/t']
IC_Arg_HelixPPII_gminus_gplus = ARG_TEMPLATES['polyproline-II']['g-/g+']
IC_Arg_HelixPPII_trans_gminus = ARG_TEMPLATES['polyproline-II']['t/g-']
IC_Arg_HelixPPII_trans_trans = ARG_TEMPLATES['polyproline-II']['t/t']
IC_Arg_HelixPPII_trans_gplus = ARG_TEMPLATES['polyproline-II']['t/g+']
IC_Arg_HelixPPII_gplus_gminus = ARG_TEMPLATES['polyproline-II']['g+/g-']
IC_Arg_HelixPPII_gplus_trans = ARG_TEMPLATES['polyproline-II']['g+/t']
IC_Arg_HelixPPII_gplus_gplus = ARG_TEMPLATES['polyproline-II']['g+/g+']
IC_Arg_Strand = ARG_TEMPLATES['beta-strand']['canonical']
IC_Arg_Strand_gminus_gminus = ARG_TEMPLATES['beta-strand']['g-/g-']
IC_Arg_Strand_gminus_trans = ARG_TEMPLATES['beta-strand']['g-/t']
IC_Arg_Strand_gminus_gplus = ARG_TEMPLATES['beta-strand']['g-/g+']
IC_Arg_Strand_trans_gminus = ARG_TEMPLATES['beta-strand']['t/g-']
IC_Arg_Strand_trans_trans = ARG_TEMPLATES['beta-strand']['t/t']
IC_Arg_Strand_trans_gplus = ARG_TEMPLATES['beta-strand']['t/g+']
IC_Arg_Strand_gplus_gminus = ARG_TEMPLATES['beta-strand']['g+/g-']
IC_Arg_Strand_gplus_trans = ARG_TEMPLATES['beta-strand']['g+/t']
IC_Arg_Strand_gplus_gplus = ARG_TEMPLATES['beta-strand']['g+/g+']
IC_Arg_StrandParallel = ARG_TEMPLATES['parallel-beta-strand']['canonical']
IC_Arg_StrandParallel_gminus_gminus = ARG_TEMPLATES['parallel-beta-strand']['g-/g-']
IC_Arg_StrandParallel_gminus_trans = ARG_TEMPLATES['parallel-beta-strand']['g-/t']
IC_Arg_StrandParallel_gminus_gplus = ARG_TEMPLATES['parallel-beta-strand']['g-/g+']
IC_Arg_StrandParallel_trans_gminus = ARG_TEMPLATES['parallel-beta-strand']['t/g-']
IC_Arg_StrandParallel_trans_trans = ARG_TEMPLATES['parallel-beta-strand']['t/t']
IC_Arg_StrandParallel_trans_gplus = ARG_TEMPLATES['parallel-beta-strand']['t/g+']
IC_Arg_StrandParallel_gplus_gminus = ARG_TEMPLATES['parallel-beta-strand']['g+/g-']
IC_Arg_StrandParallel_gplus_trans = ARG_TEMPLATES['parallel-beta-strand']['g+/t']
IC_Arg_StrandParallel_gplus_gplus = ARG_TEMPLATES['parallel-beta-strand']['g+/g+']
IC_Arg_StrandAntiParallel = ARG_TEMPLATES['antiparallel-beta-strand']['canonical']
IC_Arg_StrandAntiParallel_gminus_gminus = ARG_TEMPLATES['antiparallel-beta-strand']['g-/g-']
IC_Arg_StrandAntiParallel_gminus_trans = ARG_TEMPLATES['antiparallel-beta-strand']['g-/t']
IC_Arg_StrandAntiParallel_gminus_gplus = ARG_TEMPLATES['antiparallel-beta-strand']['g-/g+']
IC_Arg_StrandAntiParallel_trans_gminus = ARG_TEMPLATES['antiparallel-beta-strand']['t/g-']
IC_Arg_StrandAntiParallel_trans_trans = ARG_TEMPLATES['antiparallel-beta-strand']['t/t']
IC_Arg_StrandAntiParallel_trans_gplus = ARG_TEMPLATES['antiparallel-beta-strand']['t/g+']
IC_Arg_StrandAntiParallel_gplus_gminus = ARG_TEMPLATES['antiparallel-beta-strand']['g+/g-']
IC_Arg_StrandAntiParallel_gplus_trans = ARG_TEMPLATES['antiparallel-beta-strand']['g+/t']
IC_Arg_StrandAntiParallel_gplus_gplus = ARG_TEMPLATES['antiparallel-beta-strand']['g+/g+']
IC_Arg_Bridge = ARG_TEMPLATES['beta-bridge']['canonical']
IC_Arg_Bridge_gminus_gminus = ARG_TEMPLATES['beta-bridge']['g-/g-']
IC_Arg_Bridge_gminus_trans = ARG_TEMPLATES['beta-bridge']['g-/t']
IC_Arg_Bridge_gminus_gplus = ARG_TEMPLATES['beta-bridge']['g-/g+']
IC_Arg_Bridge_trans_gminus = ARG_TEMPLATES['beta-bridge']['t/g-']
IC_Arg_Bridge_trans_trans = ARG_TEMPLATES['beta-bridge']['t/t']
IC_Arg_Bridge_trans_gplus = ARG_TEMPLATES['beta-bridge']['t/g+']
IC_Arg_Bridge_gplus_gminus = ARG_TEMPLATES['beta-bridge']['g+/g-']
IC_Arg_Bridge_gplus_trans = ARG_TEMPLATES['beta-bridge']['g+/t']
IC_Arg_Bridge_gplus_gplus = ARG_TEMPLATES['beta-bridge']['g+/g+']
IC_Arg_Turn = ARG_TEMPLATES['turn']['canonical']
IC_Arg_Turn_gminus_gminus = ARG_TEMPLATES['turn']['g-/g-']
IC_Arg_Turn_gminus_trans = ARG_TEMPLATES['turn']['g-/t']
IC_Arg_Turn_gminus_gplus = ARG_TEMPLATES['turn']['g-/g+']
IC_Arg_Turn_trans_gminus = ARG_TEMPLATES['turn']['t/g-']
IC_Arg_Turn_trans_trans = ARG_TEMPLATES['turn']['t/t']
IC_Arg_Turn_trans_gplus = ARG_TEMPLATES['turn']['t/g+']
IC_Arg_Turn_gplus_gminus = ARG_TEMPLATES['turn']['g+/g-']
IC_Arg_Turn_gplus_trans = ARG_TEMPLATES['turn']['g+/t']
IC_Arg_Turn_gplus_gplus = ARG_TEMPLATES['turn']['g+/g+']
IC_Arg_Bend = ARG_TEMPLATES['bend']['canonical']
IC_Arg_Bend_gminus_gminus = ARG_TEMPLATES['bend']['g-/g-']
IC_Arg_Bend_gminus_trans = ARG_TEMPLATES['bend']['g-/t']
IC_Arg_Bend_gminus_gplus = ARG_TEMPLATES['bend']['g-/g+']
IC_Arg_Bend_trans_gminus = ARG_TEMPLATES['bend']['t/g-']
IC_Arg_Bend_trans_trans = ARG_TEMPLATES['bend']['t/t']
IC_Arg_Bend_trans_gplus = ARG_TEMPLATES['bend']['t/g+']
IC_Arg_Bend_gplus_gminus = ARG_TEMPLATES['bend']['g+/g-']
IC_Arg_Bend_gplus_trans = ARG_TEMPLATES['bend']['g+/t']
IC_Arg_Bend_gplus_gplus = ARG_TEMPLATES['bend']['g+/g+']
IC_Arg_Coil = ARG_TEMPLATES['coil']['canonical']
IC_Arg_Coil_gminus_gminus = ARG_TEMPLATES['coil']['g-/g-']
IC_Arg_Coil_gminus_trans = ARG_TEMPLATES['coil']['g-/t']
IC_Arg_Coil_gminus_gplus = ARG_TEMPLATES['coil']['g-/g+']
IC_Arg_Coil_trans_gminus = ARG_TEMPLATES['coil']['t/g-']
IC_Arg_Coil_trans_trans = ARG_TEMPLATES['coil']['t/t']
IC_Arg_Coil_trans_gplus = ARG_TEMPLATES['coil']['t/g+']
IC_Arg_Coil_gplus_gminus = ARG_TEMPLATES['coil']['g+/g-']
IC_Arg_Coil_gplus_trans = ARG_TEMPLATES['coil']['g+/t']
IC_Arg_Coil_gplus_gplus = ARG_TEMPLATES['coil']['g+/g+']
IC_Arg_CisPeptide = ARG_TEMPLATES['cis-peptide-bond']['canonical']
IC_Arg_CisPeptide_gminus_gminus = ARG_TEMPLATES['cis-peptide-bond']['g-/g-']
IC_Arg_CisPeptide_gminus_trans = ARG_TEMPLATES['cis-peptide-bond']['g-/t']
IC_Arg_CisPeptide_gminus_gplus = ARG_TEMPLATES['cis-peptide-bond']['g-/g+']
IC_Arg_CisPeptide_trans_gminus = ARG_TEMPLATES['cis-peptide-bond']['t/g-']
IC_Arg_CisPeptide_trans_trans = ARG_TEMPLATES['cis-peptide-bond']['t/t']
IC_Arg_CisPeptide_trans_gplus = ARG_TEMPLATES['cis-peptide-bond']['t/g+']
IC_Arg_CisPeptide_gplus_gminus = ARG_TEMPLATES['cis-peptide-bond']['g+/g-']
IC_Arg_CisPeptide_gplus_trans = ARG_TEMPLATES['cis-peptide-bond']['g+/t']
IC_Arg_CisPeptide_gplus_gplus = ARG_TEMPLATES['cis-peptide-bond']['g+/g+']

__all__ = [
    'IC_Arg_Bend',
    'IC_Arg_Bend_gminus_gminus',
    'IC_Arg_Bend_gminus_gplus',
    'IC_Arg_Bend_gminus_trans',
    'IC_Arg_Bend_gplus_gminus',
    'IC_Arg_Bend_gplus_gplus',
    'IC_Arg_Bend_gplus_trans',
    'IC_Arg_Bend_trans_gminus',
    'IC_Arg_Bend_trans_gplus',
    'IC_Arg_Bend_trans_trans',
    'IC_Arg_Bridge',
    'IC_Arg_Bridge_gminus_gminus',
    'IC_Arg_Bridge_gminus_gplus',
    'IC_Arg_Bridge_gminus_trans',
    'IC_Arg_Bridge_gplus_gminus',
    'IC_Arg_Bridge_gplus_gplus',
    'IC_Arg_Bridge_gplus_trans',
    'IC_Arg_Bridge_trans_gminus',
    'IC_Arg_Bridge_trans_gplus',
    'IC_Arg_Bridge_trans_trans',
    'IC_Arg_CisPeptide',
    'IC_Arg_CisPeptide_gminus_gminus',
    'IC_Arg_CisPeptide_gminus_gplus',
    'IC_Arg_CisPeptide_gminus_trans',
    'IC_Arg_CisPeptide_gplus_gminus',
    'IC_Arg_CisPeptide_gplus_gplus',
    'IC_Arg_CisPeptide_gplus_trans',
    'IC_Arg_CisPeptide_trans_gminus',
    'IC_Arg_CisPeptide_trans_gplus',
    'IC_Arg_CisPeptide_trans_trans',
    'IC_Arg_Coil',
    'IC_Arg_Coil_gminus_gminus',
    'IC_Arg_Coil_gminus_gplus',
    'IC_Arg_Coil_gminus_trans',
    'IC_Arg_Coil_gplus_gminus',
    'IC_Arg_Coil_gplus_gplus',
    'IC_Arg_Coil_gplus_trans',
    'IC_Arg_Coil_trans_gminus',
    'IC_Arg_Coil_trans_gplus',
    'IC_Arg_Coil_trans_trans',
    'IC_Arg_Helix310',
    'IC_Arg_Helix310_gminus_gminus',
    'IC_Arg_Helix310_gminus_gplus',
    'IC_Arg_Helix310_gminus_trans',
    'IC_Arg_Helix310_gplus_gminus',
    'IC_Arg_Helix310_gplus_gplus',
    'IC_Arg_Helix310_gplus_trans',
    'IC_Arg_Helix310_trans_gminus',
    'IC_Arg_Helix310_trans_gplus',
    'IC_Arg_Helix310_trans_trans',
    'IC_Arg_HelixAlpha',
    'IC_Arg_HelixAlpha_gminus_gminus',
    'IC_Arg_HelixAlpha_gminus_gplus',
    'IC_Arg_HelixAlpha_gminus_trans',
    'IC_Arg_HelixAlpha_gplus_gminus',
    'IC_Arg_HelixAlpha_gplus_gplus',
    'IC_Arg_HelixAlpha_gplus_trans',
    'IC_Arg_HelixAlpha_trans_gminus',
    'IC_Arg_HelixAlpha_trans_gplus',
    'IC_Arg_HelixAlpha_trans_trans',
    'IC_Arg_HelixPPII',
    'IC_Arg_HelixPPII_gminus_gminus',
    'IC_Arg_HelixPPII_gminus_gplus',
    'IC_Arg_HelixPPII_gminus_trans',
    'IC_Arg_HelixPPII_gplus_gminus',
    'IC_Arg_HelixPPII_gplus_gplus',
    'IC_Arg_HelixPPII_gplus_trans',
    'IC_Arg_HelixPPII_trans_gminus',
    'IC_Arg_HelixPPII_trans_gplus',
    'IC_Arg_HelixPPII_trans_trans',
    'IC_Arg_HelixPi',
    'IC_Arg_HelixPi_gminus_gminus',
    'IC_Arg_HelixPi_gminus_gplus',
    'IC_Arg_HelixPi_gminus_trans',
    'IC_Arg_HelixPi_gplus_gminus',
    'IC_Arg_HelixPi_gplus_gplus',
    'IC_Arg_HelixPi_gplus_trans',
    'IC_Arg_HelixPi_trans_gminus',
    'IC_Arg_HelixPi_trans_gplus',
    'IC_Arg_HelixPi_trans_trans',
    'IC_Arg_Strand',
    'IC_Arg_StrandAntiParallel',
    'IC_Arg_StrandAntiParallel_gminus_gminus',
    'IC_Arg_StrandAntiParallel_gminus_gplus',
    'IC_Arg_StrandAntiParallel_gminus_trans',
    'IC_Arg_StrandAntiParallel_gplus_gminus',
    'IC_Arg_StrandAntiParallel_gplus_gplus',
    'IC_Arg_StrandAntiParallel_gplus_trans',
    'IC_Arg_StrandAntiParallel_trans_gminus',
    'IC_Arg_StrandAntiParallel_trans_gplus',
    'IC_Arg_StrandAntiParallel_trans_trans',
    'IC_Arg_StrandParallel',
    'IC_Arg_StrandParallel_gminus_gminus',
    'IC_Arg_StrandParallel_gminus_gplus',
    'IC_Arg_StrandParallel_gminus_trans',
    'IC_Arg_StrandParallel_gplus_gminus',
    'IC_Arg_StrandParallel_gplus_gplus',
    'IC_Arg_StrandParallel_gplus_trans',
    'IC_Arg_StrandParallel_trans_gminus',
    'IC_Arg_StrandParallel_trans_gplus',
    'IC_Arg_StrandParallel_trans_trans',
    'IC_Arg_Strand_gminus_gminus',
    'IC_Arg_Strand_gminus_gplus',
    'IC_Arg_Strand_gminus_trans',
    'IC_Arg_Strand_gplus_gminus',
    'IC_Arg_Strand_gplus_gplus',
    'IC_Arg_Strand_gplus_trans',
    'IC_Arg_Strand_trans_gminus',
    'IC_Arg_Strand_trans_gplus',
    'IC_Arg_Strand_trans_trans',
    'IC_Arg_Turn',
    'IC_Arg_Turn_gminus_gminus',
    'IC_Arg_Turn_gminus_gplus',
    'IC_Arg_Turn_gminus_trans',
    'IC_Arg_Turn_gplus_gminus',
    'IC_Arg_Turn_gplus_gplus',
    'IC_Arg_Turn_gplus_trans',
    'IC_Arg_Turn_trans_gminus',
    'IC_Arg_Turn_trans_gplus',
    'IC_Arg_Turn_trans_trans',
]
