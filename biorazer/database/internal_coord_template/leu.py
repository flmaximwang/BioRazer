# -*- coding: utf-8 -*-
"""Ideal {secondary-structure x rotamer} {names[aa]} {var} ``InternalCoord`` templates.

Each template is the ideal heavy-atom topology for {aa} anchored at
{{N, CA, C}}, with the side chain grown to ``rotamer`` and the
``phi``/``psi``/``omega`` of ``ss`` carried on the instance.  See
:mod:`biorazer.database.internal_coord_template` and
:mod:`._builder` for the build rule and the per-conformer caveat.
"""

from biorazer.database.internal_coord_template._builder import make_residue_templates

RESN = "LEU"
#: {ss: {rotamer: :class:`~biorazer.structure.objects.InternalCoord`}}
LEU_TEMPLATES = make_residue_templates(RESN)

# === GENERATED IC_* named templates - do not edit below ===
# Regenerate with: /opt/envs/BioRazer/bin/python scripts/generate_internal_coord_template_named.py

IC_Leu_HelixAlpha = LEU_TEMPLATES['alpha-helix']['canonical']
IC_Leu_HelixAlpha_gminus_gminus = LEU_TEMPLATES['alpha-helix']['g-/g-']
IC_Leu_HelixAlpha_gminus_trans = LEU_TEMPLATES['alpha-helix']['g-/t']
IC_Leu_HelixAlpha_gminus_gplus = LEU_TEMPLATES['alpha-helix']['g-/g+']
IC_Leu_HelixAlpha_trans_gminus = LEU_TEMPLATES['alpha-helix']['t/g-']
IC_Leu_HelixAlpha_trans_trans = LEU_TEMPLATES['alpha-helix']['t/t']
IC_Leu_HelixAlpha_trans_gplus = LEU_TEMPLATES['alpha-helix']['t/g+']
IC_Leu_HelixAlpha_gplus_gminus = LEU_TEMPLATES['alpha-helix']['g+/g-']
IC_Leu_HelixAlpha_gplus_trans = LEU_TEMPLATES['alpha-helix']['g+/t']
IC_Leu_HelixAlpha_gplus_gplus = LEU_TEMPLATES['alpha-helix']['g+/g+']
IC_Leu_Helix310 = LEU_TEMPLATES['3-10-helix']['canonical']
IC_Leu_Helix310_gminus_gminus = LEU_TEMPLATES['3-10-helix']['g-/g-']
IC_Leu_Helix310_gminus_trans = LEU_TEMPLATES['3-10-helix']['g-/t']
IC_Leu_Helix310_gminus_gplus = LEU_TEMPLATES['3-10-helix']['g-/g+']
IC_Leu_Helix310_trans_gminus = LEU_TEMPLATES['3-10-helix']['t/g-']
IC_Leu_Helix310_trans_trans = LEU_TEMPLATES['3-10-helix']['t/t']
IC_Leu_Helix310_trans_gplus = LEU_TEMPLATES['3-10-helix']['t/g+']
IC_Leu_Helix310_gplus_gminus = LEU_TEMPLATES['3-10-helix']['g+/g-']
IC_Leu_Helix310_gplus_trans = LEU_TEMPLATES['3-10-helix']['g+/t']
IC_Leu_Helix310_gplus_gplus = LEU_TEMPLATES['3-10-helix']['g+/g+']
IC_Leu_HelixPi = LEU_TEMPLATES['pi-helix']['canonical']
IC_Leu_HelixPi_gminus_gminus = LEU_TEMPLATES['pi-helix']['g-/g-']
IC_Leu_HelixPi_gminus_trans = LEU_TEMPLATES['pi-helix']['g-/t']
IC_Leu_HelixPi_gminus_gplus = LEU_TEMPLATES['pi-helix']['g-/g+']
IC_Leu_HelixPi_trans_gminus = LEU_TEMPLATES['pi-helix']['t/g-']
IC_Leu_HelixPi_trans_trans = LEU_TEMPLATES['pi-helix']['t/t']
IC_Leu_HelixPi_trans_gplus = LEU_TEMPLATES['pi-helix']['t/g+']
IC_Leu_HelixPi_gplus_gminus = LEU_TEMPLATES['pi-helix']['g+/g-']
IC_Leu_HelixPi_gplus_trans = LEU_TEMPLATES['pi-helix']['g+/t']
IC_Leu_HelixPi_gplus_gplus = LEU_TEMPLATES['pi-helix']['g+/g+']
IC_Leu_HelixPPII = LEU_TEMPLATES['polyproline-II']['canonical']
IC_Leu_HelixPPII_gminus_gminus = LEU_TEMPLATES['polyproline-II']['g-/g-']
IC_Leu_HelixPPII_gminus_trans = LEU_TEMPLATES['polyproline-II']['g-/t']
IC_Leu_HelixPPII_gminus_gplus = LEU_TEMPLATES['polyproline-II']['g-/g+']
IC_Leu_HelixPPII_trans_gminus = LEU_TEMPLATES['polyproline-II']['t/g-']
IC_Leu_HelixPPII_trans_trans = LEU_TEMPLATES['polyproline-II']['t/t']
IC_Leu_HelixPPII_trans_gplus = LEU_TEMPLATES['polyproline-II']['t/g+']
IC_Leu_HelixPPII_gplus_gminus = LEU_TEMPLATES['polyproline-II']['g+/g-']
IC_Leu_HelixPPII_gplus_trans = LEU_TEMPLATES['polyproline-II']['g+/t']
IC_Leu_HelixPPII_gplus_gplus = LEU_TEMPLATES['polyproline-II']['g+/g+']
IC_Leu_Strand = LEU_TEMPLATES['beta-strand']['canonical']
IC_Leu_Strand_gminus_gminus = LEU_TEMPLATES['beta-strand']['g-/g-']
IC_Leu_Strand_gminus_trans = LEU_TEMPLATES['beta-strand']['g-/t']
IC_Leu_Strand_gminus_gplus = LEU_TEMPLATES['beta-strand']['g-/g+']
IC_Leu_Strand_trans_gminus = LEU_TEMPLATES['beta-strand']['t/g-']
IC_Leu_Strand_trans_trans = LEU_TEMPLATES['beta-strand']['t/t']
IC_Leu_Strand_trans_gplus = LEU_TEMPLATES['beta-strand']['t/g+']
IC_Leu_Strand_gplus_gminus = LEU_TEMPLATES['beta-strand']['g+/g-']
IC_Leu_Strand_gplus_trans = LEU_TEMPLATES['beta-strand']['g+/t']
IC_Leu_Strand_gplus_gplus = LEU_TEMPLATES['beta-strand']['g+/g+']
IC_Leu_StrandParallel = LEU_TEMPLATES['parallel-beta-strand']['canonical']
IC_Leu_StrandParallel_gminus_gminus = LEU_TEMPLATES['parallel-beta-strand']['g-/g-']
IC_Leu_StrandParallel_gminus_trans = LEU_TEMPLATES['parallel-beta-strand']['g-/t']
IC_Leu_StrandParallel_gminus_gplus = LEU_TEMPLATES['parallel-beta-strand']['g-/g+']
IC_Leu_StrandParallel_trans_gminus = LEU_TEMPLATES['parallel-beta-strand']['t/g-']
IC_Leu_StrandParallel_trans_trans = LEU_TEMPLATES['parallel-beta-strand']['t/t']
IC_Leu_StrandParallel_trans_gplus = LEU_TEMPLATES['parallel-beta-strand']['t/g+']
IC_Leu_StrandParallel_gplus_gminus = LEU_TEMPLATES['parallel-beta-strand']['g+/g-']
IC_Leu_StrandParallel_gplus_trans = LEU_TEMPLATES['parallel-beta-strand']['g+/t']
IC_Leu_StrandParallel_gplus_gplus = LEU_TEMPLATES['parallel-beta-strand']['g+/g+']
IC_Leu_StrandAntiParallel = LEU_TEMPLATES['antiparallel-beta-strand']['canonical']
IC_Leu_StrandAntiParallel_gminus_gminus = LEU_TEMPLATES['antiparallel-beta-strand']['g-/g-']
IC_Leu_StrandAntiParallel_gminus_trans = LEU_TEMPLATES['antiparallel-beta-strand']['g-/t']
IC_Leu_StrandAntiParallel_gminus_gplus = LEU_TEMPLATES['antiparallel-beta-strand']['g-/g+']
IC_Leu_StrandAntiParallel_trans_gminus = LEU_TEMPLATES['antiparallel-beta-strand']['t/g-']
IC_Leu_StrandAntiParallel_trans_trans = LEU_TEMPLATES['antiparallel-beta-strand']['t/t']
IC_Leu_StrandAntiParallel_trans_gplus = LEU_TEMPLATES['antiparallel-beta-strand']['t/g+']
IC_Leu_StrandAntiParallel_gplus_gminus = LEU_TEMPLATES['antiparallel-beta-strand']['g+/g-']
IC_Leu_StrandAntiParallel_gplus_trans = LEU_TEMPLATES['antiparallel-beta-strand']['g+/t']
IC_Leu_StrandAntiParallel_gplus_gplus = LEU_TEMPLATES['antiparallel-beta-strand']['g+/g+']
IC_Leu_Bridge = LEU_TEMPLATES['beta-bridge']['canonical']
IC_Leu_Bridge_gminus_gminus = LEU_TEMPLATES['beta-bridge']['g-/g-']
IC_Leu_Bridge_gminus_trans = LEU_TEMPLATES['beta-bridge']['g-/t']
IC_Leu_Bridge_gminus_gplus = LEU_TEMPLATES['beta-bridge']['g-/g+']
IC_Leu_Bridge_trans_gminus = LEU_TEMPLATES['beta-bridge']['t/g-']
IC_Leu_Bridge_trans_trans = LEU_TEMPLATES['beta-bridge']['t/t']
IC_Leu_Bridge_trans_gplus = LEU_TEMPLATES['beta-bridge']['t/g+']
IC_Leu_Bridge_gplus_gminus = LEU_TEMPLATES['beta-bridge']['g+/g-']
IC_Leu_Bridge_gplus_trans = LEU_TEMPLATES['beta-bridge']['g+/t']
IC_Leu_Bridge_gplus_gplus = LEU_TEMPLATES['beta-bridge']['g+/g+']
IC_Leu_Turn = LEU_TEMPLATES['turn']['canonical']
IC_Leu_Turn_gminus_gminus = LEU_TEMPLATES['turn']['g-/g-']
IC_Leu_Turn_gminus_trans = LEU_TEMPLATES['turn']['g-/t']
IC_Leu_Turn_gminus_gplus = LEU_TEMPLATES['turn']['g-/g+']
IC_Leu_Turn_trans_gminus = LEU_TEMPLATES['turn']['t/g-']
IC_Leu_Turn_trans_trans = LEU_TEMPLATES['turn']['t/t']
IC_Leu_Turn_trans_gplus = LEU_TEMPLATES['turn']['t/g+']
IC_Leu_Turn_gplus_gminus = LEU_TEMPLATES['turn']['g+/g-']
IC_Leu_Turn_gplus_trans = LEU_TEMPLATES['turn']['g+/t']
IC_Leu_Turn_gplus_gplus = LEU_TEMPLATES['turn']['g+/g+']
IC_Leu_Bend = LEU_TEMPLATES['bend']['canonical']
IC_Leu_Bend_gminus_gminus = LEU_TEMPLATES['bend']['g-/g-']
IC_Leu_Bend_gminus_trans = LEU_TEMPLATES['bend']['g-/t']
IC_Leu_Bend_gminus_gplus = LEU_TEMPLATES['bend']['g-/g+']
IC_Leu_Bend_trans_gminus = LEU_TEMPLATES['bend']['t/g-']
IC_Leu_Bend_trans_trans = LEU_TEMPLATES['bend']['t/t']
IC_Leu_Bend_trans_gplus = LEU_TEMPLATES['bend']['t/g+']
IC_Leu_Bend_gplus_gminus = LEU_TEMPLATES['bend']['g+/g-']
IC_Leu_Bend_gplus_trans = LEU_TEMPLATES['bend']['g+/t']
IC_Leu_Bend_gplus_gplus = LEU_TEMPLATES['bend']['g+/g+']
IC_Leu_Coil = LEU_TEMPLATES['coil']['canonical']
IC_Leu_Coil_gminus_gminus = LEU_TEMPLATES['coil']['g-/g-']
IC_Leu_Coil_gminus_trans = LEU_TEMPLATES['coil']['g-/t']
IC_Leu_Coil_gminus_gplus = LEU_TEMPLATES['coil']['g-/g+']
IC_Leu_Coil_trans_gminus = LEU_TEMPLATES['coil']['t/g-']
IC_Leu_Coil_trans_trans = LEU_TEMPLATES['coil']['t/t']
IC_Leu_Coil_trans_gplus = LEU_TEMPLATES['coil']['t/g+']
IC_Leu_Coil_gplus_gminus = LEU_TEMPLATES['coil']['g+/g-']
IC_Leu_Coil_gplus_trans = LEU_TEMPLATES['coil']['g+/t']
IC_Leu_Coil_gplus_gplus = LEU_TEMPLATES['coil']['g+/g+']
IC_Leu_CisPeptide = LEU_TEMPLATES['cis-peptide-bond']['canonical']
IC_Leu_CisPeptide_gminus_gminus = LEU_TEMPLATES['cis-peptide-bond']['g-/g-']
IC_Leu_CisPeptide_gminus_trans = LEU_TEMPLATES['cis-peptide-bond']['g-/t']
IC_Leu_CisPeptide_gminus_gplus = LEU_TEMPLATES['cis-peptide-bond']['g-/g+']
IC_Leu_CisPeptide_trans_gminus = LEU_TEMPLATES['cis-peptide-bond']['t/g-']
IC_Leu_CisPeptide_trans_trans = LEU_TEMPLATES['cis-peptide-bond']['t/t']
IC_Leu_CisPeptide_trans_gplus = LEU_TEMPLATES['cis-peptide-bond']['t/g+']
IC_Leu_CisPeptide_gplus_gminus = LEU_TEMPLATES['cis-peptide-bond']['g+/g-']
IC_Leu_CisPeptide_gplus_trans = LEU_TEMPLATES['cis-peptide-bond']['g+/t']
IC_Leu_CisPeptide_gplus_gplus = LEU_TEMPLATES['cis-peptide-bond']['g+/g+']

__all__ = [
    'IC_Leu_Bend',
    'IC_Leu_Bend_gminus_gminus',
    'IC_Leu_Bend_gminus_gplus',
    'IC_Leu_Bend_gminus_trans',
    'IC_Leu_Bend_gplus_gminus',
    'IC_Leu_Bend_gplus_gplus',
    'IC_Leu_Bend_gplus_trans',
    'IC_Leu_Bend_trans_gminus',
    'IC_Leu_Bend_trans_gplus',
    'IC_Leu_Bend_trans_trans',
    'IC_Leu_Bridge',
    'IC_Leu_Bridge_gminus_gminus',
    'IC_Leu_Bridge_gminus_gplus',
    'IC_Leu_Bridge_gminus_trans',
    'IC_Leu_Bridge_gplus_gminus',
    'IC_Leu_Bridge_gplus_gplus',
    'IC_Leu_Bridge_gplus_trans',
    'IC_Leu_Bridge_trans_gminus',
    'IC_Leu_Bridge_trans_gplus',
    'IC_Leu_Bridge_trans_trans',
    'IC_Leu_CisPeptide',
    'IC_Leu_CisPeptide_gminus_gminus',
    'IC_Leu_CisPeptide_gminus_gplus',
    'IC_Leu_CisPeptide_gminus_trans',
    'IC_Leu_CisPeptide_gplus_gminus',
    'IC_Leu_CisPeptide_gplus_gplus',
    'IC_Leu_CisPeptide_gplus_trans',
    'IC_Leu_CisPeptide_trans_gminus',
    'IC_Leu_CisPeptide_trans_gplus',
    'IC_Leu_CisPeptide_trans_trans',
    'IC_Leu_Coil',
    'IC_Leu_Coil_gminus_gminus',
    'IC_Leu_Coil_gminus_gplus',
    'IC_Leu_Coil_gminus_trans',
    'IC_Leu_Coil_gplus_gminus',
    'IC_Leu_Coil_gplus_gplus',
    'IC_Leu_Coil_gplus_trans',
    'IC_Leu_Coil_trans_gminus',
    'IC_Leu_Coil_trans_gplus',
    'IC_Leu_Coil_trans_trans',
    'IC_Leu_Helix310',
    'IC_Leu_Helix310_gminus_gminus',
    'IC_Leu_Helix310_gminus_gplus',
    'IC_Leu_Helix310_gminus_trans',
    'IC_Leu_Helix310_gplus_gminus',
    'IC_Leu_Helix310_gplus_gplus',
    'IC_Leu_Helix310_gplus_trans',
    'IC_Leu_Helix310_trans_gminus',
    'IC_Leu_Helix310_trans_gplus',
    'IC_Leu_Helix310_trans_trans',
    'IC_Leu_HelixAlpha',
    'IC_Leu_HelixAlpha_gminus_gminus',
    'IC_Leu_HelixAlpha_gminus_gplus',
    'IC_Leu_HelixAlpha_gminus_trans',
    'IC_Leu_HelixAlpha_gplus_gminus',
    'IC_Leu_HelixAlpha_gplus_gplus',
    'IC_Leu_HelixAlpha_gplus_trans',
    'IC_Leu_HelixAlpha_trans_gminus',
    'IC_Leu_HelixAlpha_trans_gplus',
    'IC_Leu_HelixAlpha_trans_trans',
    'IC_Leu_HelixPPII',
    'IC_Leu_HelixPPII_gminus_gminus',
    'IC_Leu_HelixPPII_gminus_gplus',
    'IC_Leu_HelixPPII_gminus_trans',
    'IC_Leu_HelixPPII_gplus_gminus',
    'IC_Leu_HelixPPII_gplus_gplus',
    'IC_Leu_HelixPPII_gplus_trans',
    'IC_Leu_HelixPPII_trans_gminus',
    'IC_Leu_HelixPPII_trans_gplus',
    'IC_Leu_HelixPPII_trans_trans',
    'IC_Leu_HelixPi',
    'IC_Leu_HelixPi_gminus_gminus',
    'IC_Leu_HelixPi_gminus_gplus',
    'IC_Leu_HelixPi_gminus_trans',
    'IC_Leu_HelixPi_gplus_gminus',
    'IC_Leu_HelixPi_gplus_gplus',
    'IC_Leu_HelixPi_gplus_trans',
    'IC_Leu_HelixPi_trans_gminus',
    'IC_Leu_HelixPi_trans_gplus',
    'IC_Leu_HelixPi_trans_trans',
    'IC_Leu_Strand',
    'IC_Leu_StrandAntiParallel',
    'IC_Leu_StrandAntiParallel_gminus_gminus',
    'IC_Leu_StrandAntiParallel_gminus_gplus',
    'IC_Leu_StrandAntiParallel_gminus_trans',
    'IC_Leu_StrandAntiParallel_gplus_gminus',
    'IC_Leu_StrandAntiParallel_gplus_gplus',
    'IC_Leu_StrandAntiParallel_gplus_trans',
    'IC_Leu_StrandAntiParallel_trans_gminus',
    'IC_Leu_StrandAntiParallel_trans_gplus',
    'IC_Leu_StrandAntiParallel_trans_trans',
    'IC_Leu_StrandParallel',
    'IC_Leu_StrandParallel_gminus_gminus',
    'IC_Leu_StrandParallel_gminus_gplus',
    'IC_Leu_StrandParallel_gminus_trans',
    'IC_Leu_StrandParallel_gplus_gminus',
    'IC_Leu_StrandParallel_gplus_gplus',
    'IC_Leu_StrandParallel_gplus_trans',
    'IC_Leu_StrandParallel_trans_gminus',
    'IC_Leu_StrandParallel_trans_gplus',
    'IC_Leu_StrandParallel_trans_trans',
    'IC_Leu_Strand_gminus_gminus',
    'IC_Leu_Strand_gminus_gplus',
    'IC_Leu_Strand_gminus_trans',
    'IC_Leu_Strand_gplus_gminus',
    'IC_Leu_Strand_gplus_gplus',
    'IC_Leu_Strand_gplus_trans',
    'IC_Leu_Strand_trans_gminus',
    'IC_Leu_Strand_trans_gplus',
    'IC_Leu_Strand_trans_trans',
    'IC_Leu_Turn',
    'IC_Leu_Turn_gminus_gminus',
    'IC_Leu_Turn_gminus_gplus',
    'IC_Leu_Turn_gminus_trans',
    'IC_Leu_Turn_gplus_gminus',
    'IC_Leu_Turn_gplus_gplus',
    'IC_Leu_Turn_gplus_trans',
    'IC_Leu_Turn_trans_gminus',
    'IC_Leu_Turn_trans_gplus',
    'IC_Leu_Turn_trans_trans',
]
