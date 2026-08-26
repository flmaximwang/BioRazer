# -*- coding: utf-8 -*-
"""Ideal {secondary-structure x rotamer} {names[aa]} {var} ``InternalCoord`` templates.

Each template is the ideal heavy-atom topology for {aa} anchored at
{{N, CA, C}}, with the side chain grown to ``rotamer`` and the
``phi``/``psi``/``omega`` of ``ss`` carried on the instance.  See
:mod:`biorazer.database.internal_coord_template` and
:mod:`._builder` for the build rule and the per-conformer caveat.
"""

from biorazer.database.internal_coord_template._builder import make_residue_templates

RESN = "ASN"
#: {ss: {rotamer: :class:`~biorazer.structure.objects.InternalCoord`}}
ASN_TEMPLATES = make_residue_templates(RESN)

# === GENERATED IC_* named templates - do not edit below ===
# Regenerate with: /opt/envs/BioRazer/bin/python scripts/generate_internal_coord_template_named.py

IC_Asn_HelixAlpha = ASN_TEMPLATES['alpha-helix']['canonical']
IC_Asn_HelixAlpha_gminus_gminus = ASN_TEMPLATES['alpha-helix']['g-/g-']
IC_Asn_HelixAlpha_gminus_trans = ASN_TEMPLATES['alpha-helix']['g-/t']
IC_Asn_HelixAlpha_gminus_gplus = ASN_TEMPLATES['alpha-helix']['g-/g+']
IC_Asn_HelixAlpha_trans_gminus = ASN_TEMPLATES['alpha-helix']['t/g-']
IC_Asn_HelixAlpha_trans_trans = ASN_TEMPLATES['alpha-helix']['t/t']
IC_Asn_HelixAlpha_trans_gplus = ASN_TEMPLATES['alpha-helix']['t/g+']
IC_Asn_HelixAlpha_gplus_gminus = ASN_TEMPLATES['alpha-helix']['g+/g-']
IC_Asn_HelixAlpha_gplus_trans = ASN_TEMPLATES['alpha-helix']['g+/t']
IC_Asn_HelixAlpha_gplus_gplus = ASN_TEMPLATES['alpha-helix']['g+/g+']
IC_Asn_Helix310 = ASN_TEMPLATES['3-10-helix']['canonical']
IC_Asn_Helix310_gminus_gminus = ASN_TEMPLATES['3-10-helix']['g-/g-']
IC_Asn_Helix310_gminus_trans = ASN_TEMPLATES['3-10-helix']['g-/t']
IC_Asn_Helix310_gminus_gplus = ASN_TEMPLATES['3-10-helix']['g-/g+']
IC_Asn_Helix310_trans_gminus = ASN_TEMPLATES['3-10-helix']['t/g-']
IC_Asn_Helix310_trans_trans = ASN_TEMPLATES['3-10-helix']['t/t']
IC_Asn_Helix310_trans_gplus = ASN_TEMPLATES['3-10-helix']['t/g+']
IC_Asn_Helix310_gplus_gminus = ASN_TEMPLATES['3-10-helix']['g+/g-']
IC_Asn_Helix310_gplus_trans = ASN_TEMPLATES['3-10-helix']['g+/t']
IC_Asn_Helix310_gplus_gplus = ASN_TEMPLATES['3-10-helix']['g+/g+']
IC_Asn_HelixPi = ASN_TEMPLATES['pi-helix']['canonical']
IC_Asn_HelixPi_gminus_gminus = ASN_TEMPLATES['pi-helix']['g-/g-']
IC_Asn_HelixPi_gminus_trans = ASN_TEMPLATES['pi-helix']['g-/t']
IC_Asn_HelixPi_gminus_gplus = ASN_TEMPLATES['pi-helix']['g-/g+']
IC_Asn_HelixPi_trans_gminus = ASN_TEMPLATES['pi-helix']['t/g-']
IC_Asn_HelixPi_trans_trans = ASN_TEMPLATES['pi-helix']['t/t']
IC_Asn_HelixPi_trans_gplus = ASN_TEMPLATES['pi-helix']['t/g+']
IC_Asn_HelixPi_gplus_gminus = ASN_TEMPLATES['pi-helix']['g+/g-']
IC_Asn_HelixPi_gplus_trans = ASN_TEMPLATES['pi-helix']['g+/t']
IC_Asn_HelixPi_gplus_gplus = ASN_TEMPLATES['pi-helix']['g+/g+']
IC_Asn_HelixPPII = ASN_TEMPLATES['polyproline-II']['canonical']
IC_Asn_HelixPPII_gminus_gminus = ASN_TEMPLATES['polyproline-II']['g-/g-']
IC_Asn_HelixPPII_gminus_trans = ASN_TEMPLATES['polyproline-II']['g-/t']
IC_Asn_HelixPPII_gminus_gplus = ASN_TEMPLATES['polyproline-II']['g-/g+']
IC_Asn_HelixPPII_trans_gminus = ASN_TEMPLATES['polyproline-II']['t/g-']
IC_Asn_HelixPPII_trans_trans = ASN_TEMPLATES['polyproline-II']['t/t']
IC_Asn_HelixPPII_trans_gplus = ASN_TEMPLATES['polyproline-II']['t/g+']
IC_Asn_HelixPPII_gplus_gminus = ASN_TEMPLATES['polyproline-II']['g+/g-']
IC_Asn_HelixPPII_gplus_trans = ASN_TEMPLATES['polyproline-II']['g+/t']
IC_Asn_HelixPPII_gplus_gplus = ASN_TEMPLATES['polyproline-II']['g+/g+']
IC_Asn_Strand = ASN_TEMPLATES['beta-strand']['canonical']
IC_Asn_Strand_gminus_gminus = ASN_TEMPLATES['beta-strand']['g-/g-']
IC_Asn_Strand_gminus_trans = ASN_TEMPLATES['beta-strand']['g-/t']
IC_Asn_Strand_gminus_gplus = ASN_TEMPLATES['beta-strand']['g-/g+']
IC_Asn_Strand_trans_gminus = ASN_TEMPLATES['beta-strand']['t/g-']
IC_Asn_Strand_trans_trans = ASN_TEMPLATES['beta-strand']['t/t']
IC_Asn_Strand_trans_gplus = ASN_TEMPLATES['beta-strand']['t/g+']
IC_Asn_Strand_gplus_gminus = ASN_TEMPLATES['beta-strand']['g+/g-']
IC_Asn_Strand_gplus_trans = ASN_TEMPLATES['beta-strand']['g+/t']
IC_Asn_Strand_gplus_gplus = ASN_TEMPLATES['beta-strand']['g+/g+']
IC_Asn_StrandParallel = ASN_TEMPLATES['parallel-beta-strand']['canonical']
IC_Asn_StrandParallel_gminus_gminus = ASN_TEMPLATES['parallel-beta-strand']['g-/g-']
IC_Asn_StrandParallel_gminus_trans = ASN_TEMPLATES['parallel-beta-strand']['g-/t']
IC_Asn_StrandParallel_gminus_gplus = ASN_TEMPLATES['parallel-beta-strand']['g-/g+']
IC_Asn_StrandParallel_trans_gminus = ASN_TEMPLATES['parallel-beta-strand']['t/g-']
IC_Asn_StrandParallel_trans_trans = ASN_TEMPLATES['parallel-beta-strand']['t/t']
IC_Asn_StrandParallel_trans_gplus = ASN_TEMPLATES['parallel-beta-strand']['t/g+']
IC_Asn_StrandParallel_gplus_gminus = ASN_TEMPLATES['parallel-beta-strand']['g+/g-']
IC_Asn_StrandParallel_gplus_trans = ASN_TEMPLATES['parallel-beta-strand']['g+/t']
IC_Asn_StrandParallel_gplus_gplus = ASN_TEMPLATES['parallel-beta-strand']['g+/g+']
IC_Asn_StrandAntiParallel = ASN_TEMPLATES['antiparallel-beta-strand']['canonical']
IC_Asn_StrandAntiParallel_gminus_gminus = ASN_TEMPLATES['antiparallel-beta-strand']['g-/g-']
IC_Asn_StrandAntiParallel_gminus_trans = ASN_TEMPLATES['antiparallel-beta-strand']['g-/t']
IC_Asn_StrandAntiParallel_gminus_gplus = ASN_TEMPLATES['antiparallel-beta-strand']['g-/g+']
IC_Asn_StrandAntiParallel_trans_gminus = ASN_TEMPLATES['antiparallel-beta-strand']['t/g-']
IC_Asn_StrandAntiParallel_trans_trans = ASN_TEMPLATES['antiparallel-beta-strand']['t/t']
IC_Asn_StrandAntiParallel_trans_gplus = ASN_TEMPLATES['antiparallel-beta-strand']['t/g+']
IC_Asn_StrandAntiParallel_gplus_gminus = ASN_TEMPLATES['antiparallel-beta-strand']['g+/g-']
IC_Asn_StrandAntiParallel_gplus_trans = ASN_TEMPLATES['antiparallel-beta-strand']['g+/t']
IC_Asn_StrandAntiParallel_gplus_gplus = ASN_TEMPLATES['antiparallel-beta-strand']['g+/g+']
IC_Asn_Bridge = ASN_TEMPLATES['beta-bridge']['canonical']
IC_Asn_Bridge_gminus_gminus = ASN_TEMPLATES['beta-bridge']['g-/g-']
IC_Asn_Bridge_gminus_trans = ASN_TEMPLATES['beta-bridge']['g-/t']
IC_Asn_Bridge_gminus_gplus = ASN_TEMPLATES['beta-bridge']['g-/g+']
IC_Asn_Bridge_trans_gminus = ASN_TEMPLATES['beta-bridge']['t/g-']
IC_Asn_Bridge_trans_trans = ASN_TEMPLATES['beta-bridge']['t/t']
IC_Asn_Bridge_trans_gplus = ASN_TEMPLATES['beta-bridge']['t/g+']
IC_Asn_Bridge_gplus_gminus = ASN_TEMPLATES['beta-bridge']['g+/g-']
IC_Asn_Bridge_gplus_trans = ASN_TEMPLATES['beta-bridge']['g+/t']
IC_Asn_Bridge_gplus_gplus = ASN_TEMPLATES['beta-bridge']['g+/g+']
IC_Asn_Turn = ASN_TEMPLATES['turn']['canonical']
IC_Asn_Turn_gminus_gminus = ASN_TEMPLATES['turn']['g-/g-']
IC_Asn_Turn_gminus_trans = ASN_TEMPLATES['turn']['g-/t']
IC_Asn_Turn_gminus_gplus = ASN_TEMPLATES['turn']['g-/g+']
IC_Asn_Turn_trans_gminus = ASN_TEMPLATES['turn']['t/g-']
IC_Asn_Turn_trans_trans = ASN_TEMPLATES['turn']['t/t']
IC_Asn_Turn_trans_gplus = ASN_TEMPLATES['turn']['t/g+']
IC_Asn_Turn_gplus_gminus = ASN_TEMPLATES['turn']['g+/g-']
IC_Asn_Turn_gplus_trans = ASN_TEMPLATES['turn']['g+/t']
IC_Asn_Turn_gplus_gplus = ASN_TEMPLATES['turn']['g+/g+']
IC_Asn_Bend = ASN_TEMPLATES['bend']['canonical']
IC_Asn_Bend_gminus_gminus = ASN_TEMPLATES['bend']['g-/g-']
IC_Asn_Bend_gminus_trans = ASN_TEMPLATES['bend']['g-/t']
IC_Asn_Bend_gminus_gplus = ASN_TEMPLATES['bend']['g-/g+']
IC_Asn_Bend_trans_gminus = ASN_TEMPLATES['bend']['t/g-']
IC_Asn_Bend_trans_trans = ASN_TEMPLATES['bend']['t/t']
IC_Asn_Bend_trans_gplus = ASN_TEMPLATES['bend']['t/g+']
IC_Asn_Bend_gplus_gminus = ASN_TEMPLATES['bend']['g+/g-']
IC_Asn_Bend_gplus_trans = ASN_TEMPLATES['bend']['g+/t']
IC_Asn_Bend_gplus_gplus = ASN_TEMPLATES['bend']['g+/g+']
IC_Asn_Coil = ASN_TEMPLATES['coil']['canonical']
IC_Asn_Coil_gminus_gminus = ASN_TEMPLATES['coil']['g-/g-']
IC_Asn_Coil_gminus_trans = ASN_TEMPLATES['coil']['g-/t']
IC_Asn_Coil_gminus_gplus = ASN_TEMPLATES['coil']['g-/g+']
IC_Asn_Coil_trans_gminus = ASN_TEMPLATES['coil']['t/g-']
IC_Asn_Coil_trans_trans = ASN_TEMPLATES['coil']['t/t']
IC_Asn_Coil_trans_gplus = ASN_TEMPLATES['coil']['t/g+']
IC_Asn_Coil_gplus_gminus = ASN_TEMPLATES['coil']['g+/g-']
IC_Asn_Coil_gplus_trans = ASN_TEMPLATES['coil']['g+/t']
IC_Asn_Coil_gplus_gplus = ASN_TEMPLATES['coil']['g+/g+']
IC_Asn_CisPeptide = ASN_TEMPLATES['cis-peptide-bond']['canonical']
IC_Asn_CisPeptide_gminus_gminus = ASN_TEMPLATES['cis-peptide-bond']['g-/g-']
IC_Asn_CisPeptide_gminus_trans = ASN_TEMPLATES['cis-peptide-bond']['g-/t']
IC_Asn_CisPeptide_gminus_gplus = ASN_TEMPLATES['cis-peptide-bond']['g-/g+']
IC_Asn_CisPeptide_trans_gminus = ASN_TEMPLATES['cis-peptide-bond']['t/g-']
IC_Asn_CisPeptide_trans_trans = ASN_TEMPLATES['cis-peptide-bond']['t/t']
IC_Asn_CisPeptide_trans_gplus = ASN_TEMPLATES['cis-peptide-bond']['t/g+']
IC_Asn_CisPeptide_gplus_gminus = ASN_TEMPLATES['cis-peptide-bond']['g+/g-']
IC_Asn_CisPeptide_gplus_trans = ASN_TEMPLATES['cis-peptide-bond']['g+/t']
IC_Asn_CisPeptide_gplus_gplus = ASN_TEMPLATES['cis-peptide-bond']['g+/g+']

__all__ = [
    'IC_Asn_Bend',
    'IC_Asn_Bend_gminus_gminus',
    'IC_Asn_Bend_gminus_gplus',
    'IC_Asn_Bend_gminus_trans',
    'IC_Asn_Bend_gplus_gminus',
    'IC_Asn_Bend_gplus_gplus',
    'IC_Asn_Bend_gplus_trans',
    'IC_Asn_Bend_trans_gminus',
    'IC_Asn_Bend_trans_gplus',
    'IC_Asn_Bend_trans_trans',
    'IC_Asn_Bridge',
    'IC_Asn_Bridge_gminus_gminus',
    'IC_Asn_Bridge_gminus_gplus',
    'IC_Asn_Bridge_gminus_trans',
    'IC_Asn_Bridge_gplus_gminus',
    'IC_Asn_Bridge_gplus_gplus',
    'IC_Asn_Bridge_gplus_trans',
    'IC_Asn_Bridge_trans_gminus',
    'IC_Asn_Bridge_trans_gplus',
    'IC_Asn_Bridge_trans_trans',
    'IC_Asn_CisPeptide',
    'IC_Asn_CisPeptide_gminus_gminus',
    'IC_Asn_CisPeptide_gminus_gplus',
    'IC_Asn_CisPeptide_gminus_trans',
    'IC_Asn_CisPeptide_gplus_gminus',
    'IC_Asn_CisPeptide_gplus_gplus',
    'IC_Asn_CisPeptide_gplus_trans',
    'IC_Asn_CisPeptide_trans_gminus',
    'IC_Asn_CisPeptide_trans_gplus',
    'IC_Asn_CisPeptide_trans_trans',
    'IC_Asn_Coil',
    'IC_Asn_Coil_gminus_gminus',
    'IC_Asn_Coil_gminus_gplus',
    'IC_Asn_Coil_gminus_trans',
    'IC_Asn_Coil_gplus_gminus',
    'IC_Asn_Coil_gplus_gplus',
    'IC_Asn_Coil_gplus_trans',
    'IC_Asn_Coil_trans_gminus',
    'IC_Asn_Coil_trans_gplus',
    'IC_Asn_Coil_trans_trans',
    'IC_Asn_Helix310',
    'IC_Asn_Helix310_gminus_gminus',
    'IC_Asn_Helix310_gminus_gplus',
    'IC_Asn_Helix310_gminus_trans',
    'IC_Asn_Helix310_gplus_gminus',
    'IC_Asn_Helix310_gplus_gplus',
    'IC_Asn_Helix310_gplus_trans',
    'IC_Asn_Helix310_trans_gminus',
    'IC_Asn_Helix310_trans_gplus',
    'IC_Asn_Helix310_trans_trans',
    'IC_Asn_HelixAlpha',
    'IC_Asn_HelixAlpha_gminus_gminus',
    'IC_Asn_HelixAlpha_gminus_gplus',
    'IC_Asn_HelixAlpha_gminus_trans',
    'IC_Asn_HelixAlpha_gplus_gminus',
    'IC_Asn_HelixAlpha_gplus_gplus',
    'IC_Asn_HelixAlpha_gplus_trans',
    'IC_Asn_HelixAlpha_trans_gminus',
    'IC_Asn_HelixAlpha_trans_gplus',
    'IC_Asn_HelixAlpha_trans_trans',
    'IC_Asn_HelixPPII',
    'IC_Asn_HelixPPII_gminus_gminus',
    'IC_Asn_HelixPPII_gminus_gplus',
    'IC_Asn_HelixPPII_gminus_trans',
    'IC_Asn_HelixPPII_gplus_gminus',
    'IC_Asn_HelixPPII_gplus_gplus',
    'IC_Asn_HelixPPII_gplus_trans',
    'IC_Asn_HelixPPII_trans_gminus',
    'IC_Asn_HelixPPII_trans_gplus',
    'IC_Asn_HelixPPII_trans_trans',
    'IC_Asn_HelixPi',
    'IC_Asn_HelixPi_gminus_gminus',
    'IC_Asn_HelixPi_gminus_gplus',
    'IC_Asn_HelixPi_gminus_trans',
    'IC_Asn_HelixPi_gplus_gminus',
    'IC_Asn_HelixPi_gplus_gplus',
    'IC_Asn_HelixPi_gplus_trans',
    'IC_Asn_HelixPi_trans_gminus',
    'IC_Asn_HelixPi_trans_gplus',
    'IC_Asn_HelixPi_trans_trans',
    'IC_Asn_Strand',
    'IC_Asn_StrandAntiParallel',
    'IC_Asn_StrandAntiParallel_gminus_gminus',
    'IC_Asn_StrandAntiParallel_gminus_gplus',
    'IC_Asn_StrandAntiParallel_gminus_trans',
    'IC_Asn_StrandAntiParallel_gplus_gminus',
    'IC_Asn_StrandAntiParallel_gplus_gplus',
    'IC_Asn_StrandAntiParallel_gplus_trans',
    'IC_Asn_StrandAntiParallel_trans_gminus',
    'IC_Asn_StrandAntiParallel_trans_gplus',
    'IC_Asn_StrandAntiParallel_trans_trans',
    'IC_Asn_StrandParallel',
    'IC_Asn_StrandParallel_gminus_gminus',
    'IC_Asn_StrandParallel_gminus_gplus',
    'IC_Asn_StrandParallel_gminus_trans',
    'IC_Asn_StrandParallel_gplus_gminus',
    'IC_Asn_StrandParallel_gplus_gplus',
    'IC_Asn_StrandParallel_gplus_trans',
    'IC_Asn_StrandParallel_trans_gminus',
    'IC_Asn_StrandParallel_trans_gplus',
    'IC_Asn_StrandParallel_trans_trans',
    'IC_Asn_Strand_gminus_gminus',
    'IC_Asn_Strand_gminus_gplus',
    'IC_Asn_Strand_gminus_trans',
    'IC_Asn_Strand_gplus_gminus',
    'IC_Asn_Strand_gplus_gplus',
    'IC_Asn_Strand_gplus_trans',
    'IC_Asn_Strand_trans_gminus',
    'IC_Asn_Strand_trans_gplus',
    'IC_Asn_Strand_trans_trans',
    'IC_Asn_Turn',
    'IC_Asn_Turn_gminus_gminus',
    'IC_Asn_Turn_gminus_gplus',
    'IC_Asn_Turn_gminus_trans',
    'IC_Asn_Turn_gplus_gminus',
    'IC_Asn_Turn_gplus_gplus',
    'IC_Asn_Turn_gplus_trans',
    'IC_Asn_Turn_trans_gminus',
    'IC_Asn_Turn_trans_gplus',
    'IC_Asn_Turn_trans_trans',
]
