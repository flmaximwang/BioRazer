# -*- coding: utf-8 -*-
"""Ideal {secondary-structure x rotamer} {names[aa]} {var} ``InternalCoord`` templates.

Each template is the ideal heavy-atom topology for {aa} anchored at
{{N, CA, C}}, with the side chain grown to ``rotamer`` and the
``phi``/``psi``/``omega`` of ``ss`` carried on the instance.  See
:mod:`biorazer.database.molecule.icoor.protein.template` and
:mod:`._builder` for the build rule and the per-conformer caveat.
"""

from biorazer.database.molecule.icoor.protein.template._builder import make_residue_templates

RESN = "SER"
#: {ss: {rotamer: :class:`~biorazer.structure.objects.InternalCoord`}}
SER_TEMPLATES = make_residue_templates(RESN)

# === GENERATED IC_* named templates - do not edit below ===
# Regenerate with: /opt/envs/BioRazer/bin/python scripts/generate_internal_coord_template_named.py

IC_Ser_HelixAlpha = SER_TEMPLATES['alpha-helix']['canonical']
IC_Ser_HelixAlpha_gminus = SER_TEMPLATES['alpha-helix']['g-']
IC_Ser_HelixAlpha_trans = SER_TEMPLATES['alpha-helix']['t']
IC_Ser_HelixAlpha_gplus = SER_TEMPLATES['alpha-helix']['g+']
IC_Ser_Helix310 = SER_TEMPLATES['3-10-helix']['canonical']
IC_Ser_Helix310_gminus = SER_TEMPLATES['3-10-helix']['g-']
IC_Ser_Helix310_trans = SER_TEMPLATES['3-10-helix']['t']
IC_Ser_Helix310_gplus = SER_TEMPLATES['3-10-helix']['g+']
IC_Ser_HelixPi = SER_TEMPLATES['pi-helix']['canonical']
IC_Ser_HelixPi_gminus = SER_TEMPLATES['pi-helix']['g-']
IC_Ser_HelixPi_trans = SER_TEMPLATES['pi-helix']['t']
IC_Ser_HelixPi_gplus = SER_TEMPLATES['pi-helix']['g+']
IC_Ser_HelixPPII = SER_TEMPLATES['polyproline-II']['canonical']
IC_Ser_HelixPPII_gminus = SER_TEMPLATES['polyproline-II']['g-']
IC_Ser_HelixPPII_trans = SER_TEMPLATES['polyproline-II']['t']
IC_Ser_HelixPPII_gplus = SER_TEMPLATES['polyproline-II']['g+']
IC_Ser_Strand = SER_TEMPLATES['beta-strand']['canonical']
IC_Ser_Strand_gminus = SER_TEMPLATES['beta-strand']['g-']
IC_Ser_Strand_trans = SER_TEMPLATES['beta-strand']['t']
IC_Ser_Strand_gplus = SER_TEMPLATES['beta-strand']['g+']
IC_Ser_StrandParallel = SER_TEMPLATES['parallel-beta-strand']['canonical']
IC_Ser_StrandParallel_gminus = SER_TEMPLATES['parallel-beta-strand']['g-']
IC_Ser_StrandParallel_trans = SER_TEMPLATES['parallel-beta-strand']['t']
IC_Ser_StrandParallel_gplus = SER_TEMPLATES['parallel-beta-strand']['g+']
IC_Ser_StrandAntiParallel = SER_TEMPLATES['antiparallel-beta-strand']['canonical']
IC_Ser_StrandAntiParallel_gminus = SER_TEMPLATES['antiparallel-beta-strand']['g-']
IC_Ser_StrandAntiParallel_trans = SER_TEMPLATES['antiparallel-beta-strand']['t']
IC_Ser_StrandAntiParallel_gplus = SER_TEMPLATES['antiparallel-beta-strand']['g+']
IC_Ser_Bridge = SER_TEMPLATES['beta-bridge']['canonical']
IC_Ser_Bridge_gminus = SER_TEMPLATES['beta-bridge']['g-']
IC_Ser_Bridge_trans = SER_TEMPLATES['beta-bridge']['t']
IC_Ser_Bridge_gplus = SER_TEMPLATES['beta-bridge']['g+']
IC_Ser_Turn = SER_TEMPLATES['turn']['canonical']
IC_Ser_Turn_gminus = SER_TEMPLATES['turn']['g-']
IC_Ser_Turn_trans = SER_TEMPLATES['turn']['t']
IC_Ser_Turn_gplus = SER_TEMPLATES['turn']['g+']
IC_Ser_Bend = SER_TEMPLATES['bend']['canonical']
IC_Ser_Bend_gminus = SER_TEMPLATES['bend']['g-']
IC_Ser_Bend_trans = SER_TEMPLATES['bend']['t']
IC_Ser_Bend_gplus = SER_TEMPLATES['bend']['g+']
IC_Ser_Coil = SER_TEMPLATES['coil']['canonical']
IC_Ser_Coil_gminus = SER_TEMPLATES['coil']['g-']
IC_Ser_Coil_trans = SER_TEMPLATES['coil']['t']
IC_Ser_Coil_gplus = SER_TEMPLATES['coil']['g+']
IC_Ser_CisPeptide = SER_TEMPLATES['cis-peptide-bond']['canonical']
IC_Ser_CisPeptide_gminus = SER_TEMPLATES['cis-peptide-bond']['g-']
IC_Ser_CisPeptide_trans = SER_TEMPLATES['cis-peptide-bond']['t']
IC_Ser_CisPeptide_gplus = SER_TEMPLATES['cis-peptide-bond']['g+']

__all__ = [
    'IC_Ser_Bend',
    'IC_Ser_Bend_gminus',
    'IC_Ser_Bend_gplus',
    'IC_Ser_Bend_trans',
    'IC_Ser_Bridge',
    'IC_Ser_Bridge_gminus',
    'IC_Ser_Bridge_gplus',
    'IC_Ser_Bridge_trans',
    'IC_Ser_CisPeptide',
    'IC_Ser_CisPeptide_gminus',
    'IC_Ser_CisPeptide_gplus',
    'IC_Ser_CisPeptide_trans',
    'IC_Ser_Coil',
    'IC_Ser_Coil_gminus',
    'IC_Ser_Coil_gplus',
    'IC_Ser_Coil_trans',
    'IC_Ser_Helix310',
    'IC_Ser_Helix310_gminus',
    'IC_Ser_Helix310_gplus',
    'IC_Ser_Helix310_trans',
    'IC_Ser_HelixAlpha',
    'IC_Ser_HelixAlpha_gminus',
    'IC_Ser_HelixAlpha_gplus',
    'IC_Ser_HelixAlpha_trans',
    'IC_Ser_HelixPPII',
    'IC_Ser_HelixPPII_gminus',
    'IC_Ser_HelixPPII_gplus',
    'IC_Ser_HelixPPII_trans',
    'IC_Ser_HelixPi',
    'IC_Ser_HelixPi_gminus',
    'IC_Ser_HelixPi_gplus',
    'IC_Ser_HelixPi_trans',
    'IC_Ser_Strand',
    'IC_Ser_StrandAntiParallel',
    'IC_Ser_StrandAntiParallel_gminus',
    'IC_Ser_StrandAntiParallel_gplus',
    'IC_Ser_StrandAntiParallel_trans',
    'IC_Ser_StrandParallel',
    'IC_Ser_StrandParallel_gminus',
    'IC_Ser_StrandParallel_gplus',
    'IC_Ser_StrandParallel_trans',
    'IC_Ser_Strand_gminus',
    'IC_Ser_Strand_gplus',
    'IC_Ser_Strand_trans',
    'IC_Ser_Turn',
    'IC_Ser_Turn_gminus',
    'IC_Ser_Turn_gplus',
    'IC_Ser_Turn_trans',
]
