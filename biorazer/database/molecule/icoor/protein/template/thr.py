# -*- coding: utf-8 -*-
"""Ideal {secondary-structure x rotamer} {names[aa]} {var} ``InternalCoord`` templates.

Each template is the ideal heavy-atom topology for {aa} anchored at
{{N, CA, C}}, with the side chain grown to ``rotamer`` and the
``phi``/``psi``/``omega`` of ``ss`` carried on the instance.  See
:mod:`biorazer.database.internal_coord_template` and
:mod:`._builder` for the build rule and the per-conformer caveat.
"""

from biorazer.database.molecule.icoor.protein.template._builder import make_residue_templates

RESN = "THR"
#: {ss: {rotamer: :class:`~biorazer.structure.objects.InternalCoord`}}
THR_TEMPLATES = make_residue_templates(RESN)

# === GENERATED IC_* named templates - do not edit below ===
# Regenerate with: /opt/envs/BioRazer/bin/python scripts/generate_internal_coord_template_named.py

IC_Thr_HelixAlpha = THR_TEMPLATES['alpha-helix']['canonical']
IC_Thr_HelixAlpha_gminus = THR_TEMPLATES['alpha-helix']['g-']
IC_Thr_HelixAlpha_trans = THR_TEMPLATES['alpha-helix']['t']
IC_Thr_HelixAlpha_gplus = THR_TEMPLATES['alpha-helix']['g+']
IC_Thr_Helix310 = THR_TEMPLATES['3-10-helix']['canonical']
IC_Thr_Helix310_gminus = THR_TEMPLATES['3-10-helix']['g-']
IC_Thr_Helix310_trans = THR_TEMPLATES['3-10-helix']['t']
IC_Thr_Helix310_gplus = THR_TEMPLATES['3-10-helix']['g+']
IC_Thr_HelixPi = THR_TEMPLATES['pi-helix']['canonical']
IC_Thr_HelixPi_gminus = THR_TEMPLATES['pi-helix']['g-']
IC_Thr_HelixPi_trans = THR_TEMPLATES['pi-helix']['t']
IC_Thr_HelixPi_gplus = THR_TEMPLATES['pi-helix']['g+']
IC_Thr_HelixPPII = THR_TEMPLATES['polyproline-II']['canonical']
IC_Thr_HelixPPII_gminus = THR_TEMPLATES['polyproline-II']['g-']
IC_Thr_HelixPPII_trans = THR_TEMPLATES['polyproline-II']['t']
IC_Thr_HelixPPII_gplus = THR_TEMPLATES['polyproline-II']['g+']
IC_Thr_Strand = THR_TEMPLATES['beta-strand']['canonical']
IC_Thr_Strand_gminus = THR_TEMPLATES['beta-strand']['g-']
IC_Thr_Strand_trans = THR_TEMPLATES['beta-strand']['t']
IC_Thr_Strand_gplus = THR_TEMPLATES['beta-strand']['g+']
IC_Thr_StrandParallel = THR_TEMPLATES['parallel-beta-strand']['canonical']
IC_Thr_StrandParallel_gminus = THR_TEMPLATES['parallel-beta-strand']['g-']
IC_Thr_StrandParallel_trans = THR_TEMPLATES['parallel-beta-strand']['t']
IC_Thr_StrandParallel_gplus = THR_TEMPLATES['parallel-beta-strand']['g+']
IC_Thr_StrandAntiParallel = THR_TEMPLATES['antiparallel-beta-strand']['canonical']
IC_Thr_StrandAntiParallel_gminus = THR_TEMPLATES['antiparallel-beta-strand']['g-']
IC_Thr_StrandAntiParallel_trans = THR_TEMPLATES['antiparallel-beta-strand']['t']
IC_Thr_StrandAntiParallel_gplus = THR_TEMPLATES['antiparallel-beta-strand']['g+']
IC_Thr_Bridge = THR_TEMPLATES['beta-bridge']['canonical']
IC_Thr_Bridge_gminus = THR_TEMPLATES['beta-bridge']['g-']
IC_Thr_Bridge_trans = THR_TEMPLATES['beta-bridge']['t']
IC_Thr_Bridge_gplus = THR_TEMPLATES['beta-bridge']['g+']
IC_Thr_Turn = THR_TEMPLATES['turn']['canonical']
IC_Thr_Turn_gminus = THR_TEMPLATES['turn']['g-']
IC_Thr_Turn_trans = THR_TEMPLATES['turn']['t']
IC_Thr_Turn_gplus = THR_TEMPLATES['turn']['g+']
IC_Thr_Bend = THR_TEMPLATES['bend']['canonical']
IC_Thr_Bend_gminus = THR_TEMPLATES['bend']['g-']
IC_Thr_Bend_trans = THR_TEMPLATES['bend']['t']
IC_Thr_Bend_gplus = THR_TEMPLATES['bend']['g+']
IC_Thr_Coil = THR_TEMPLATES['coil']['canonical']
IC_Thr_Coil_gminus = THR_TEMPLATES['coil']['g-']
IC_Thr_Coil_trans = THR_TEMPLATES['coil']['t']
IC_Thr_Coil_gplus = THR_TEMPLATES['coil']['g+']
IC_Thr_CisPeptide = THR_TEMPLATES['cis-peptide-bond']['canonical']
IC_Thr_CisPeptide_gminus = THR_TEMPLATES['cis-peptide-bond']['g-']
IC_Thr_CisPeptide_trans = THR_TEMPLATES['cis-peptide-bond']['t']
IC_Thr_CisPeptide_gplus = THR_TEMPLATES['cis-peptide-bond']['g+']

__all__ = [
    'IC_Thr_Bend',
    'IC_Thr_Bend_gminus',
    'IC_Thr_Bend_gplus',
    'IC_Thr_Bend_trans',
    'IC_Thr_Bridge',
    'IC_Thr_Bridge_gminus',
    'IC_Thr_Bridge_gplus',
    'IC_Thr_Bridge_trans',
    'IC_Thr_CisPeptide',
    'IC_Thr_CisPeptide_gminus',
    'IC_Thr_CisPeptide_gplus',
    'IC_Thr_CisPeptide_trans',
    'IC_Thr_Coil',
    'IC_Thr_Coil_gminus',
    'IC_Thr_Coil_gplus',
    'IC_Thr_Coil_trans',
    'IC_Thr_Helix310',
    'IC_Thr_Helix310_gminus',
    'IC_Thr_Helix310_gplus',
    'IC_Thr_Helix310_trans',
    'IC_Thr_HelixAlpha',
    'IC_Thr_HelixAlpha_gminus',
    'IC_Thr_HelixAlpha_gplus',
    'IC_Thr_HelixAlpha_trans',
    'IC_Thr_HelixPPII',
    'IC_Thr_HelixPPII_gminus',
    'IC_Thr_HelixPPII_gplus',
    'IC_Thr_HelixPPII_trans',
    'IC_Thr_HelixPi',
    'IC_Thr_HelixPi_gminus',
    'IC_Thr_HelixPi_gplus',
    'IC_Thr_HelixPi_trans',
    'IC_Thr_Strand',
    'IC_Thr_StrandAntiParallel',
    'IC_Thr_StrandAntiParallel_gminus',
    'IC_Thr_StrandAntiParallel_gplus',
    'IC_Thr_StrandAntiParallel_trans',
    'IC_Thr_StrandParallel',
    'IC_Thr_StrandParallel_gminus',
    'IC_Thr_StrandParallel_gplus',
    'IC_Thr_StrandParallel_trans',
    'IC_Thr_Strand_gminus',
    'IC_Thr_Strand_gplus',
    'IC_Thr_Strand_trans',
    'IC_Thr_Turn',
    'IC_Thr_Turn_gminus',
    'IC_Thr_Turn_gplus',
    'IC_Thr_Turn_trans',
]
