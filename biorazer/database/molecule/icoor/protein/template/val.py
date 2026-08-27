# -*- coding: utf-8 -*-
"""Ideal {secondary-structure x rotamer} {names[aa]} {var} ``InternalCoord`` templates.

Each template is the ideal heavy-atom topology for {aa} anchored at
{{N, CA, C}}, with the side chain grown to ``rotamer`` and the
``phi``/``psi``/``omega`` of ``ss`` carried on the instance.  See
:mod:`biorazer.database.molecule.icoor.protein.template` and
:mod:`._builder` for the build rule and the per-conformer caveat.
"""

from biorazer.database.molecule.icoor.protein.template._builder import build_template

# === GENERATED IC_* named templates - do not edit below ===
# Regenerate with: /opt/envs/BioRazer/bin/python scripts/generate_internal_coord_template_named.py

IC_Val_HelixAlpha = build_template('VAL', 'alpha-helix', 'canonical')
IC_Val_HelixAlpha_gminus = build_template('VAL', 'alpha-helix', 'g-')
IC_Val_HelixAlpha_trans = build_template('VAL', 'alpha-helix', 't')
IC_Val_HelixAlpha_gplus = build_template('VAL', 'alpha-helix', 'g+')
IC_Val_Helix310 = build_template('VAL', '3-10-helix', 'canonical')
IC_Val_Helix310_gminus = build_template('VAL', '3-10-helix', 'g-')
IC_Val_Helix310_trans = build_template('VAL', '3-10-helix', 't')
IC_Val_Helix310_gplus = build_template('VAL', '3-10-helix', 'g+')
IC_Val_HelixPi = build_template('VAL', 'pi-helix', 'canonical')
IC_Val_HelixPi_gminus = build_template('VAL', 'pi-helix', 'g-')
IC_Val_HelixPi_trans = build_template('VAL', 'pi-helix', 't')
IC_Val_HelixPi_gplus = build_template('VAL', 'pi-helix', 'g+')
IC_Val_HelixPPII = build_template('VAL', 'polyproline-II', 'canonical')
IC_Val_HelixPPII_gminus = build_template('VAL', 'polyproline-II', 'g-')
IC_Val_HelixPPII_trans = build_template('VAL', 'polyproline-II', 't')
IC_Val_HelixPPII_gplus = build_template('VAL', 'polyproline-II', 'g+')
IC_Val_Strand = build_template('VAL', 'beta-strand', 'canonical')
IC_Val_Strand_gminus = build_template('VAL', 'beta-strand', 'g-')
IC_Val_Strand_trans = build_template('VAL', 'beta-strand', 't')
IC_Val_Strand_gplus = build_template('VAL', 'beta-strand', 'g+')
IC_Val_StrandParallel = build_template('VAL', 'parallel-beta-strand', 'canonical')
IC_Val_StrandParallel_gminus = build_template('VAL', 'parallel-beta-strand', 'g-')
IC_Val_StrandParallel_trans = build_template('VAL', 'parallel-beta-strand', 't')
IC_Val_StrandParallel_gplus = build_template('VAL', 'parallel-beta-strand', 'g+')
IC_Val_StrandAntiParallel = build_template('VAL', 'antiparallel-beta-strand', 'canonical')
IC_Val_StrandAntiParallel_gminus = build_template('VAL', 'antiparallel-beta-strand', 'g-')
IC_Val_StrandAntiParallel_trans = build_template('VAL', 'antiparallel-beta-strand', 't')
IC_Val_StrandAntiParallel_gplus = build_template('VAL', 'antiparallel-beta-strand', 'g+')
IC_Val_Bridge = build_template('VAL', 'beta-bridge', 'canonical')
IC_Val_Bridge_gminus = build_template('VAL', 'beta-bridge', 'g-')
IC_Val_Bridge_trans = build_template('VAL', 'beta-bridge', 't')
IC_Val_Bridge_gplus = build_template('VAL', 'beta-bridge', 'g+')
IC_Val_Turn = build_template('VAL', 'turn', 'canonical')
IC_Val_Turn_gminus = build_template('VAL', 'turn', 'g-')
IC_Val_Turn_trans = build_template('VAL', 'turn', 't')
IC_Val_Turn_gplus = build_template('VAL', 'turn', 'g+')
IC_Val_Bend = build_template('VAL', 'bend', 'canonical')
IC_Val_Bend_gminus = build_template('VAL', 'bend', 'g-')
IC_Val_Bend_trans = build_template('VAL', 'bend', 't')
IC_Val_Bend_gplus = build_template('VAL', 'bend', 'g+')
IC_Val_Coil = build_template('VAL', 'coil', 'canonical')
IC_Val_Coil_gminus = build_template('VAL', 'coil', 'g-')
IC_Val_Coil_trans = build_template('VAL', 'coil', 't')
IC_Val_Coil_gplus = build_template('VAL', 'coil', 'g+')
IC_Val_CisPeptide = build_template('VAL', 'cis-peptide-bond', 'canonical')
IC_Val_CisPeptide_gminus = build_template('VAL', 'cis-peptide-bond', 'g-')
IC_Val_CisPeptide_trans = build_template('VAL', 'cis-peptide-bond', 't')
IC_Val_CisPeptide_gplus = build_template('VAL', 'cis-peptide-bond', 'g+')

__all__ = [
    'IC_Val_Bend',
    'IC_Val_Bend_gminus',
    'IC_Val_Bend_gplus',
    'IC_Val_Bend_trans',
    'IC_Val_Bridge',
    'IC_Val_Bridge_gminus',
    'IC_Val_Bridge_gplus',
    'IC_Val_Bridge_trans',
    'IC_Val_CisPeptide',
    'IC_Val_CisPeptide_gminus',
    'IC_Val_CisPeptide_gplus',
    'IC_Val_CisPeptide_trans',
    'IC_Val_Coil',
    'IC_Val_Coil_gminus',
    'IC_Val_Coil_gplus',
    'IC_Val_Coil_trans',
    'IC_Val_Helix310',
    'IC_Val_Helix310_gminus',
    'IC_Val_Helix310_gplus',
    'IC_Val_Helix310_trans',
    'IC_Val_HelixAlpha',
    'IC_Val_HelixAlpha_gminus',
    'IC_Val_HelixAlpha_gplus',
    'IC_Val_HelixAlpha_trans',
    'IC_Val_HelixPPII',
    'IC_Val_HelixPPII_gminus',
    'IC_Val_HelixPPII_gplus',
    'IC_Val_HelixPPII_trans',
    'IC_Val_HelixPi',
    'IC_Val_HelixPi_gminus',
    'IC_Val_HelixPi_gplus',
    'IC_Val_HelixPi_trans',
    'IC_Val_Strand',
    'IC_Val_StrandAntiParallel',
    'IC_Val_StrandAntiParallel_gminus',
    'IC_Val_StrandAntiParallel_gplus',
    'IC_Val_StrandAntiParallel_trans',
    'IC_Val_StrandParallel',
    'IC_Val_StrandParallel_gminus',
    'IC_Val_StrandParallel_gplus',
    'IC_Val_StrandParallel_trans',
    'IC_Val_Strand_gminus',
    'IC_Val_Strand_gplus',
    'IC_Val_Strand_trans',
    'IC_Val_Turn',
    'IC_Val_Turn_gminus',
    'IC_Val_Turn_gplus',
    'IC_Val_Turn_trans',
]
