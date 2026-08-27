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

IC_Thr_HelixAlpha = build_template('THR', 'alpha-helix', 'canonical')
IC_Thr_HelixAlpha_gminus = build_template('THR', 'alpha-helix', 'g-')
IC_Thr_HelixAlpha_trans = build_template('THR', 'alpha-helix', 't')
IC_Thr_HelixAlpha_gplus = build_template('THR', 'alpha-helix', 'g+')
IC_Thr_Helix310 = build_template('THR', '3-10-helix', 'canonical')
IC_Thr_Helix310_gminus = build_template('THR', '3-10-helix', 'g-')
IC_Thr_Helix310_trans = build_template('THR', '3-10-helix', 't')
IC_Thr_Helix310_gplus = build_template('THR', '3-10-helix', 'g+')
IC_Thr_HelixPi = build_template('THR', 'pi-helix', 'canonical')
IC_Thr_HelixPi_gminus = build_template('THR', 'pi-helix', 'g-')
IC_Thr_HelixPi_trans = build_template('THR', 'pi-helix', 't')
IC_Thr_HelixPi_gplus = build_template('THR', 'pi-helix', 'g+')
IC_Thr_HelixPPII = build_template('THR', 'polyproline-II', 'canonical')
IC_Thr_HelixPPII_gminus = build_template('THR', 'polyproline-II', 'g-')
IC_Thr_HelixPPII_trans = build_template('THR', 'polyproline-II', 't')
IC_Thr_HelixPPII_gplus = build_template('THR', 'polyproline-II', 'g+')
IC_Thr_Strand = build_template('THR', 'beta-strand', 'canonical')
IC_Thr_Strand_gminus = build_template('THR', 'beta-strand', 'g-')
IC_Thr_Strand_trans = build_template('THR', 'beta-strand', 't')
IC_Thr_Strand_gplus = build_template('THR', 'beta-strand', 'g+')
IC_Thr_StrandParallel = build_template('THR', 'parallel-beta-strand', 'canonical')
IC_Thr_StrandParallel_gminus = build_template('THR', 'parallel-beta-strand', 'g-')
IC_Thr_StrandParallel_trans = build_template('THR', 'parallel-beta-strand', 't')
IC_Thr_StrandParallel_gplus = build_template('THR', 'parallel-beta-strand', 'g+')
IC_Thr_StrandAntiParallel = build_template('THR', 'antiparallel-beta-strand', 'canonical')
IC_Thr_StrandAntiParallel_gminus = build_template('THR', 'antiparallel-beta-strand', 'g-')
IC_Thr_StrandAntiParallel_trans = build_template('THR', 'antiparallel-beta-strand', 't')
IC_Thr_StrandAntiParallel_gplus = build_template('THR', 'antiparallel-beta-strand', 'g+')
IC_Thr_Bridge = build_template('THR', 'beta-bridge', 'canonical')
IC_Thr_Bridge_gminus = build_template('THR', 'beta-bridge', 'g-')
IC_Thr_Bridge_trans = build_template('THR', 'beta-bridge', 't')
IC_Thr_Bridge_gplus = build_template('THR', 'beta-bridge', 'g+')
IC_Thr_Turn = build_template('THR', 'turn', 'canonical')
IC_Thr_Turn_gminus = build_template('THR', 'turn', 'g-')
IC_Thr_Turn_trans = build_template('THR', 'turn', 't')
IC_Thr_Turn_gplus = build_template('THR', 'turn', 'g+')
IC_Thr_Bend = build_template('THR', 'bend', 'canonical')
IC_Thr_Bend_gminus = build_template('THR', 'bend', 'g-')
IC_Thr_Bend_trans = build_template('THR', 'bend', 't')
IC_Thr_Bend_gplus = build_template('THR', 'bend', 'g+')
IC_Thr_Coil = build_template('THR', 'coil', 'canonical')
IC_Thr_Coil_gminus = build_template('THR', 'coil', 'g-')
IC_Thr_Coil_trans = build_template('THR', 'coil', 't')
IC_Thr_Coil_gplus = build_template('THR', 'coil', 'g+')
IC_Thr_CisPeptide = build_template('THR', 'cis-peptide-bond', 'canonical')
IC_Thr_CisPeptide_gminus = build_template('THR', 'cis-peptide-bond', 'g-')
IC_Thr_CisPeptide_trans = build_template('THR', 'cis-peptide-bond', 't')
IC_Thr_CisPeptide_gplus = build_template('THR', 'cis-peptide-bond', 'g+')

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
