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

IC_Trp_HelixAlpha = build_template('TRP', 'alpha-helix', 'canonical')
IC_Trp_HelixAlpha_gminus_gminus = build_template('TRP', 'alpha-helix', 'g-/g-')
IC_Trp_HelixAlpha_gminus_trans = build_template('TRP', 'alpha-helix', 'g-/t')
IC_Trp_HelixAlpha_gminus_gplus = build_template('TRP', 'alpha-helix', 'g-/g+')
IC_Trp_HelixAlpha_trans_gminus = build_template('TRP', 'alpha-helix', 't/g-')
IC_Trp_HelixAlpha_trans_trans = build_template('TRP', 'alpha-helix', 't/t')
IC_Trp_HelixAlpha_trans_gplus = build_template('TRP', 'alpha-helix', 't/g+')
IC_Trp_HelixAlpha_gplus_gminus = build_template('TRP', 'alpha-helix', 'g+/g-')
IC_Trp_HelixAlpha_gplus_trans = build_template('TRP', 'alpha-helix', 'g+/t')
IC_Trp_HelixAlpha_gplus_gplus = build_template('TRP', 'alpha-helix', 'g+/g+')
IC_Trp_Helix310 = build_template('TRP', '3-10-helix', 'canonical')
IC_Trp_Helix310_gminus_gminus = build_template('TRP', '3-10-helix', 'g-/g-')
IC_Trp_Helix310_gminus_trans = build_template('TRP', '3-10-helix', 'g-/t')
IC_Trp_Helix310_gminus_gplus = build_template('TRP', '3-10-helix', 'g-/g+')
IC_Trp_Helix310_trans_gminus = build_template('TRP', '3-10-helix', 't/g-')
IC_Trp_Helix310_trans_trans = build_template('TRP', '3-10-helix', 't/t')
IC_Trp_Helix310_trans_gplus = build_template('TRP', '3-10-helix', 't/g+')
IC_Trp_Helix310_gplus_gminus = build_template('TRP', '3-10-helix', 'g+/g-')
IC_Trp_Helix310_gplus_trans = build_template('TRP', '3-10-helix', 'g+/t')
IC_Trp_Helix310_gplus_gplus = build_template('TRP', '3-10-helix', 'g+/g+')
IC_Trp_HelixPi = build_template('TRP', 'pi-helix', 'canonical')
IC_Trp_HelixPi_gminus_gminus = build_template('TRP', 'pi-helix', 'g-/g-')
IC_Trp_HelixPi_gminus_trans = build_template('TRP', 'pi-helix', 'g-/t')
IC_Trp_HelixPi_gminus_gplus = build_template('TRP', 'pi-helix', 'g-/g+')
IC_Trp_HelixPi_trans_gminus = build_template('TRP', 'pi-helix', 't/g-')
IC_Trp_HelixPi_trans_trans = build_template('TRP', 'pi-helix', 't/t')
IC_Trp_HelixPi_trans_gplus = build_template('TRP', 'pi-helix', 't/g+')
IC_Trp_HelixPi_gplus_gminus = build_template('TRP', 'pi-helix', 'g+/g-')
IC_Trp_HelixPi_gplus_trans = build_template('TRP', 'pi-helix', 'g+/t')
IC_Trp_HelixPi_gplus_gplus = build_template('TRP', 'pi-helix', 'g+/g+')
IC_Trp_HelixPPII = build_template('TRP', 'polyproline-II', 'canonical')
IC_Trp_HelixPPII_gminus_gminus = build_template('TRP', 'polyproline-II', 'g-/g-')
IC_Trp_HelixPPII_gminus_trans = build_template('TRP', 'polyproline-II', 'g-/t')
IC_Trp_HelixPPII_gminus_gplus = build_template('TRP', 'polyproline-II', 'g-/g+')
IC_Trp_HelixPPII_trans_gminus = build_template('TRP', 'polyproline-II', 't/g-')
IC_Trp_HelixPPII_trans_trans = build_template('TRP', 'polyproline-II', 't/t')
IC_Trp_HelixPPII_trans_gplus = build_template('TRP', 'polyproline-II', 't/g+')
IC_Trp_HelixPPII_gplus_gminus = build_template('TRP', 'polyproline-II', 'g+/g-')
IC_Trp_HelixPPII_gplus_trans = build_template('TRP', 'polyproline-II', 'g+/t')
IC_Trp_HelixPPII_gplus_gplus = build_template('TRP', 'polyproline-II', 'g+/g+')
IC_Trp_Strand = build_template('TRP', 'beta-strand', 'canonical')
IC_Trp_Strand_gminus_gminus = build_template('TRP', 'beta-strand', 'g-/g-')
IC_Trp_Strand_gminus_trans = build_template('TRP', 'beta-strand', 'g-/t')
IC_Trp_Strand_gminus_gplus = build_template('TRP', 'beta-strand', 'g-/g+')
IC_Trp_Strand_trans_gminus = build_template('TRP', 'beta-strand', 't/g-')
IC_Trp_Strand_trans_trans = build_template('TRP', 'beta-strand', 't/t')
IC_Trp_Strand_trans_gplus = build_template('TRP', 'beta-strand', 't/g+')
IC_Trp_Strand_gplus_gminus = build_template('TRP', 'beta-strand', 'g+/g-')
IC_Trp_Strand_gplus_trans = build_template('TRP', 'beta-strand', 'g+/t')
IC_Trp_Strand_gplus_gplus = build_template('TRP', 'beta-strand', 'g+/g+')
IC_Trp_StrandParallel = build_template('TRP', 'parallel-beta-strand', 'canonical')
IC_Trp_StrandParallel_gminus_gminus = build_template('TRP', 'parallel-beta-strand', 'g-/g-')
IC_Trp_StrandParallel_gminus_trans = build_template('TRP', 'parallel-beta-strand', 'g-/t')
IC_Trp_StrandParallel_gminus_gplus = build_template('TRP', 'parallel-beta-strand', 'g-/g+')
IC_Trp_StrandParallel_trans_gminus = build_template('TRP', 'parallel-beta-strand', 't/g-')
IC_Trp_StrandParallel_trans_trans = build_template('TRP', 'parallel-beta-strand', 't/t')
IC_Trp_StrandParallel_trans_gplus = build_template('TRP', 'parallel-beta-strand', 't/g+')
IC_Trp_StrandParallel_gplus_gminus = build_template('TRP', 'parallel-beta-strand', 'g+/g-')
IC_Trp_StrandParallel_gplus_trans = build_template('TRP', 'parallel-beta-strand', 'g+/t')
IC_Trp_StrandParallel_gplus_gplus = build_template('TRP', 'parallel-beta-strand', 'g+/g+')
IC_Trp_StrandAntiParallel = build_template('TRP', 'antiparallel-beta-strand', 'canonical')
IC_Trp_StrandAntiParallel_gminus_gminus = build_template('TRP', 'antiparallel-beta-strand', 'g-/g-')
IC_Trp_StrandAntiParallel_gminus_trans = build_template('TRP', 'antiparallel-beta-strand', 'g-/t')
IC_Trp_StrandAntiParallel_gminus_gplus = build_template('TRP', 'antiparallel-beta-strand', 'g-/g+')
IC_Trp_StrandAntiParallel_trans_gminus = build_template('TRP', 'antiparallel-beta-strand', 't/g-')
IC_Trp_StrandAntiParallel_trans_trans = build_template('TRP', 'antiparallel-beta-strand', 't/t')
IC_Trp_StrandAntiParallel_trans_gplus = build_template('TRP', 'antiparallel-beta-strand', 't/g+')
IC_Trp_StrandAntiParallel_gplus_gminus = build_template('TRP', 'antiparallel-beta-strand', 'g+/g-')
IC_Trp_StrandAntiParallel_gplus_trans = build_template('TRP', 'antiparallel-beta-strand', 'g+/t')
IC_Trp_StrandAntiParallel_gplus_gplus = build_template('TRP', 'antiparallel-beta-strand', 'g+/g+')
IC_Trp_Bridge = build_template('TRP', 'beta-bridge', 'canonical')
IC_Trp_Bridge_gminus_gminus = build_template('TRP', 'beta-bridge', 'g-/g-')
IC_Trp_Bridge_gminus_trans = build_template('TRP', 'beta-bridge', 'g-/t')
IC_Trp_Bridge_gminus_gplus = build_template('TRP', 'beta-bridge', 'g-/g+')
IC_Trp_Bridge_trans_gminus = build_template('TRP', 'beta-bridge', 't/g-')
IC_Trp_Bridge_trans_trans = build_template('TRP', 'beta-bridge', 't/t')
IC_Trp_Bridge_trans_gplus = build_template('TRP', 'beta-bridge', 't/g+')
IC_Trp_Bridge_gplus_gminus = build_template('TRP', 'beta-bridge', 'g+/g-')
IC_Trp_Bridge_gplus_trans = build_template('TRP', 'beta-bridge', 'g+/t')
IC_Trp_Bridge_gplus_gplus = build_template('TRP', 'beta-bridge', 'g+/g+')
IC_Trp_Turn = build_template('TRP', 'turn', 'canonical')
IC_Trp_Turn_gminus_gminus = build_template('TRP', 'turn', 'g-/g-')
IC_Trp_Turn_gminus_trans = build_template('TRP', 'turn', 'g-/t')
IC_Trp_Turn_gminus_gplus = build_template('TRP', 'turn', 'g-/g+')
IC_Trp_Turn_trans_gminus = build_template('TRP', 'turn', 't/g-')
IC_Trp_Turn_trans_trans = build_template('TRP', 'turn', 't/t')
IC_Trp_Turn_trans_gplus = build_template('TRP', 'turn', 't/g+')
IC_Trp_Turn_gplus_gminus = build_template('TRP', 'turn', 'g+/g-')
IC_Trp_Turn_gplus_trans = build_template('TRP', 'turn', 'g+/t')
IC_Trp_Turn_gplus_gplus = build_template('TRP', 'turn', 'g+/g+')
IC_Trp_Bend = build_template('TRP', 'bend', 'canonical')
IC_Trp_Bend_gminus_gminus = build_template('TRP', 'bend', 'g-/g-')
IC_Trp_Bend_gminus_trans = build_template('TRP', 'bend', 'g-/t')
IC_Trp_Bend_gminus_gplus = build_template('TRP', 'bend', 'g-/g+')
IC_Trp_Bend_trans_gminus = build_template('TRP', 'bend', 't/g-')
IC_Trp_Bend_trans_trans = build_template('TRP', 'bend', 't/t')
IC_Trp_Bend_trans_gplus = build_template('TRP', 'bend', 't/g+')
IC_Trp_Bend_gplus_gminus = build_template('TRP', 'bend', 'g+/g-')
IC_Trp_Bend_gplus_trans = build_template('TRP', 'bend', 'g+/t')
IC_Trp_Bend_gplus_gplus = build_template('TRP', 'bend', 'g+/g+')
IC_Trp_Coil = build_template('TRP', 'coil', 'canonical')
IC_Trp_Coil_gminus_gminus = build_template('TRP', 'coil', 'g-/g-')
IC_Trp_Coil_gminus_trans = build_template('TRP', 'coil', 'g-/t')
IC_Trp_Coil_gminus_gplus = build_template('TRP', 'coil', 'g-/g+')
IC_Trp_Coil_trans_gminus = build_template('TRP', 'coil', 't/g-')
IC_Trp_Coil_trans_trans = build_template('TRP', 'coil', 't/t')
IC_Trp_Coil_trans_gplus = build_template('TRP', 'coil', 't/g+')
IC_Trp_Coil_gplus_gminus = build_template('TRP', 'coil', 'g+/g-')
IC_Trp_Coil_gplus_trans = build_template('TRP', 'coil', 'g+/t')
IC_Trp_Coil_gplus_gplus = build_template('TRP', 'coil', 'g+/g+')
IC_Trp_CisPeptide = build_template('TRP', 'cis-peptide-bond', 'canonical')
IC_Trp_CisPeptide_gminus_gminus = build_template('TRP', 'cis-peptide-bond', 'g-/g-')
IC_Trp_CisPeptide_gminus_trans = build_template('TRP', 'cis-peptide-bond', 'g-/t')
IC_Trp_CisPeptide_gminus_gplus = build_template('TRP', 'cis-peptide-bond', 'g-/g+')
IC_Trp_CisPeptide_trans_gminus = build_template('TRP', 'cis-peptide-bond', 't/g-')
IC_Trp_CisPeptide_trans_trans = build_template('TRP', 'cis-peptide-bond', 't/t')
IC_Trp_CisPeptide_trans_gplus = build_template('TRP', 'cis-peptide-bond', 't/g+')
IC_Trp_CisPeptide_gplus_gminus = build_template('TRP', 'cis-peptide-bond', 'g+/g-')
IC_Trp_CisPeptide_gplus_trans = build_template('TRP', 'cis-peptide-bond', 'g+/t')
IC_Trp_CisPeptide_gplus_gplus = build_template('TRP', 'cis-peptide-bond', 'g+/g+')

__all__ = [
    'IC_Trp_Bend',
    'IC_Trp_Bend_gminus_gminus',
    'IC_Trp_Bend_gminus_gplus',
    'IC_Trp_Bend_gminus_trans',
    'IC_Trp_Bend_gplus_gminus',
    'IC_Trp_Bend_gplus_gplus',
    'IC_Trp_Bend_gplus_trans',
    'IC_Trp_Bend_trans_gminus',
    'IC_Trp_Bend_trans_gplus',
    'IC_Trp_Bend_trans_trans',
    'IC_Trp_Bridge',
    'IC_Trp_Bridge_gminus_gminus',
    'IC_Trp_Bridge_gminus_gplus',
    'IC_Trp_Bridge_gminus_trans',
    'IC_Trp_Bridge_gplus_gminus',
    'IC_Trp_Bridge_gplus_gplus',
    'IC_Trp_Bridge_gplus_trans',
    'IC_Trp_Bridge_trans_gminus',
    'IC_Trp_Bridge_trans_gplus',
    'IC_Trp_Bridge_trans_trans',
    'IC_Trp_CisPeptide',
    'IC_Trp_CisPeptide_gminus_gminus',
    'IC_Trp_CisPeptide_gminus_gplus',
    'IC_Trp_CisPeptide_gminus_trans',
    'IC_Trp_CisPeptide_gplus_gminus',
    'IC_Trp_CisPeptide_gplus_gplus',
    'IC_Trp_CisPeptide_gplus_trans',
    'IC_Trp_CisPeptide_trans_gminus',
    'IC_Trp_CisPeptide_trans_gplus',
    'IC_Trp_CisPeptide_trans_trans',
    'IC_Trp_Coil',
    'IC_Trp_Coil_gminus_gminus',
    'IC_Trp_Coil_gminus_gplus',
    'IC_Trp_Coil_gminus_trans',
    'IC_Trp_Coil_gplus_gminus',
    'IC_Trp_Coil_gplus_gplus',
    'IC_Trp_Coil_gplus_trans',
    'IC_Trp_Coil_trans_gminus',
    'IC_Trp_Coil_trans_gplus',
    'IC_Trp_Coil_trans_trans',
    'IC_Trp_Helix310',
    'IC_Trp_Helix310_gminus_gminus',
    'IC_Trp_Helix310_gminus_gplus',
    'IC_Trp_Helix310_gminus_trans',
    'IC_Trp_Helix310_gplus_gminus',
    'IC_Trp_Helix310_gplus_gplus',
    'IC_Trp_Helix310_gplus_trans',
    'IC_Trp_Helix310_trans_gminus',
    'IC_Trp_Helix310_trans_gplus',
    'IC_Trp_Helix310_trans_trans',
    'IC_Trp_HelixAlpha',
    'IC_Trp_HelixAlpha_gminus_gminus',
    'IC_Trp_HelixAlpha_gminus_gplus',
    'IC_Trp_HelixAlpha_gminus_trans',
    'IC_Trp_HelixAlpha_gplus_gminus',
    'IC_Trp_HelixAlpha_gplus_gplus',
    'IC_Trp_HelixAlpha_gplus_trans',
    'IC_Trp_HelixAlpha_trans_gminus',
    'IC_Trp_HelixAlpha_trans_gplus',
    'IC_Trp_HelixAlpha_trans_trans',
    'IC_Trp_HelixPPII',
    'IC_Trp_HelixPPII_gminus_gminus',
    'IC_Trp_HelixPPII_gminus_gplus',
    'IC_Trp_HelixPPII_gminus_trans',
    'IC_Trp_HelixPPII_gplus_gminus',
    'IC_Trp_HelixPPII_gplus_gplus',
    'IC_Trp_HelixPPII_gplus_trans',
    'IC_Trp_HelixPPII_trans_gminus',
    'IC_Trp_HelixPPII_trans_gplus',
    'IC_Trp_HelixPPII_trans_trans',
    'IC_Trp_HelixPi',
    'IC_Trp_HelixPi_gminus_gminus',
    'IC_Trp_HelixPi_gminus_gplus',
    'IC_Trp_HelixPi_gminus_trans',
    'IC_Trp_HelixPi_gplus_gminus',
    'IC_Trp_HelixPi_gplus_gplus',
    'IC_Trp_HelixPi_gplus_trans',
    'IC_Trp_HelixPi_trans_gminus',
    'IC_Trp_HelixPi_trans_gplus',
    'IC_Trp_HelixPi_trans_trans',
    'IC_Trp_Strand',
    'IC_Trp_StrandAntiParallel',
    'IC_Trp_StrandAntiParallel_gminus_gminus',
    'IC_Trp_StrandAntiParallel_gminus_gplus',
    'IC_Trp_StrandAntiParallel_gminus_trans',
    'IC_Trp_StrandAntiParallel_gplus_gminus',
    'IC_Trp_StrandAntiParallel_gplus_gplus',
    'IC_Trp_StrandAntiParallel_gplus_trans',
    'IC_Trp_StrandAntiParallel_trans_gminus',
    'IC_Trp_StrandAntiParallel_trans_gplus',
    'IC_Trp_StrandAntiParallel_trans_trans',
    'IC_Trp_StrandParallel',
    'IC_Trp_StrandParallel_gminus_gminus',
    'IC_Trp_StrandParallel_gminus_gplus',
    'IC_Trp_StrandParallel_gminus_trans',
    'IC_Trp_StrandParallel_gplus_gminus',
    'IC_Trp_StrandParallel_gplus_gplus',
    'IC_Trp_StrandParallel_gplus_trans',
    'IC_Trp_StrandParallel_trans_gminus',
    'IC_Trp_StrandParallel_trans_gplus',
    'IC_Trp_StrandParallel_trans_trans',
    'IC_Trp_Strand_gminus_gminus',
    'IC_Trp_Strand_gminus_gplus',
    'IC_Trp_Strand_gminus_trans',
    'IC_Trp_Strand_gplus_gminus',
    'IC_Trp_Strand_gplus_gplus',
    'IC_Trp_Strand_gplus_trans',
    'IC_Trp_Strand_trans_gminus',
    'IC_Trp_Strand_trans_gplus',
    'IC_Trp_Strand_trans_trans',
    'IC_Trp_Turn',
    'IC_Trp_Turn_gminus_gminus',
    'IC_Trp_Turn_gminus_gplus',
    'IC_Trp_Turn_gminus_trans',
    'IC_Trp_Turn_gplus_gminus',
    'IC_Trp_Turn_gplus_gplus',
    'IC_Trp_Turn_gplus_trans',
    'IC_Trp_Turn_trans_gminus',
    'IC_Trp_Turn_trans_gplus',
    'IC_Trp_Turn_trans_trans',
]
