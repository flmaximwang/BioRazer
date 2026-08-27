# -*- coding: utf-8 -*-
"""Ideal {secondary-structure x rotamer} {names[aa]} {var} ``InternalCoord`` templates.

Each template is the ideal heavy-atom topology for {aa} anchored at
{{N, CA, C}}, with the side chain grown to ``rotamer`` and the
``phi``/``psi``/``omega`` of ``ss`` carried on the instance.  See
:mod:`biorazer.database.molecule.icoor.protein.template` and
:mod:`._builder` for the build rule and the per-conformer caveat.
"""

from biorazer.database.molecule.icoor.protein.template._builder import build_template

RESN = "TYR"

# === GENERATED IC_* named templates - do not edit below ===
# Regenerate with: /opt/envs/BioRazer/bin/python scripts/generate_internal_coord_template_named.py

IC_Tyr_HelixAlpha = build_template('TYR', 'alpha-helix', 'canonical')
IC_Tyr_HelixAlpha_gminus_gminus = build_template('TYR', 'alpha-helix', 'g-/g-')
IC_Tyr_HelixAlpha_gminus_trans = build_template('TYR', 'alpha-helix', 'g-/t')
IC_Tyr_HelixAlpha_gminus_gplus = build_template('TYR', 'alpha-helix', 'g-/g+')
IC_Tyr_HelixAlpha_trans_gminus = build_template('TYR', 'alpha-helix', 't/g-')
IC_Tyr_HelixAlpha_trans_trans = build_template('TYR', 'alpha-helix', 't/t')
IC_Tyr_HelixAlpha_trans_gplus = build_template('TYR', 'alpha-helix', 't/g+')
IC_Tyr_HelixAlpha_gplus_gminus = build_template('TYR', 'alpha-helix', 'g+/g-')
IC_Tyr_HelixAlpha_gplus_trans = build_template('TYR', 'alpha-helix', 'g+/t')
IC_Tyr_HelixAlpha_gplus_gplus = build_template('TYR', 'alpha-helix', 'g+/g+')
IC_Tyr_Helix310 = build_template('TYR', '3-10-helix', 'canonical')
IC_Tyr_Helix310_gminus_gminus = build_template('TYR', '3-10-helix', 'g-/g-')
IC_Tyr_Helix310_gminus_trans = build_template('TYR', '3-10-helix', 'g-/t')
IC_Tyr_Helix310_gminus_gplus = build_template('TYR', '3-10-helix', 'g-/g+')
IC_Tyr_Helix310_trans_gminus = build_template('TYR', '3-10-helix', 't/g-')
IC_Tyr_Helix310_trans_trans = build_template('TYR', '3-10-helix', 't/t')
IC_Tyr_Helix310_trans_gplus = build_template('TYR', '3-10-helix', 't/g+')
IC_Tyr_Helix310_gplus_gminus = build_template('TYR', '3-10-helix', 'g+/g-')
IC_Tyr_Helix310_gplus_trans = build_template('TYR', '3-10-helix', 'g+/t')
IC_Tyr_Helix310_gplus_gplus = build_template('TYR', '3-10-helix', 'g+/g+')
IC_Tyr_HelixPi = build_template('TYR', 'pi-helix', 'canonical')
IC_Tyr_HelixPi_gminus_gminus = build_template('TYR', 'pi-helix', 'g-/g-')
IC_Tyr_HelixPi_gminus_trans = build_template('TYR', 'pi-helix', 'g-/t')
IC_Tyr_HelixPi_gminus_gplus = build_template('TYR', 'pi-helix', 'g-/g+')
IC_Tyr_HelixPi_trans_gminus = build_template('TYR', 'pi-helix', 't/g-')
IC_Tyr_HelixPi_trans_trans = build_template('TYR', 'pi-helix', 't/t')
IC_Tyr_HelixPi_trans_gplus = build_template('TYR', 'pi-helix', 't/g+')
IC_Tyr_HelixPi_gplus_gminus = build_template('TYR', 'pi-helix', 'g+/g-')
IC_Tyr_HelixPi_gplus_trans = build_template('TYR', 'pi-helix', 'g+/t')
IC_Tyr_HelixPi_gplus_gplus = build_template('TYR', 'pi-helix', 'g+/g+')
IC_Tyr_HelixPPII = build_template('TYR', 'polyproline-II', 'canonical')
IC_Tyr_HelixPPII_gminus_gminus = build_template('TYR', 'polyproline-II', 'g-/g-')
IC_Tyr_HelixPPII_gminus_trans = build_template('TYR', 'polyproline-II', 'g-/t')
IC_Tyr_HelixPPII_gminus_gplus = build_template('TYR', 'polyproline-II', 'g-/g+')
IC_Tyr_HelixPPII_trans_gminus = build_template('TYR', 'polyproline-II', 't/g-')
IC_Tyr_HelixPPII_trans_trans = build_template('TYR', 'polyproline-II', 't/t')
IC_Tyr_HelixPPII_trans_gplus = build_template('TYR', 'polyproline-II', 't/g+')
IC_Tyr_HelixPPII_gplus_gminus = build_template('TYR', 'polyproline-II', 'g+/g-')
IC_Tyr_HelixPPII_gplus_trans = build_template('TYR', 'polyproline-II', 'g+/t')
IC_Tyr_HelixPPII_gplus_gplus = build_template('TYR', 'polyproline-II', 'g+/g+')
IC_Tyr_Strand = build_template('TYR', 'beta-strand', 'canonical')
IC_Tyr_Strand_gminus_gminus = build_template('TYR', 'beta-strand', 'g-/g-')
IC_Tyr_Strand_gminus_trans = build_template('TYR', 'beta-strand', 'g-/t')
IC_Tyr_Strand_gminus_gplus = build_template('TYR', 'beta-strand', 'g-/g+')
IC_Tyr_Strand_trans_gminus = build_template('TYR', 'beta-strand', 't/g-')
IC_Tyr_Strand_trans_trans = build_template('TYR', 'beta-strand', 't/t')
IC_Tyr_Strand_trans_gplus = build_template('TYR', 'beta-strand', 't/g+')
IC_Tyr_Strand_gplus_gminus = build_template('TYR', 'beta-strand', 'g+/g-')
IC_Tyr_Strand_gplus_trans = build_template('TYR', 'beta-strand', 'g+/t')
IC_Tyr_Strand_gplus_gplus = build_template('TYR', 'beta-strand', 'g+/g+')
IC_Tyr_StrandParallel = build_template('TYR', 'parallel-beta-strand', 'canonical')
IC_Tyr_StrandParallel_gminus_gminus = build_template('TYR', 'parallel-beta-strand', 'g-/g-')
IC_Tyr_StrandParallel_gminus_trans = build_template('TYR', 'parallel-beta-strand', 'g-/t')
IC_Tyr_StrandParallel_gminus_gplus = build_template('TYR', 'parallel-beta-strand', 'g-/g+')
IC_Tyr_StrandParallel_trans_gminus = build_template('TYR', 'parallel-beta-strand', 't/g-')
IC_Tyr_StrandParallel_trans_trans = build_template('TYR', 'parallel-beta-strand', 't/t')
IC_Tyr_StrandParallel_trans_gplus = build_template('TYR', 'parallel-beta-strand', 't/g+')
IC_Tyr_StrandParallel_gplus_gminus = build_template('TYR', 'parallel-beta-strand', 'g+/g-')
IC_Tyr_StrandParallel_gplus_trans = build_template('TYR', 'parallel-beta-strand', 'g+/t')
IC_Tyr_StrandParallel_gplus_gplus = build_template('TYR', 'parallel-beta-strand', 'g+/g+')
IC_Tyr_StrandAntiParallel = build_template('TYR', 'antiparallel-beta-strand', 'canonical')
IC_Tyr_StrandAntiParallel_gminus_gminus = build_template('TYR', 'antiparallel-beta-strand', 'g-/g-')
IC_Tyr_StrandAntiParallel_gminus_trans = build_template('TYR', 'antiparallel-beta-strand', 'g-/t')
IC_Tyr_StrandAntiParallel_gminus_gplus = build_template('TYR', 'antiparallel-beta-strand', 'g-/g+')
IC_Tyr_StrandAntiParallel_trans_gminus = build_template('TYR', 'antiparallel-beta-strand', 't/g-')
IC_Tyr_StrandAntiParallel_trans_trans = build_template('TYR', 'antiparallel-beta-strand', 't/t')
IC_Tyr_StrandAntiParallel_trans_gplus = build_template('TYR', 'antiparallel-beta-strand', 't/g+')
IC_Tyr_StrandAntiParallel_gplus_gminus = build_template('TYR', 'antiparallel-beta-strand', 'g+/g-')
IC_Tyr_StrandAntiParallel_gplus_trans = build_template('TYR', 'antiparallel-beta-strand', 'g+/t')
IC_Tyr_StrandAntiParallel_gplus_gplus = build_template('TYR', 'antiparallel-beta-strand', 'g+/g+')
IC_Tyr_Bridge = build_template('TYR', 'beta-bridge', 'canonical')
IC_Tyr_Bridge_gminus_gminus = build_template('TYR', 'beta-bridge', 'g-/g-')
IC_Tyr_Bridge_gminus_trans = build_template('TYR', 'beta-bridge', 'g-/t')
IC_Tyr_Bridge_gminus_gplus = build_template('TYR', 'beta-bridge', 'g-/g+')
IC_Tyr_Bridge_trans_gminus = build_template('TYR', 'beta-bridge', 't/g-')
IC_Tyr_Bridge_trans_trans = build_template('TYR', 'beta-bridge', 't/t')
IC_Tyr_Bridge_trans_gplus = build_template('TYR', 'beta-bridge', 't/g+')
IC_Tyr_Bridge_gplus_gminus = build_template('TYR', 'beta-bridge', 'g+/g-')
IC_Tyr_Bridge_gplus_trans = build_template('TYR', 'beta-bridge', 'g+/t')
IC_Tyr_Bridge_gplus_gplus = build_template('TYR', 'beta-bridge', 'g+/g+')
IC_Tyr_Turn = build_template('TYR', 'turn', 'canonical')
IC_Tyr_Turn_gminus_gminus = build_template('TYR', 'turn', 'g-/g-')
IC_Tyr_Turn_gminus_trans = build_template('TYR', 'turn', 'g-/t')
IC_Tyr_Turn_gminus_gplus = build_template('TYR', 'turn', 'g-/g+')
IC_Tyr_Turn_trans_gminus = build_template('TYR', 'turn', 't/g-')
IC_Tyr_Turn_trans_trans = build_template('TYR', 'turn', 't/t')
IC_Tyr_Turn_trans_gplus = build_template('TYR', 'turn', 't/g+')
IC_Tyr_Turn_gplus_gminus = build_template('TYR', 'turn', 'g+/g-')
IC_Tyr_Turn_gplus_trans = build_template('TYR', 'turn', 'g+/t')
IC_Tyr_Turn_gplus_gplus = build_template('TYR', 'turn', 'g+/g+')
IC_Tyr_Bend = build_template('TYR', 'bend', 'canonical')
IC_Tyr_Bend_gminus_gminus = build_template('TYR', 'bend', 'g-/g-')
IC_Tyr_Bend_gminus_trans = build_template('TYR', 'bend', 'g-/t')
IC_Tyr_Bend_gminus_gplus = build_template('TYR', 'bend', 'g-/g+')
IC_Tyr_Bend_trans_gminus = build_template('TYR', 'bend', 't/g-')
IC_Tyr_Bend_trans_trans = build_template('TYR', 'bend', 't/t')
IC_Tyr_Bend_trans_gplus = build_template('TYR', 'bend', 't/g+')
IC_Tyr_Bend_gplus_gminus = build_template('TYR', 'bend', 'g+/g-')
IC_Tyr_Bend_gplus_trans = build_template('TYR', 'bend', 'g+/t')
IC_Tyr_Bend_gplus_gplus = build_template('TYR', 'bend', 'g+/g+')
IC_Tyr_Coil = build_template('TYR', 'coil', 'canonical')
IC_Tyr_Coil_gminus_gminus = build_template('TYR', 'coil', 'g-/g-')
IC_Tyr_Coil_gminus_trans = build_template('TYR', 'coil', 'g-/t')
IC_Tyr_Coil_gminus_gplus = build_template('TYR', 'coil', 'g-/g+')
IC_Tyr_Coil_trans_gminus = build_template('TYR', 'coil', 't/g-')
IC_Tyr_Coil_trans_trans = build_template('TYR', 'coil', 't/t')
IC_Tyr_Coil_trans_gplus = build_template('TYR', 'coil', 't/g+')
IC_Tyr_Coil_gplus_gminus = build_template('TYR', 'coil', 'g+/g-')
IC_Tyr_Coil_gplus_trans = build_template('TYR', 'coil', 'g+/t')
IC_Tyr_Coil_gplus_gplus = build_template('TYR', 'coil', 'g+/g+')
IC_Tyr_CisPeptide = build_template('TYR', 'cis-peptide-bond', 'canonical')
IC_Tyr_CisPeptide_gminus_gminus = build_template('TYR', 'cis-peptide-bond', 'g-/g-')
IC_Tyr_CisPeptide_gminus_trans = build_template('TYR', 'cis-peptide-bond', 'g-/t')
IC_Tyr_CisPeptide_gminus_gplus = build_template('TYR', 'cis-peptide-bond', 'g-/g+')
IC_Tyr_CisPeptide_trans_gminus = build_template('TYR', 'cis-peptide-bond', 't/g-')
IC_Tyr_CisPeptide_trans_trans = build_template('TYR', 'cis-peptide-bond', 't/t')
IC_Tyr_CisPeptide_trans_gplus = build_template('TYR', 'cis-peptide-bond', 't/g+')
IC_Tyr_CisPeptide_gplus_gminus = build_template('TYR', 'cis-peptide-bond', 'g+/g-')
IC_Tyr_CisPeptide_gplus_trans = build_template('TYR', 'cis-peptide-bond', 'g+/t')
IC_Tyr_CisPeptide_gplus_gplus = build_template('TYR', 'cis-peptide-bond', 'g+/g+')

__all__ = [
    'IC_Tyr_Bend',
    'IC_Tyr_Bend_gminus_gminus',
    'IC_Tyr_Bend_gminus_gplus',
    'IC_Tyr_Bend_gminus_trans',
    'IC_Tyr_Bend_gplus_gminus',
    'IC_Tyr_Bend_gplus_gplus',
    'IC_Tyr_Bend_gplus_trans',
    'IC_Tyr_Bend_trans_gminus',
    'IC_Tyr_Bend_trans_gplus',
    'IC_Tyr_Bend_trans_trans',
    'IC_Tyr_Bridge',
    'IC_Tyr_Bridge_gminus_gminus',
    'IC_Tyr_Bridge_gminus_gplus',
    'IC_Tyr_Bridge_gminus_trans',
    'IC_Tyr_Bridge_gplus_gminus',
    'IC_Tyr_Bridge_gplus_gplus',
    'IC_Tyr_Bridge_gplus_trans',
    'IC_Tyr_Bridge_trans_gminus',
    'IC_Tyr_Bridge_trans_gplus',
    'IC_Tyr_Bridge_trans_trans',
    'IC_Tyr_CisPeptide',
    'IC_Tyr_CisPeptide_gminus_gminus',
    'IC_Tyr_CisPeptide_gminus_gplus',
    'IC_Tyr_CisPeptide_gminus_trans',
    'IC_Tyr_CisPeptide_gplus_gminus',
    'IC_Tyr_CisPeptide_gplus_gplus',
    'IC_Tyr_CisPeptide_gplus_trans',
    'IC_Tyr_CisPeptide_trans_gminus',
    'IC_Tyr_CisPeptide_trans_gplus',
    'IC_Tyr_CisPeptide_trans_trans',
    'IC_Tyr_Coil',
    'IC_Tyr_Coil_gminus_gminus',
    'IC_Tyr_Coil_gminus_gplus',
    'IC_Tyr_Coil_gminus_trans',
    'IC_Tyr_Coil_gplus_gminus',
    'IC_Tyr_Coil_gplus_gplus',
    'IC_Tyr_Coil_gplus_trans',
    'IC_Tyr_Coil_trans_gminus',
    'IC_Tyr_Coil_trans_gplus',
    'IC_Tyr_Coil_trans_trans',
    'IC_Tyr_Helix310',
    'IC_Tyr_Helix310_gminus_gminus',
    'IC_Tyr_Helix310_gminus_gplus',
    'IC_Tyr_Helix310_gminus_trans',
    'IC_Tyr_Helix310_gplus_gminus',
    'IC_Tyr_Helix310_gplus_gplus',
    'IC_Tyr_Helix310_gplus_trans',
    'IC_Tyr_Helix310_trans_gminus',
    'IC_Tyr_Helix310_trans_gplus',
    'IC_Tyr_Helix310_trans_trans',
    'IC_Tyr_HelixAlpha',
    'IC_Tyr_HelixAlpha_gminus_gminus',
    'IC_Tyr_HelixAlpha_gminus_gplus',
    'IC_Tyr_HelixAlpha_gminus_trans',
    'IC_Tyr_HelixAlpha_gplus_gminus',
    'IC_Tyr_HelixAlpha_gplus_gplus',
    'IC_Tyr_HelixAlpha_gplus_trans',
    'IC_Tyr_HelixAlpha_trans_gminus',
    'IC_Tyr_HelixAlpha_trans_gplus',
    'IC_Tyr_HelixAlpha_trans_trans',
    'IC_Tyr_HelixPPII',
    'IC_Tyr_HelixPPII_gminus_gminus',
    'IC_Tyr_HelixPPII_gminus_gplus',
    'IC_Tyr_HelixPPII_gminus_trans',
    'IC_Tyr_HelixPPII_gplus_gminus',
    'IC_Tyr_HelixPPII_gplus_gplus',
    'IC_Tyr_HelixPPII_gplus_trans',
    'IC_Tyr_HelixPPII_trans_gminus',
    'IC_Tyr_HelixPPII_trans_gplus',
    'IC_Tyr_HelixPPII_trans_trans',
    'IC_Tyr_HelixPi',
    'IC_Tyr_HelixPi_gminus_gminus',
    'IC_Tyr_HelixPi_gminus_gplus',
    'IC_Tyr_HelixPi_gminus_trans',
    'IC_Tyr_HelixPi_gplus_gminus',
    'IC_Tyr_HelixPi_gplus_gplus',
    'IC_Tyr_HelixPi_gplus_trans',
    'IC_Tyr_HelixPi_trans_gminus',
    'IC_Tyr_HelixPi_trans_gplus',
    'IC_Tyr_HelixPi_trans_trans',
    'IC_Tyr_Strand',
    'IC_Tyr_StrandAntiParallel',
    'IC_Tyr_StrandAntiParallel_gminus_gminus',
    'IC_Tyr_StrandAntiParallel_gminus_gplus',
    'IC_Tyr_StrandAntiParallel_gminus_trans',
    'IC_Tyr_StrandAntiParallel_gplus_gminus',
    'IC_Tyr_StrandAntiParallel_gplus_gplus',
    'IC_Tyr_StrandAntiParallel_gplus_trans',
    'IC_Tyr_StrandAntiParallel_trans_gminus',
    'IC_Tyr_StrandAntiParallel_trans_gplus',
    'IC_Tyr_StrandAntiParallel_trans_trans',
    'IC_Tyr_StrandParallel',
    'IC_Tyr_StrandParallel_gminus_gminus',
    'IC_Tyr_StrandParallel_gminus_gplus',
    'IC_Tyr_StrandParallel_gminus_trans',
    'IC_Tyr_StrandParallel_gplus_gminus',
    'IC_Tyr_StrandParallel_gplus_gplus',
    'IC_Tyr_StrandParallel_gplus_trans',
    'IC_Tyr_StrandParallel_trans_gminus',
    'IC_Tyr_StrandParallel_trans_gplus',
    'IC_Tyr_StrandParallel_trans_trans',
    'IC_Tyr_Strand_gminus_gminus',
    'IC_Tyr_Strand_gminus_gplus',
    'IC_Tyr_Strand_gminus_trans',
    'IC_Tyr_Strand_gplus_gminus',
    'IC_Tyr_Strand_gplus_gplus',
    'IC_Tyr_Strand_gplus_trans',
    'IC_Tyr_Strand_trans_gminus',
    'IC_Tyr_Strand_trans_gplus',
    'IC_Tyr_Strand_trans_trans',
    'IC_Tyr_Turn',
    'IC_Tyr_Turn_gminus_gminus',
    'IC_Tyr_Turn_gminus_gplus',
    'IC_Tyr_Turn_gminus_trans',
    'IC_Tyr_Turn_gplus_gminus',
    'IC_Tyr_Turn_gplus_gplus',
    'IC_Tyr_Turn_gplus_trans',
    'IC_Tyr_Turn_trans_gminus',
    'IC_Tyr_Turn_trans_gplus',
    'IC_Tyr_Turn_trans_trans',
]
