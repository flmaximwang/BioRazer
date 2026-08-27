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

IC_Leu_HelixAlpha = build_template('LEU', 'alpha-helix', 'canonical')
IC_Leu_HelixAlpha_gminus_gminus = build_template('LEU', 'alpha-helix', 'g-/g-')
IC_Leu_HelixAlpha_gminus_trans = build_template('LEU', 'alpha-helix', 'g-/t')
IC_Leu_HelixAlpha_gminus_gplus = build_template('LEU', 'alpha-helix', 'g-/g+')
IC_Leu_HelixAlpha_trans_gminus = build_template('LEU', 'alpha-helix', 't/g-')
IC_Leu_HelixAlpha_trans_trans = build_template('LEU', 'alpha-helix', 't/t')
IC_Leu_HelixAlpha_trans_gplus = build_template('LEU', 'alpha-helix', 't/g+')
IC_Leu_HelixAlpha_gplus_gminus = build_template('LEU', 'alpha-helix', 'g+/g-')
IC_Leu_HelixAlpha_gplus_trans = build_template('LEU', 'alpha-helix', 'g+/t')
IC_Leu_HelixAlpha_gplus_gplus = build_template('LEU', 'alpha-helix', 'g+/g+')
IC_Leu_Helix310 = build_template('LEU', '3-10-helix', 'canonical')
IC_Leu_Helix310_gminus_gminus = build_template('LEU', '3-10-helix', 'g-/g-')
IC_Leu_Helix310_gminus_trans = build_template('LEU', '3-10-helix', 'g-/t')
IC_Leu_Helix310_gminus_gplus = build_template('LEU', '3-10-helix', 'g-/g+')
IC_Leu_Helix310_trans_gminus = build_template('LEU', '3-10-helix', 't/g-')
IC_Leu_Helix310_trans_trans = build_template('LEU', '3-10-helix', 't/t')
IC_Leu_Helix310_trans_gplus = build_template('LEU', '3-10-helix', 't/g+')
IC_Leu_Helix310_gplus_gminus = build_template('LEU', '3-10-helix', 'g+/g-')
IC_Leu_Helix310_gplus_trans = build_template('LEU', '3-10-helix', 'g+/t')
IC_Leu_Helix310_gplus_gplus = build_template('LEU', '3-10-helix', 'g+/g+')
IC_Leu_HelixPi = build_template('LEU', 'pi-helix', 'canonical')
IC_Leu_HelixPi_gminus_gminus = build_template('LEU', 'pi-helix', 'g-/g-')
IC_Leu_HelixPi_gminus_trans = build_template('LEU', 'pi-helix', 'g-/t')
IC_Leu_HelixPi_gminus_gplus = build_template('LEU', 'pi-helix', 'g-/g+')
IC_Leu_HelixPi_trans_gminus = build_template('LEU', 'pi-helix', 't/g-')
IC_Leu_HelixPi_trans_trans = build_template('LEU', 'pi-helix', 't/t')
IC_Leu_HelixPi_trans_gplus = build_template('LEU', 'pi-helix', 't/g+')
IC_Leu_HelixPi_gplus_gminus = build_template('LEU', 'pi-helix', 'g+/g-')
IC_Leu_HelixPi_gplus_trans = build_template('LEU', 'pi-helix', 'g+/t')
IC_Leu_HelixPi_gplus_gplus = build_template('LEU', 'pi-helix', 'g+/g+')
IC_Leu_HelixPPII = build_template('LEU', 'polyproline-II', 'canonical')
IC_Leu_HelixPPII_gminus_gminus = build_template('LEU', 'polyproline-II', 'g-/g-')
IC_Leu_HelixPPII_gminus_trans = build_template('LEU', 'polyproline-II', 'g-/t')
IC_Leu_HelixPPII_gminus_gplus = build_template('LEU', 'polyproline-II', 'g-/g+')
IC_Leu_HelixPPII_trans_gminus = build_template('LEU', 'polyproline-II', 't/g-')
IC_Leu_HelixPPII_trans_trans = build_template('LEU', 'polyproline-II', 't/t')
IC_Leu_HelixPPII_trans_gplus = build_template('LEU', 'polyproline-II', 't/g+')
IC_Leu_HelixPPII_gplus_gminus = build_template('LEU', 'polyproline-II', 'g+/g-')
IC_Leu_HelixPPII_gplus_trans = build_template('LEU', 'polyproline-II', 'g+/t')
IC_Leu_HelixPPII_gplus_gplus = build_template('LEU', 'polyproline-II', 'g+/g+')
IC_Leu_Strand = build_template('LEU', 'beta-strand', 'canonical')
IC_Leu_Strand_gminus_gminus = build_template('LEU', 'beta-strand', 'g-/g-')
IC_Leu_Strand_gminus_trans = build_template('LEU', 'beta-strand', 'g-/t')
IC_Leu_Strand_gminus_gplus = build_template('LEU', 'beta-strand', 'g-/g+')
IC_Leu_Strand_trans_gminus = build_template('LEU', 'beta-strand', 't/g-')
IC_Leu_Strand_trans_trans = build_template('LEU', 'beta-strand', 't/t')
IC_Leu_Strand_trans_gplus = build_template('LEU', 'beta-strand', 't/g+')
IC_Leu_Strand_gplus_gminus = build_template('LEU', 'beta-strand', 'g+/g-')
IC_Leu_Strand_gplus_trans = build_template('LEU', 'beta-strand', 'g+/t')
IC_Leu_Strand_gplus_gplus = build_template('LEU', 'beta-strand', 'g+/g+')
IC_Leu_StrandParallel = build_template('LEU', 'parallel-beta-strand', 'canonical')
IC_Leu_StrandParallel_gminus_gminus = build_template('LEU', 'parallel-beta-strand', 'g-/g-')
IC_Leu_StrandParallel_gminus_trans = build_template('LEU', 'parallel-beta-strand', 'g-/t')
IC_Leu_StrandParallel_gminus_gplus = build_template('LEU', 'parallel-beta-strand', 'g-/g+')
IC_Leu_StrandParallel_trans_gminus = build_template('LEU', 'parallel-beta-strand', 't/g-')
IC_Leu_StrandParallel_trans_trans = build_template('LEU', 'parallel-beta-strand', 't/t')
IC_Leu_StrandParallel_trans_gplus = build_template('LEU', 'parallel-beta-strand', 't/g+')
IC_Leu_StrandParallel_gplus_gminus = build_template('LEU', 'parallel-beta-strand', 'g+/g-')
IC_Leu_StrandParallel_gplus_trans = build_template('LEU', 'parallel-beta-strand', 'g+/t')
IC_Leu_StrandParallel_gplus_gplus = build_template('LEU', 'parallel-beta-strand', 'g+/g+')
IC_Leu_StrandAntiParallel = build_template('LEU', 'antiparallel-beta-strand', 'canonical')
IC_Leu_StrandAntiParallel_gminus_gminus = build_template('LEU', 'antiparallel-beta-strand', 'g-/g-')
IC_Leu_StrandAntiParallel_gminus_trans = build_template('LEU', 'antiparallel-beta-strand', 'g-/t')
IC_Leu_StrandAntiParallel_gminus_gplus = build_template('LEU', 'antiparallel-beta-strand', 'g-/g+')
IC_Leu_StrandAntiParallel_trans_gminus = build_template('LEU', 'antiparallel-beta-strand', 't/g-')
IC_Leu_StrandAntiParallel_trans_trans = build_template('LEU', 'antiparallel-beta-strand', 't/t')
IC_Leu_StrandAntiParallel_trans_gplus = build_template('LEU', 'antiparallel-beta-strand', 't/g+')
IC_Leu_StrandAntiParallel_gplus_gminus = build_template('LEU', 'antiparallel-beta-strand', 'g+/g-')
IC_Leu_StrandAntiParallel_gplus_trans = build_template('LEU', 'antiparallel-beta-strand', 'g+/t')
IC_Leu_StrandAntiParallel_gplus_gplus = build_template('LEU', 'antiparallel-beta-strand', 'g+/g+')
IC_Leu_Bridge = build_template('LEU', 'beta-bridge', 'canonical')
IC_Leu_Bridge_gminus_gminus = build_template('LEU', 'beta-bridge', 'g-/g-')
IC_Leu_Bridge_gminus_trans = build_template('LEU', 'beta-bridge', 'g-/t')
IC_Leu_Bridge_gminus_gplus = build_template('LEU', 'beta-bridge', 'g-/g+')
IC_Leu_Bridge_trans_gminus = build_template('LEU', 'beta-bridge', 't/g-')
IC_Leu_Bridge_trans_trans = build_template('LEU', 'beta-bridge', 't/t')
IC_Leu_Bridge_trans_gplus = build_template('LEU', 'beta-bridge', 't/g+')
IC_Leu_Bridge_gplus_gminus = build_template('LEU', 'beta-bridge', 'g+/g-')
IC_Leu_Bridge_gplus_trans = build_template('LEU', 'beta-bridge', 'g+/t')
IC_Leu_Bridge_gplus_gplus = build_template('LEU', 'beta-bridge', 'g+/g+')
IC_Leu_Turn = build_template('LEU', 'turn', 'canonical')
IC_Leu_Turn_gminus_gminus = build_template('LEU', 'turn', 'g-/g-')
IC_Leu_Turn_gminus_trans = build_template('LEU', 'turn', 'g-/t')
IC_Leu_Turn_gminus_gplus = build_template('LEU', 'turn', 'g-/g+')
IC_Leu_Turn_trans_gminus = build_template('LEU', 'turn', 't/g-')
IC_Leu_Turn_trans_trans = build_template('LEU', 'turn', 't/t')
IC_Leu_Turn_trans_gplus = build_template('LEU', 'turn', 't/g+')
IC_Leu_Turn_gplus_gminus = build_template('LEU', 'turn', 'g+/g-')
IC_Leu_Turn_gplus_trans = build_template('LEU', 'turn', 'g+/t')
IC_Leu_Turn_gplus_gplus = build_template('LEU', 'turn', 'g+/g+')
IC_Leu_Bend = build_template('LEU', 'bend', 'canonical')
IC_Leu_Bend_gminus_gminus = build_template('LEU', 'bend', 'g-/g-')
IC_Leu_Bend_gminus_trans = build_template('LEU', 'bend', 'g-/t')
IC_Leu_Bend_gminus_gplus = build_template('LEU', 'bend', 'g-/g+')
IC_Leu_Bend_trans_gminus = build_template('LEU', 'bend', 't/g-')
IC_Leu_Bend_trans_trans = build_template('LEU', 'bend', 't/t')
IC_Leu_Bend_trans_gplus = build_template('LEU', 'bend', 't/g+')
IC_Leu_Bend_gplus_gminus = build_template('LEU', 'bend', 'g+/g-')
IC_Leu_Bend_gplus_trans = build_template('LEU', 'bend', 'g+/t')
IC_Leu_Bend_gplus_gplus = build_template('LEU', 'bend', 'g+/g+')
IC_Leu_Coil = build_template('LEU', 'coil', 'canonical')
IC_Leu_Coil_gminus_gminus = build_template('LEU', 'coil', 'g-/g-')
IC_Leu_Coil_gminus_trans = build_template('LEU', 'coil', 'g-/t')
IC_Leu_Coil_gminus_gplus = build_template('LEU', 'coil', 'g-/g+')
IC_Leu_Coil_trans_gminus = build_template('LEU', 'coil', 't/g-')
IC_Leu_Coil_trans_trans = build_template('LEU', 'coil', 't/t')
IC_Leu_Coil_trans_gplus = build_template('LEU', 'coil', 't/g+')
IC_Leu_Coil_gplus_gminus = build_template('LEU', 'coil', 'g+/g-')
IC_Leu_Coil_gplus_trans = build_template('LEU', 'coil', 'g+/t')
IC_Leu_Coil_gplus_gplus = build_template('LEU', 'coil', 'g+/g+')
IC_Leu_CisPeptide = build_template('LEU', 'cis-peptide-bond', 'canonical')
IC_Leu_CisPeptide_gminus_gminus = build_template('LEU', 'cis-peptide-bond', 'g-/g-')
IC_Leu_CisPeptide_gminus_trans = build_template('LEU', 'cis-peptide-bond', 'g-/t')
IC_Leu_CisPeptide_gminus_gplus = build_template('LEU', 'cis-peptide-bond', 'g-/g+')
IC_Leu_CisPeptide_trans_gminus = build_template('LEU', 'cis-peptide-bond', 't/g-')
IC_Leu_CisPeptide_trans_trans = build_template('LEU', 'cis-peptide-bond', 't/t')
IC_Leu_CisPeptide_trans_gplus = build_template('LEU', 'cis-peptide-bond', 't/g+')
IC_Leu_CisPeptide_gplus_gminus = build_template('LEU', 'cis-peptide-bond', 'g+/g-')
IC_Leu_CisPeptide_gplus_trans = build_template('LEU', 'cis-peptide-bond', 'g+/t')
IC_Leu_CisPeptide_gplus_gplus = build_template('LEU', 'cis-peptide-bond', 'g+/g+')

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
