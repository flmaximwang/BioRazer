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

IC_Arg_HelixAlpha = build_template('ARG', 'alpha-helix', 'canonical')
IC_Arg_HelixAlpha_gminus_gminus = build_template('ARG', 'alpha-helix', 'g-/g-')
IC_Arg_HelixAlpha_gminus_trans = build_template('ARG', 'alpha-helix', 'g-/t')
IC_Arg_HelixAlpha_gminus_gplus = build_template('ARG', 'alpha-helix', 'g-/g+')
IC_Arg_HelixAlpha_trans_gminus = build_template('ARG', 'alpha-helix', 't/g-')
IC_Arg_HelixAlpha_trans_trans = build_template('ARG', 'alpha-helix', 't/t')
IC_Arg_HelixAlpha_trans_gplus = build_template('ARG', 'alpha-helix', 't/g+')
IC_Arg_HelixAlpha_gplus_gminus = build_template('ARG', 'alpha-helix', 'g+/g-')
IC_Arg_HelixAlpha_gplus_trans = build_template('ARG', 'alpha-helix', 'g+/t')
IC_Arg_HelixAlpha_gplus_gplus = build_template('ARG', 'alpha-helix', 'g+/g+')
IC_Arg_Helix310 = build_template('ARG', '3-10-helix', 'canonical')
IC_Arg_Helix310_gminus_gminus = build_template('ARG', '3-10-helix', 'g-/g-')
IC_Arg_Helix310_gminus_trans = build_template('ARG', '3-10-helix', 'g-/t')
IC_Arg_Helix310_gminus_gplus = build_template('ARG', '3-10-helix', 'g-/g+')
IC_Arg_Helix310_trans_gminus = build_template('ARG', '3-10-helix', 't/g-')
IC_Arg_Helix310_trans_trans = build_template('ARG', '3-10-helix', 't/t')
IC_Arg_Helix310_trans_gplus = build_template('ARG', '3-10-helix', 't/g+')
IC_Arg_Helix310_gplus_gminus = build_template('ARG', '3-10-helix', 'g+/g-')
IC_Arg_Helix310_gplus_trans = build_template('ARG', '3-10-helix', 'g+/t')
IC_Arg_Helix310_gplus_gplus = build_template('ARG', '3-10-helix', 'g+/g+')
IC_Arg_HelixPi = build_template('ARG', 'pi-helix', 'canonical')
IC_Arg_HelixPi_gminus_gminus = build_template('ARG', 'pi-helix', 'g-/g-')
IC_Arg_HelixPi_gminus_trans = build_template('ARG', 'pi-helix', 'g-/t')
IC_Arg_HelixPi_gminus_gplus = build_template('ARG', 'pi-helix', 'g-/g+')
IC_Arg_HelixPi_trans_gminus = build_template('ARG', 'pi-helix', 't/g-')
IC_Arg_HelixPi_trans_trans = build_template('ARG', 'pi-helix', 't/t')
IC_Arg_HelixPi_trans_gplus = build_template('ARG', 'pi-helix', 't/g+')
IC_Arg_HelixPi_gplus_gminus = build_template('ARG', 'pi-helix', 'g+/g-')
IC_Arg_HelixPi_gplus_trans = build_template('ARG', 'pi-helix', 'g+/t')
IC_Arg_HelixPi_gplus_gplus = build_template('ARG', 'pi-helix', 'g+/g+')
IC_Arg_HelixPPII = build_template('ARG', 'polyproline-II', 'canonical')
IC_Arg_HelixPPII_gminus_gminus = build_template('ARG', 'polyproline-II', 'g-/g-')
IC_Arg_HelixPPII_gminus_trans = build_template('ARG', 'polyproline-II', 'g-/t')
IC_Arg_HelixPPII_gminus_gplus = build_template('ARG', 'polyproline-II', 'g-/g+')
IC_Arg_HelixPPII_trans_gminus = build_template('ARG', 'polyproline-II', 't/g-')
IC_Arg_HelixPPII_trans_trans = build_template('ARG', 'polyproline-II', 't/t')
IC_Arg_HelixPPII_trans_gplus = build_template('ARG', 'polyproline-II', 't/g+')
IC_Arg_HelixPPII_gplus_gminus = build_template('ARG', 'polyproline-II', 'g+/g-')
IC_Arg_HelixPPII_gplus_trans = build_template('ARG', 'polyproline-II', 'g+/t')
IC_Arg_HelixPPII_gplus_gplus = build_template('ARG', 'polyproline-II', 'g+/g+')
IC_Arg_Strand = build_template('ARG', 'beta-strand', 'canonical')
IC_Arg_Strand_gminus_gminus = build_template('ARG', 'beta-strand', 'g-/g-')
IC_Arg_Strand_gminus_trans = build_template('ARG', 'beta-strand', 'g-/t')
IC_Arg_Strand_gminus_gplus = build_template('ARG', 'beta-strand', 'g-/g+')
IC_Arg_Strand_trans_gminus = build_template('ARG', 'beta-strand', 't/g-')
IC_Arg_Strand_trans_trans = build_template('ARG', 'beta-strand', 't/t')
IC_Arg_Strand_trans_gplus = build_template('ARG', 'beta-strand', 't/g+')
IC_Arg_Strand_gplus_gminus = build_template('ARG', 'beta-strand', 'g+/g-')
IC_Arg_Strand_gplus_trans = build_template('ARG', 'beta-strand', 'g+/t')
IC_Arg_Strand_gplus_gplus = build_template('ARG', 'beta-strand', 'g+/g+')
IC_Arg_StrandParallel = build_template('ARG', 'parallel-beta-strand', 'canonical')
IC_Arg_StrandParallel_gminus_gminus = build_template('ARG', 'parallel-beta-strand', 'g-/g-')
IC_Arg_StrandParallel_gminus_trans = build_template('ARG', 'parallel-beta-strand', 'g-/t')
IC_Arg_StrandParallel_gminus_gplus = build_template('ARG', 'parallel-beta-strand', 'g-/g+')
IC_Arg_StrandParallel_trans_gminus = build_template('ARG', 'parallel-beta-strand', 't/g-')
IC_Arg_StrandParallel_trans_trans = build_template('ARG', 'parallel-beta-strand', 't/t')
IC_Arg_StrandParallel_trans_gplus = build_template('ARG', 'parallel-beta-strand', 't/g+')
IC_Arg_StrandParallel_gplus_gminus = build_template('ARG', 'parallel-beta-strand', 'g+/g-')
IC_Arg_StrandParallel_gplus_trans = build_template('ARG', 'parallel-beta-strand', 'g+/t')
IC_Arg_StrandParallel_gplus_gplus = build_template('ARG', 'parallel-beta-strand', 'g+/g+')
IC_Arg_StrandAntiParallel = build_template('ARG', 'antiparallel-beta-strand', 'canonical')
IC_Arg_StrandAntiParallel_gminus_gminus = build_template('ARG', 'antiparallel-beta-strand', 'g-/g-')
IC_Arg_StrandAntiParallel_gminus_trans = build_template('ARG', 'antiparallel-beta-strand', 'g-/t')
IC_Arg_StrandAntiParallel_gminus_gplus = build_template('ARG', 'antiparallel-beta-strand', 'g-/g+')
IC_Arg_StrandAntiParallel_trans_gminus = build_template('ARG', 'antiparallel-beta-strand', 't/g-')
IC_Arg_StrandAntiParallel_trans_trans = build_template('ARG', 'antiparallel-beta-strand', 't/t')
IC_Arg_StrandAntiParallel_trans_gplus = build_template('ARG', 'antiparallel-beta-strand', 't/g+')
IC_Arg_StrandAntiParallel_gplus_gminus = build_template('ARG', 'antiparallel-beta-strand', 'g+/g-')
IC_Arg_StrandAntiParallel_gplus_trans = build_template('ARG', 'antiparallel-beta-strand', 'g+/t')
IC_Arg_StrandAntiParallel_gplus_gplus = build_template('ARG', 'antiparallel-beta-strand', 'g+/g+')
IC_Arg_Bridge = build_template('ARG', 'beta-bridge', 'canonical')
IC_Arg_Bridge_gminus_gminus = build_template('ARG', 'beta-bridge', 'g-/g-')
IC_Arg_Bridge_gminus_trans = build_template('ARG', 'beta-bridge', 'g-/t')
IC_Arg_Bridge_gminus_gplus = build_template('ARG', 'beta-bridge', 'g-/g+')
IC_Arg_Bridge_trans_gminus = build_template('ARG', 'beta-bridge', 't/g-')
IC_Arg_Bridge_trans_trans = build_template('ARG', 'beta-bridge', 't/t')
IC_Arg_Bridge_trans_gplus = build_template('ARG', 'beta-bridge', 't/g+')
IC_Arg_Bridge_gplus_gminus = build_template('ARG', 'beta-bridge', 'g+/g-')
IC_Arg_Bridge_gplus_trans = build_template('ARG', 'beta-bridge', 'g+/t')
IC_Arg_Bridge_gplus_gplus = build_template('ARG', 'beta-bridge', 'g+/g+')
IC_Arg_Turn = build_template('ARG', 'turn', 'canonical')
IC_Arg_Turn_gminus_gminus = build_template('ARG', 'turn', 'g-/g-')
IC_Arg_Turn_gminus_trans = build_template('ARG', 'turn', 'g-/t')
IC_Arg_Turn_gminus_gplus = build_template('ARG', 'turn', 'g-/g+')
IC_Arg_Turn_trans_gminus = build_template('ARG', 'turn', 't/g-')
IC_Arg_Turn_trans_trans = build_template('ARG', 'turn', 't/t')
IC_Arg_Turn_trans_gplus = build_template('ARG', 'turn', 't/g+')
IC_Arg_Turn_gplus_gminus = build_template('ARG', 'turn', 'g+/g-')
IC_Arg_Turn_gplus_trans = build_template('ARG', 'turn', 'g+/t')
IC_Arg_Turn_gplus_gplus = build_template('ARG', 'turn', 'g+/g+')
IC_Arg_Bend = build_template('ARG', 'bend', 'canonical')
IC_Arg_Bend_gminus_gminus = build_template('ARG', 'bend', 'g-/g-')
IC_Arg_Bend_gminus_trans = build_template('ARG', 'bend', 'g-/t')
IC_Arg_Bend_gminus_gplus = build_template('ARG', 'bend', 'g-/g+')
IC_Arg_Bend_trans_gminus = build_template('ARG', 'bend', 't/g-')
IC_Arg_Bend_trans_trans = build_template('ARG', 'bend', 't/t')
IC_Arg_Bend_trans_gplus = build_template('ARG', 'bend', 't/g+')
IC_Arg_Bend_gplus_gminus = build_template('ARG', 'bend', 'g+/g-')
IC_Arg_Bend_gplus_trans = build_template('ARG', 'bend', 'g+/t')
IC_Arg_Bend_gplus_gplus = build_template('ARG', 'bend', 'g+/g+')
IC_Arg_Coil = build_template('ARG', 'coil', 'canonical')
IC_Arg_Coil_gminus_gminus = build_template('ARG', 'coil', 'g-/g-')
IC_Arg_Coil_gminus_trans = build_template('ARG', 'coil', 'g-/t')
IC_Arg_Coil_gminus_gplus = build_template('ARG', 'coil', 'g-/g+')
IC_Arg_Coil_trans_gminus = build_template('ARG', 'coil', 't/g-')
IC_Arg_Coil_trans_trans = build_template('ARG', 'coil', 't/t')
IC_Arg_Coil_trans_gplus = build_template('ARG', 'coil', 't/g+')
IC_Arg_Coil_gplus_gminus = build_template('ARG', 'coil', 'g+/g-')
IC_Arg_Coil_gplus_trans = build_template('ARG', 'coil', 'g+/t')
IC_Arg_Coil_gplus_gplus = build_template('ARG', 'coil', 'g+/g+')
IC_Arg_CisPeptide = build_template('ARG', 'cis-peptide-bond', 'canonical')
IC_Arg_CisPeptide_gminus_gminus = build_template('ARG', 'cis-peptide-bond', 'g-/g-')
IC_Arg_CisPeptide_gminus_trans = build_template('ARG', 'cis-peptide-bond', 'g-/t')
IC_Arg_CisPeptide_gminus_gplus = build_template('ARG', 'cis-peptide-bond', 'g-/g+')
IC_Arg_CisPeptide_trans_gminus = build_template('ARG', 'cis-peptide-bond', 't/g-')
IC_Arg_CisPeptide_trans_trans = build_template('ARG', 'cis-peptide-bond', 't/t')
IC_Arg_CisPeptide_trans_gplus = build_template('ARG', 'cis-peptide-bond', 't/g+')
IC_Arg_CisPeptide_gplus_gminus = build_template('ARG', 'cis-peptide-bond', 'g+/g-')
IC_Arg_CisPeptide_gplus_trans = build_template('ARG', 'cis-peptide-bond', 'g+/t')
IC_Arg_CisPeptide_gplus_gplus = build_template('ARG', 'cis-peptide-bond', 'g+/g+')

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
