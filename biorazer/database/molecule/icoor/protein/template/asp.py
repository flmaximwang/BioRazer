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

IC_Asp_HelixAlpha = build_template('ASP', 'alpha-helix', 'canonical')
IC_Asp_HelixAlpha_gminus_gminus = build_template('ASP', 'alpha-helix', 'g-/g-')
IC_Asp_HelixAlpha_gminus_trans = build_template('ASP', 'alpha-helix', 'g-/t')
IC_Asp_HelixAlpha_gminus_gplus = build_template('ASP', 'alpha-helix', 'g-/g+')
IC_Asp_HelixAlpha_trans_gminus = build_template('ASP', 'alpha-helix', 't/g-')
IC_Asp_HelixAlpha_trans_trans = build_template('ASP', 'alpha-helix', 't/t')
IC_Asp_HelixAlpha_trans_gplus = build_template('ASP', 'alpha-helix', 't/g+')
IC_Asp_HelixAlpha_gplus_gminus = build_template('ASP', 'alpha-helix', 'g+/g-')
IC_Asp_HelixAlpha_gplus_trans = build_template('ASP', 'alpha-helix', 'g+/t')
IC_Asp_HelixAlpha_gplus_gplus = build_template('ASP', 'alpha-helix', 'g+/g+')
IC_Asp_Helix310 = build_template('ASP', '3-10-helix', 'canonical')
IC_Asp_Helix310_gminus_gminus = build_template('ASP', '3-10-helix', 'g-/g-')
IC_Asp_Helix310_gminus_trans = build_template('ASP', '3-10-helix', 'g-/t')
IC_Asp_Helix310_gminus_gplus = build_template('ASP', '3-10-helix', 'g-/g+')
IC_Asp_Helix310_trans_gminus = build_template('ASP', '3-10-helix', 't/g-')
IC_Asp_Helix310_trans_trans = build_template('ASP', '3-10-helix', 't/t')
IC_Asp_Helix310_trans_gplus = build_template('ASP', '3-10-helix', 't/g+')
IC_Asp_Helix310_gplus_gminus = build_template('ASP', '3-10-helix', 'g+/g-')
IC_Asp_Helix310_gplus_trans = build_template('ASP', '3-10-helix', 'g+/t')
IC_Asp_Helix310_gplus_gplus = build_template('ASP', '3-10-helix', 'g+/g+')
IC_Asp_HelixPi = build_template('ASP', 'pi-helix', 'canonical')
IC_Asp_HelixPi_gminus_gminus = build_template('ASP', 'pi-helix', 'g-/g-')
IC_Asp_HelixPi_gminus_trans = build_template('ASP', 'pi-helix', 'g-/t')
IC_Asp_HelixPi_gminus_gplus = build_template('ASP', 'pi-helix', 'g-/g+')
IC_Asp_HelixPi_trans_gminus = build_template('ASP', 'pi-helix', 't/g-')
IC_Asp_HelixPi_trans_trans = build_template('ASP', 'pi-helix', 't/t')
IC_Asp_HelixPi_trans_gplus = build_template('ASP', 'pi-helix', 't/g+')
IC_Asp_HelixPi_gplus_gminus = build_template('ASP', 'pi-helix', 'g+/g-')
IC_Asp_HelixPi_gplus_trans = build_template('ASP', 'pi-helix', 'g+/t')
IC_Asp_HelixPi_gplus_gplus = build_template('ASP', 'pi-helix', 'g+/g+')
IC_Asp_HelixPPII = build_template('ASP', 'polyproline-II', 'canonical')
IC_Asp_HelixPPII_gminus_gminus = build_template('ASP', 'polyproline-II', 'g-/g-')
IC_Asp_HelixPPII_gminus_trans = build_template('ASP', 'polyproline-II', 'g-/t')
IC_Asp_HelixPPII_gminus_gplus = build_template('ASP', 'polyproline-II', 'g-/g+')
IC_Asp_HelixPPII_trans_gminus = build_template('ASP', 'polyproline-II', 't/g-')
IC_Asp_HelixPPII_trans_trans = build_template('ASP', 'polyproline-II', 't/t')
IC_Asp_HelixPPII_trans_gplus = build_template('ASP', 'polyproline-II', 't/g+')
IC_Asp_HelixPPII_gplus_gminus = build_template('ASP', 'polyproline-II', 'g+/g-')
IC_Asp_HelixPPII_gplus_trans = build_template('ASP', 'polyproline-II', 'g+/t')
IC_Asp_HelixPPII_gplus_gplus = build_template('ASP', 'polyproline-II', 'g+/g+')
IC_Asp_Strand = build_template('ASP', 'beta-strand', 'canonical')
IC_Asp_Strand_gminus_gminus = build_template('ASP', 'beta-strand', 'g-/g-')
IC_Asp_Strand_gminus_trans = build_template('ASP', 'beta-strand', 'g-/t')
IC_Asp_Strand_gminus_gplus = build_template('ASP', 'beta-strand', 'g-/g+')
IC_Asp_Strand_trans_gminus = build_template('ASP', 'beta-strand', 't/g-')
IC_Asp_Strand_trans_trans = build_template('ASP', 'beta-strand', 't/t')
IC_Asp_Strand_trans_gplus = build_template('ASP', 'beta-strand', 't/g+')
IC_Asp_Strand_gplus_gminus = build_template('ASP', 'beta-strand', 'g+/g-')
IC_Asp_Strand_gplus_trans = build_template('ASP', 'beta-strand', 'g+/t')
IC_Asp_Strand_gplus_gplus = build_template('ASP', 'beta-strand', 'g+/g+')
IC_Asp_StrandParallel = build_template('ASP', 'parallel-beta-strand', 'canonical')
IC_Asp_StrandParallel_gminus_gminus = build_template('ASP', 'parallel-beta-strand', 'g-/g-')
IC_Asp_StrandParallel_gminus_trans = build_template('ASP', 'parallel-beta-strand', 'g-/t')
IC_Asp_StrandParallel_gminus_gplus = build_template('ASP', 'parallel-beta-strand', 'g-/g+')
IC_Asp_StrandParallel_trans_gminus = build_template('ASP', 'parallel-beta-strand', 't/g-')
IC_Asp_StrandParallel_trans_trans = build_template('ASP', 'parallel-beta-strand', 't/t')
IC_Asp_StrandParallel_trans_gplus = build_template('ASP', 'parallel-beta-strand', 't/g+')
IC_Asp_StrandParallel_gplus_gminus = build_template('ASP', 'parallel-beta-strand', 'g+/g-')
IC_Asp_StrandParallel_gplus_trans = build_template('ASP', 'parallel-beta-strand', 'g+/t')
IC_Asp_StrandParallel_gplus_gplus = build_template('ASP', 'parallel-beta-strand', 'g+/g+')
IC_Asp_StrandAntiParallel = build_template('ASP', 'antiparallel-beta-strand', 'canonical')
IC_Asp_StrandAntiParallel_gminus_gminus = build_template('ASP', 'antiparallel-beta-strand', 'g-/g-')
IC_Asp_StrandAntiParallel_gminus_trans = build_template('ASP', 'antiparallel-beta-strand', 'g-/t')
IC_Asp_StrandAntiParallel_gminus_gplus = build_template('ASP', 'antiparallel-beta-strand', 'g-/g+')
IC_Asp_StrandAntiParallel_trans_gminus = build_template('ASP', 'antiparallel-beta-strand', 't/g-')
IC_Asp_StrandAntiParallel_trans_trans = build_template('ASP', 'antiparallel-beta-strand', 't/t')
IC_Asp_StrandAntiParallel_trans_gplus = build_template('ASP', 'antiparallel-beta-strand', 't/g+')
IC_Asp_StrandAntiParallel_gplus_gminus = build_template('ASP', 'antiparallel-beta-strand', 'g+/g-')
IC_Asp_StrandAntiParallel_gplus_trans = build_template('ASP', 'antiparallel-beta-strand', 'g+/t')
IC_Asp_StrandAntiParallel_gplus_gplus = build_template('ASP', 'antiparallel-beta-strand', 'g+/g+')
IC_Asp_Bridge = build_template('ASP', 'beta-bridge', 'canonical')
IC_Asp_Bridge_gminus_gminus = build_template('ASP', 'beta-bridge', 'g-/g-')
IC_Asp_Bridge_gminus_trans = build_template('ASP', 'beta-bridge', 'g-/t')
IC_Asp_Bridge_gminus_gplus = build_template('ASP', 'beta-bridge', 'g-/g+')
IC_Asp_Bridge_trans_gminus = build_template('ASP', 'beta-bridge', 't/g-')
IC_Asp_Bridge_trans_trans = build_template('ASP', 'beta-bridge', 't/t')
IC_Asp_Bridge_trans_gplus = build_template('ASP', 'beta-bridge', 't/g+')
IC_Asp_Bridge_gplus_gminus = build_template('ASP', 'beta-bridge', 'g+/g-')
IC_Asp_Bridge_gplus_trans = build_template('ASP', 'beta-bridge', 'g+/t')
IC_Asp_Bridge_gplus_gplus = build_template('ASP', 'beta-bridge', 'g+/g+')
IC_Asp_Turn = build_template('ASP', 'turn', 'canonical')
IC_Asp_Turn_gminus_gminus = build_template('ASP', 'turn', 'g-/g-')
IC_Asp_Turn_gminus_trans = build_template('ASP', 'turn', 'g-/t')
IC_Asp_Turn_gminus_gplus = build_template('ASP', 'turn', 'g-/g+')
IC_Asp_Turn_trans_gminus = build_template('ASP', 'turn', 't/g-')
IC_Asp_Turn_trans_trans = build_template('ASP', 'turn', 't/t')
IC_Asp_Turn_trans_gplus = build_template('ASP', 'turn', 't/g+')
IC_Asp_Turn_gplus_gminus = build_template('ASP', 'turn', 'g+/g-')
IC_Asp_Turn_gplus_trans = build_template('ASP', 'turn', 'g+/t')
IC_Asp_Turn_gplus_gplus = build_template('ASP', 'turn', 'g+/g+')
IC_Asp_Bend = build_template('ASP', 'bend', 'canonical')
IC_Asp_Bend_gminus_gminus = build_template('ASP', 'bend', 'g-/g-')
IC_Asp_Bend_gminus_trans = build_template('ASP', 'bend', 'g-/t')
IC_Asp_Bend_gminus_gplus = build_template('ASP', 'bend', 'g-/g+')
IC_Asp_Bend_trans_gminus = build_template('ASP', 'bend', 't/g-')
IC_Asp_Bend_trans_trans = build_template('ASP', 'bend', 't/t')
IC_Asp_Bend_trans_gplus = build_template('ASP', 'bend', 't/g+')
IC_Asp_Bend_gplus_gminus = build_template('ASP', 'bend', 'g+/g-')
IC_Asp_Bend_gplus_trans = build_template('ASP', 'bend', 'g+/t')
IC_Asp_Bend_gplus_gplus = build_template('ASP', 'bend', 'g+/g+')
IC_Asp_Coil = build_template('ASP', 'coil', 'canonical')
IC_Asp_Coil_gminus_gminus = build_template('ASP', 'coil', 'g-/g-')
IC_Asp_Coil_gminus_trans = build_template('ASP', 'coil', 'g-/t')
IC_Asp_Coil_gminus_gplus = build_template('ASP', 'coil', 'g-/g+')
IC_Asp_Coil_trans_gminus = build_template('ASP', 'coil', 't/g-')
IC_Asp_Coil_trans_trans = build_template('ASP', 'coil', 't/t')
IC_Asp_Coil_trans_gplus = build_template('ASP', 'coil', 't/g+')
IC_Asp_Coil_gplus_gminus = build_template('ASP', 'coil', 'g+/g-')
IC_Asp_Coil_gplus_trans = build_template('ASP', 'coil', 'g+/t')
IC_Asp_Coil_gplus_gplus = build_template('ASP', 'coil', 'g+/g+')
IC_Asp_CisPeptide = build_template('ASP', 'cis-peptide-bond', 'canonical')
IC_Asp_CisPeptide_gminus_gminus = build_template('ASP', 'cis-peptide-bond', 'g-/g-')
IC_Asp_CisPeptide_gminus_trans = build_template('ASP', 'cis-peptide-bond', 'g-/t')
IC_Asp_CisPeptide_gminus_gplus = build_template('ASP', 'cis-peptide-bond', 'g-/g+')
IC_Asp_CisPeptide_trans_gminus = build_template('ASP', 'cis-peptide-bond', 't/g-')
IC_Asp_CisPeptide_trans_trans = build_template('ASP', 'cis-peptide-bond', 't/t')
IC_Asp_CisPeptide_trans_gplus = build_template('ASP', 'cis-peptide-bond', 't/g+')
IC_Asp_CisPeptide_gplus_gminus = build_template('ASP', 'cis-peptide-bond', 'g+/g-')
IC_Asp_CisPeptide_gplus_trans = build_template('ASP', 'cis-peptide-bond', 'g+/t')
IC_Asp_CisPeptide_gplus_gplus = build_template('ASP', 'cis-peptide-bond', 'g+/g+')

__all__ = [
    'IC_Asp_Bend',
    'IC_Asp_Bend_gminus_gminus',
    'IC_Asp_Bend_gminus_gplus',
    'IC_Asp_Bend_gminus_trans',
    'IC_Asp_Bend_gplus_gminus',
    'IC_Asp_Bend_gplus_gplus',
    'IC_Asp_Bend_gplus_trans',
    'IC_Asp_Bend_trans_gminus',
    'IC_Asp_Bend_trans_gplus',
    'IC_Asp_Bend_trans_trans',
    'IC_Asp_Bridge',
    'IC_Asp_Bridge_gminus_gminus',
    'IC_Asp_Bridge_gminus_gplus',
    'IC_Asp_Bridge_gminus_trans',
    'IC_Asp_Bridge_gplus_gminus',
    'IC_Asp_Bridge_gplus_gplus',
    'IC_Asp_Bridge_gplus_trans',
    'IC_Asp_Bridge_trans_gminus',
    'IC_Asp_Bridge_trans_gplus',
    'IC_Asp_Bridge_trans_trans',
    'IC_Asp_CisPeptide',
    'IC_Asp_CisPeptide_gminus_gminus',
    'IC_Asp_CisPeptide_gminus_gplus',
    'IC_Asp_CisPeptide_gminus_trans',
    'IC_Asp_CisPeptide_gplus_gminus',
    'IC_Asp_CisPeptide_gplus_gplus',
    'IC_Asp_CisPeptide_gplus_trans',
    'IC_Asp_CisPeptide_trans_gminus',
    'IC_Asp_CisPeptide_trans_gplus',
    'IC_Asp_CisPeptide_trans_trans',
    'IC_Asp_Coil',
    'IC_Asp_Coil_gminus_gminus',
    'IC_Asp_Coil_gminus_gplus',
    'IC_Asp_Coil_gminus_trans',
    'IC_Asp_Coil_gplus_gminus',
    'IC_Asp_Coil_gplus_gplus',
    'IC_Asp_Coil_gplus_trans',
    'IC_Asp_Coil_trans_gminus',
    'IC_Asp_Coil_trans_gplus',
    'IC_Asp_Coil_trans_trans',
    'IC_Asp_Helix310',
    'IC_Asp_Helix310_gminus_gminus',
    'IC_Asp_Helix310_gminus_gplus',
    'IC_Asp_Helix310_gminus_trans',
    'IC_Asp_Helix310_gplus_gminus',
    'IC_Asp_Helix310_gplus_gplus',
    'IC_Asp_Helix310_gplus_trans',
    'IC_Asp_Helix310_trans_gminus',
    'IC_Asp_Helix310_trans_gplus',
    'IC_Asp_Helix310_trans_trans',
    'IC_Asp_HelixAlpha',
    'IC_Asp_HelixAlpha_gminus_gminus',
    'IC_Asp_HelixAlpha_gminus_gplus',
    'IC_Asp_HelixAlpha_gminus_trans',
    'IC_Asp_HelixAlpha_gplus_gminus',
    'IC_Asp_HelixAlpha_gplus_gplus',
    'IC_Asp_HelixAlpha_gplus_trans',
    'IC_Asp_HelixAlpha_trans_gminus',
    'IC_Asp_HelixAlpha_trans_gplus',
    'IC_Asp_HelixAlpha_trans_trans',
    'IC_Asp_HelixPPII',
    'IC_Asp_HelixPPII_gminus_gminus',
    'IC_Asp_HelixPPII_gminus_gplus',
    'IC_Asp_HelixPPII_gminus_trans',
    'IC_Asp_HelixPPII_gplus_gminus',
    'IC_Asp_HelixPPII_gplus_gplus',
    'IC_Asp_HelixPPII_gplus_trans',
    'IC_Asp_HelixPPII_trans_gminus',
    'IC_Asp_HelixPPII_trans_gplus',
    'IC_Asp_HelixPPII_trans_trans',
    'IC_Asp_HelixPi',
    'IC_Asp_HelixPi_gminus_gminus',
    'IC_Asp_HelixPi_gminus_gplus',
    'IC_Asp_HelixPi_gminus_trans',
    'IC_Asp_HelixPi_gplus_gminus',
    'IC_Asp_HelixPi_gplus_gplus',
    'IC_Asp_HelixPi_gplus_trans',
    'IC_Asp_HelixPi_trans_gminus',
    'IC_Asp_HelixPi_trans_gplus',
    'IC_Asp_HelixPi_trans_trans',
    'IC_Asp_Strand',
    'IC_Asp_StrandAntiParallel',
    'IC_Asp_StrandAntiParallel_gminus_gminus',
    'IC_Asp_StrandAntiParallel_gminus_gplus',
    'IC_Asp_StrandAntiParallel_gminus_trans',
    'IC_Asp_StrandAntiParallel_gplus_gminus',
    'IC_Asp_StrandAntiParallel_gplus_gplus',
    'IC_Asp_StrandAntiParallel_gplus_trans',
    'IC_Asp_StrandAntiParallel_trans_gminus',
    'IC_Asp_StrandAntiParallel_trans_gplus',
    'IC_Asp_StrandAntiParallel_trans_trans',
    'IC_Asp_StrandParallel',
    'IC_Asp_StrandParallel_gminus_gminus',
    'IC_Asp_StrandParallel_gminus_gplus',
    'IC_Asp_StrandParallel_gminus_trans',
    'IC_Asp_StrandParallel_gplus_gminus',
    'IC_Asp_StrandParallel_gplus_gplus',
    'IC_Asp_StrandParallel_gplus_trans',
    'IC_Asp_StrandParallel_trans_gminus',
    'IC_Asp_StrandParallel_trans_gplus',
    'IC_Asp_StrandParallel_trans_trans',
    'IC_Asp_Strand_gminus_gminus',
    'IC_Asp_Strand_gminus_gplus',
    'IC_Asp_Strand_gminus_trans',
    'IC_Asp_Strand_gplus_gminus',
    'IC_Asp_Strand_gplus_gplus',
    'IC_Asp_Strand_gplus_trans',
    'IC_Asp_Strand_trans_gminus',
    'IC_Asp_Strand_trans_gplus',
    'IC_Asp_Strand_trans_trans',
    'IC_Asp_Turn',
    'IC_Asp_Turn_gminus_gminus',
    'IC_Asp_Turn_gminus_gplus',
    'IC_Asp_Turn_gminus_trans',
    'IC_Asp_Turn_gplus_gminus',
    'IC_Asp_Turn_gplus_gplus',
    'IC_Asp_Turn_gplus_trans',
    'IC_Asp_Turn_trans_gminus',
    'IC_Asp_Turn_trans_gplus',
    'IC_Asp_Turn_trans_trans',
]
