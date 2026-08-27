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

IC_Gln_HelixAlpha = build_template('GLN', 'alpha-helix', 'canonical')
IC_Gln_HelixAlpha_gminus_gminus = build_template('GLN', 'alpha-helix', 'g-/g-')
IC_Gln_HelixAlpha_gminus_trans = build_template('GLN', 'alpha-helix', 'g-/t')
IC_Gln_HelixAlpha_gminus_gplus = build_template('GLN', 'alpha-helix', 'g-/g+')
IC_Gln_HelixAlpha_trans_gminus = build_template('GLN', 'alpha-helix', 't/g-')
IC_Gln_HelixAlpha_trans_trans = build_template('GLN', 'alpha-helix', 't/t')
IC_Gln_HelixAlpha_trans_gplus = build_template('GLN', 'alpha-helix', 't/g+')
IC_Gln_HelixAlpha_gplus_gminus = build_template('GLN', 'alpha-helix', 'g+/g-')
IC_Gln_HelixAlpha_gplus_trans = build_template('GLN', 'alpha-helix', 'g+/t')
IC_Gln_HelixAlpha_gplus_gplus = build_template('GLN', 'alpha-helix', 'g+/g+')
IC_Gln_Helix310 = build_template('GLN', '3-10-helix', 'canonical')
IC_Gln_Helix310_gminus_gminus = build_template('GLN', '3-10-helix', 'g-/g-')
IC_Gln_Helix310_gminus_trans = build_template('GLN', '3-10-helix', 'g-/t')
IC_Gln_Helix310_gminus_gplus = build_template('GLN', '3-10-helix', 'g-/g+')
IC_Gln_Helix310_trans_gminus = build_template('GLN', '3-10-helix', 't/g-')
IC_Gln_Helix310_trans_trans = build_template('GLN', '3-10-helix', 't/t')
IC_Gln_Helix310_trans_gplus = build_template('GLN', '3-10-helix', 't/g+')
IC_Gln_Helix310_gplus_gminus = build_template('GLN', '3-10-helix', 'g+/g-')
IC_Gln_Helix310_gplus_trans = build_template('GLN', '3-10-helix', 'g+/t')
IC_Gln_Helix310_gplus_gplus = build_template('GLN', '3-10-helix', 'g+/g+')
IC_Gln_HelixPi = build_template('GLN', 'pi-helix', 'canonical')
IC_Gln_HelixPi_gminus_gminus = build_template('GLN', 'pi-helix', 'g-/g-')
IC_Gln_HelixPi_gminus_trans = build_template('GLN', 'pi-helix', 'g-/t')
IC_Gln_HelixPi_gminus_gplus = build_template('GLN', 'pi-helix', 'g-/g+')
IC_Gln_HelixPi_trans_gminus = build_template('GLN', 'pi-helix', 't/g-')
IC_Gln_HelixPi_trans_trans = build_template('GLN', 'pi-helix', 't/t')
IC_Gln_HelixPi_trans_gplus = build_template('GLN', 'pi-helix', 't/g+')
IC_Gln_HelixPi_gplus_gminus = build_template('GLN', 'pi-helix', 'g+/g-')
IC_Gln_HelixPi_gplus_trans = build_template('GLN', 'pi-helix', 'g+/t')
IC_Gln_HelixPi_gplus_gplus = build_template('GLN', 'pi-helix', 'g+/g+')
IC_Gln_HelixPPII = build_template('GLN', 'polyproline-II', 'canonical')
IC_Gln_HelixPPII_gminus_gminus = build_template('GLN', 'polyproline-II', 'g-/g-')
IC_Gln_HelixPPII_gminus_trans = build_template('GLN', 'polyproline-II', 'g-/t')
IC_Gln_HelixPPII_gminus_gplus = build_template('GLN', 'polyproline-II', 'g-/g+')
IC_Gln_HelixPPII_trans_gminus = build_template('GLN', 'polyproline-II', 't/g-')
IC_Gln_HelixPPII_trans_trans = build_template('GLN', 'polyproline-II', 't/t')
IC_Gln_HelixPPII_trans_gplus = build_template('GLN', 'polyproline-II', 't/g+')
IC_Gln_HelixPPII_gplus_gminus = build_template('GLN', 'polyproline-II', 'g+/g-')
IC_Gln_HelixPPII_gplus_trans = build_template('GLN', 'polyproline-II', 'g+/t')
IC_Gln_HelixPPII_gplus_gplus = build_template('GLN', 'polyproline-II', 'g+/g+')
IC_Gln_Strand = build_template('GLN', 'beta-strand', 'canonical')
IC_Gln_Strand_gminus_gminus = build_template('GLN', 'beta-strand', 'g-/g-')
IC_Gln_Strand_gminus_trans = build_template('GLN', 'beta-strand', 'g-/t')
IC_Gln_Strand_gminus_gplus = build_template('GLN', 'beta-strand', 'g-/g+')
IC_Gln_Strand_trans_gminus = build_template('GLN', 'beta-strand', 't/g-')
IC_Gln_Strand_trans_trans = build_template('GLN', 'beta-strand', 't/t')
IC_Gln_Strand_trans_gplus = build_template('GLN', 'beta-strand', 't/g+')
IC_Gln_Strand_gplus_gminus = build_template('GLN', 'beta-strand', 'g+/g-')
IC_Gln_Strand_gplus_trans = build_template('GLN', 'beta-strand', 'g+/t')
IC_Gln_Strand_gplus_gplus = build_template('GLN', 'beta-strand', 'g+/g+')
IC_Gln_StrandParallel = build_template('GLN', 'parallel-beta-strand', 'canonical')
IC_Gln_StrandParallel_gminus_gminus = build_template('GLN', 'parallel-beta-strand', 'g-/g-')
IC_Gln_StrandParallel_gminus_trans = build_template('GLN', 'parallel-beta-strand', 'g-/t')
IC_Gln_StrandParallel_gminus_gplus = build_template('GLN', 'parallel-beta-strand', 'g-/g+')
IC_Gln_StrandParallel_trans_gminus = build_template('GLN', 'parallel-beta-strand', 't/g-')
IC_Gln_StrandParallel_trans_trans = build_template('GLN', 'parallel-beta-strand', 't/t')
IC_Gln_StrandParallel_trans_gplus = build_template('GLN', 'parallel-beta-strand', 't/g+')
IC_Gln_StrandParallel_gplus_gminus = build_template('GLN', 'parallel-beta-strand', 'g+/g-')
IC_Gln_StrandParallel_gplus_trans = build_template('GLN', 'parallel-beta-strand', 'g+/t')
IC_Gln_StrandParallel_gplus_gplus = build_template('GLN', 'parallel-beta-strand', 'g+/g+')
IC_Gln_StrandAntiParallel = build_template('GLN', 'antiparallel-beta-strand', 'canonical')
IC_Gln_StrandAntiParallel_gminus_gminus = build_template('GLN', 'antiparallel-beta-strand', 'g-/g-')
IC_Gln_StrandAntiParallel_gminus_trans = build_template('GLN', 'antiparallel-beta-strand', 'g-/t')
IC_Gln_StrandAntiParallel_gminus_gplus = build_template('GLN', 'antiparallel-beta-strand', 'g-/g+')
IC_Gln_StrandAntiParallel_trans_gminus = build_template('GLN', 'antiparallel-beta-strand', 't/g-')
IC_Gln_StrandAntiParallel_trans_trans = build_template('GLN', 'antiparallel-beta-strand', 't/t')
IC_Gln_StrandAntiParallel_trans_gplus = build_template('GLN', 'antiparallel-beta-strand', 't/g+')
IC_Gln_StrandAntiParallel_gplus_gminus = build_template('GLN', 'antiparallel-beta-strand', 'g+/g-')
IC_Gln_StrandAntiParallel_gplus_trans = build_template('GLN', 'antiparallel-beta-strand', 'g+/t')
IC_Gln_StrandAntiParallel_gplus_gplus = build_template('GLN', 'antiparallel-beta-strand', 'g+/g+')
IC_Gln_Bridge = build_template('GLN', 'beta-bridge', 'canonical')
IC_Gln_Bridge_gminus_gminus = build_template('GLN', 'beta-bridge', 'g-/g-')
IC_Gln_Bridge_gminus_trans = build_template('GLN', 'beta-bridge', 'g-/t')
IC_Gln_Bridge_gminus_gplus = build_template('GLN', 'beta-bridge', 'g-/g+')
IC_Gln_Bridge_trans_gminus = build_template('GLN', 'beta-bridge', 't/g-')
IC_Gln_Bridge_trans_trans = build_template('GLN', 'beta-bridge', 't/t')
IC_Gln_Bridge_trans_gplus = build_template('GLN', 'beta-bridge', 't/g+')
IC_Gln_Bridge_gplus_gminus = build_template('GLN', 'beta-bridge', 'g+/g-')
IC_Gln_Bridge_gplus_trans = build_template('GLN', 'beta-bridge', 'g+/t')
IC_Gln_Bridge_gplus_gplus = build_template('GLN', 'beta-bridge', 'g+/g+')
IC_Gln_Turn = build_template('GLN', 'turn', 'canonical')
IC_Gln_Turn_gminus_gminus = build_template('GLN', 'turn', 'g-/g-')
IC_Gln_Turn_gminus_trans = build_template('GLN', 'turn', 'g-/t')
IC_Gln_Turn_gminus_gplus = build_template('GLN', 'turn', 'g-/g+')
IC_Gln_Turn_trans_gminus = build_template('GLN', 'turn', 't/g-')
IC_Gln_Turn_trans_trans = build_template('GLN', 'turn', 't/t')
IC_Gln_Turn_trans_gplus = build_template('GLN', 'turn', 't/g+')
IC_Gln_Turn_gplus_gminus = build_template('GLN', 'turn', 'g+/g-')
IC_Gln_Turn_gplus_trans = build_template('GLN', 'turn', 'g+/t')
IC_Gln_Turn_gplus_gplus = build_template('GLN', 'turn', 'g+/g+')
IC_Gln_Bend = build_template('GLN', 'bend', 'canonical')
IC_Gln_Bend_gminus_gminus = build_template('GLN', 'bend', 'g-/g-')
IC_Gln_Bend_gminus_trans = build_template('GLN', 'bend', 'g-/t')
IC_Gln_Bend_gminus_gplus = build_template('GLN', 'bend', 'g-/g+')
IC_Gln_Bend_trans_gminus = build_template('GLN', 'bend', 't/g-')
IC_Gln_Bend_trans_trans = build_template('GLN', 'bend', 't/t')
IC_Gln_Bend_trans_gplus = build_template('GLN', 'bend', 't/g+')
IC_Gln_Bend_gplus_gminus = build_template('GLN', 'bend', 'g+/g-')
IC_Gln_Bend_gplus_trans = build_template('GLN', 'bend', 'g+/t')
IC_Gln_Bend_gplus_gplus = build_template('GLN', 'bend', 'g+/g+')
IC_Gln_Coil = build_template('GLN', 'coil', 'canonical')
IC_Gln_Coil_gminus_gminus = build_template('GLN', 'coil', 'g-/g-')
IC_Gln_Coil_gminus_trans = build_template('GLN', 'coil', 'g-/t')
IC_Gln_Coil_gminus_gplus = build_template('GLN', 'coil', 'g-/g+')
IC_Gln_Coil_trans_gminus = build_template('GLN', 'coil', 't/g-')
IC_Gln_Coil_trans_trans = build_template('GLN', 'coil', 't/t')
IC_Gln_Coil_trans_gplus = build_template('GLN', 'coil', 't/g+')
IC_Gln_Coil_gplus_gminus = build_template('GLN', 'coil', 'g+/g-')
IC_Gln_Coil_gplus_trans = build_template('GLN', 'coil', 'g+/t')
IC_Gln_Coil_gplus_gplus = build_template('GLN', 'coil', 'g+/g+')
IC_Gln_CisPeptide = build_template('GLN', 'cis-peptide-bond', 'canonical')
IC_Gln_CisPeptide_gminus_gminus = build_template('GLN', 'cis-peptide-bond', 'g-/g-')
IC_Gln_CisPeptide_gminus_trans = build_template('GLN', 'cis-peptide-bond', 'g-/t')
IC_Gln_CisPeptide_gminus_gplus = build_template('GLN', 'cis-peptide-bond', 'g-/g+')
IC_Gln_CisPeptide_trans_gminus = build_template('GLN', 'cis-peptide-bond', 't/g-')
IC_Gln_CisPeptide_trans_trans = build_template('GLN', 'cis-peptide-bond', 't/t')
IC_Gln_CisPeptide_trans_gplus = build_template('GLN', 'cis-peptide-bond', 't/g+')
IC_Gln_CisPeptide_gplus_gminus = build_template('GLN', 'cis-peptide-bond', 'g+/g-')
IC_Gln_CisPeptide_gplus_trans = build_template('GLN', 'cis-peptide-bond', 'g+/t')
IC_Gln_CisPeptide_gplus_gplus = build_template('GLN', 'cis-peptide-bond', 'g+/g+')

__all__ = [
    'IC_Gln_Bend',
    'IC_Gln_Bend_gminus_gminus',
    'IC_Gln_Bend_gminus_gplus',
    'IC_Gln_Bend_gminus_trans',
    'IC_Gln_Bend_gplus_gminus',
    'IC_Gln_Bend_gplus_gplus',
    'IC_Gln_Bend_gplus_trans',
    'IC_Gln_Bend_trans_gminus',
    'IC_Gln_Bend_trans_gplus',
    'IC_Gln_Bend_trans_trans',
    'IC_Gln_Bridge',
    'IC_Gln_Bridge_gminus_gminus',
    'IC_Gln_Bridge_gminus_gplus',
    'IC_Gln_Bridge_gminus_trans',
    'IC_Gln_Bridge_gplus_gminus',
    'IC_Gln_Bridge_gplus_gplus',
    'IC_Gln_Bridge_gplus_trans',
    'IC_Gln_Bridge_trans_gminus',
    'IC_Gln_Bridge_trans_gplus',
    'IC_Gln_Bridge_trans_trans',
    'IC_Gln_CisPeptide',
    'IC_Gln_CisPeptide_gminus_gminus',
    'IC_Gln_CisPeptide_gminus_gplus',
    'IC_Gln_CisPeptide_gminus_trans',
    'IC_Gln_CisPeptide_gplus_gminus',
    'IC_Gln_CisPeptide_gplus_gplus',
    'IC_Gln_CisPeptide_gplus_trans',
    'IC_Gln_CisPeptide_trans_gminus',
    'IC_Gln_CisPeptide_trans_gplus',
    'IC_Gln_CisPeptide_trans_trans',
    'IC_Gln_Coil',
    'IC_Gln_Coil_gminus_gminus',
    'IC_Gln_Coil_gminus_gplus',
    'IC_Gln_Coil_gminus_trans',
    'IC_Gln_Coil_gplus_gminus',
    'IC_Gln_Coil_gplus_gplus',
    'IC_Gln_Coil_gplus_trans',
    'IC_Gln_Coil_trans_gminus',
    'IC_Gln_Coil_trans_gplus',
    'IC_Gln_Coil_trans_trans',
    'IC_Gln_Helix310',
    'IC_Gln_Helix310_gminus_gminus',
    'IC_Gln_Helix310_gminus_gplus',
    'IC_Gln_Helix310_gminus_trans',
    'IC_Gln_Helix310_gplus_gminus',
    'IC_Gln_Helix310_gplus_gplus',
    'IC_Gln_Helix310_gplus_trans',
    'IC_Gln_Helix310_trans_gminus',
    'IC_Gln_Helix310_trans_gplus',
    'IC_Gln_Helix310_trans_trans',
    'IC_Gln_HelixAlpha',
    'IC_Gln_HelixAlpha_gminus_gminus',
    'IC_Gln_HelixAlpha_gminus_gplus',
    'IC_Gln_HelixAlpha_gminus_trans',
    'IC_Gln_HelixAlpha_gplus_gminus',
    'IC_Gln_HelixAlpha_gplus_gplus',
    'IC_Gln_HelixAlpha_gplus_trans',
    'IC_Gln_HelixAlpha_trans_gminus',
    'IC_Gln_HelixAlpha_trans_gplus',
    'IC_Gln_HelixAlpha_trans_trans',
    'IC_Gln_HelixPPII',
    'IC_Gln_HelixPPII_gminus_gminus',
    'IC_Gln_HelixPPII_gminus_gplus',
    'IC_Gln_HelixPPII_gminus_trans',
    'IC_Gln_HelixPPII_gplus_gminus',
    'IC_Gln_HelixPPII_gplus_gplus',
    'IC_Gln_HelixPPII_gplus_trans',
    'IC_Gln_HelixPPII_trans_gminus',
    'IC_Gln_HelixPPII_trans_gplus',
    'IC_Gln_HelixPPII_trans_trans',
    'IC_Gln_HelixPi',
    'IC_Gln_HelixPi_gminus_gminus',
    'IC_Gln_HelixPi_gminus_gplus',
    'IC_Gln_HelixPi_gminus_trans',
    'IC_Gln_HelixPi_gplus_gminus',
    'IC_Gln_HelixPi_gplus_gplus',
    'IC_Gln_HelixPi_gplus_trans',
    'IC_Gln_HelixPi_trans_gminus',
    'IC_Gln_HelixPi_trans_gplus',
    'IC_Gln_HelixPi_trans_trans',
    'IC_Gln_Strand',
    'IC_Gln_StrandAntiParallel',
    'IC_Gln_StrandAntiParallel_gminus_gminus',
    'IC_Gln_StrandAntiParallel_gminus_gplus',
    'IC_Gln_StrandAntiParallel_gminus_trans',
    'IC_Gln_StrandAntiParallel_gplus_gminus',
    'IC_Gln_StrandAntiParallel_gplus_gplus',
    'IC_Gln_StrandAntiParallel_gplus_trans',
    'IC_Gln_StrandAntiParallel_trans_gminus',
    'IC_Gln_StrandAntiParallel_trans_gplus',
    'IC_Gln_StrandAntiParallel_trans_trans',
    'IC_Gln_StrandParallel',
    'IC_Gln_StrandParallel_gminus_gminus',
    'IC_Gln_StrandParallel_gminus_gplus',
    'IC_Gln_StrandParallel_gminus_trans',
    'IC_Gln_StrandParallel_gplus_gminus',
    'IC_Gln_StrandParallel_gplus_gplus',
    'IC_Gln_StrandParallel_gplus_trans',
    'IC_Gln_StrandParallel_trans_gminus',
    'IC_Gln_StrandParallel_trans_gplus',
    'IC_Gln_StrandParallel_trans_trans',
    'IC_Gln_Strand_gminus_gminus',
    'IC_Gln_Strand_gminus_gplus',
    'IC_Gln_Strand_gminus_trans',
    'IC_Gln_Strand_gplus_gminus',
    'IC_Gln_Strand_gplus_gplus',
    'IC_Gln_Strand_gplus_trans',
    'IC_Gln_Strand_trans_gminus',
    'IC_Gln_Strand_trans_gplus',
    'IC_Gln_Strand_trans_trans',
    'IC_Gln_Turn',
    'IC_Gln_Turn_gminus_gminus',
    'IC_Gln_Turn_gminus_gplus',
    'IC_Gln_Turn_gminus_trans',
    'IC_Gln_Turn_gplus_gminus',
    'IC_Gln_Turn_gplus_gplus',
    'IC_Gln_Turn_gplus_trans',
    'IC_Gln_Turn_trans_gminus',
    'IC_Gln_Turn_trans_gplus',
    'IC_Gln_Turn_trans_trans',
]
