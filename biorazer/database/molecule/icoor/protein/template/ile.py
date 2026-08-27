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

IC_Ile_HelixAlpha = build_template('ILE', 'alpha-helix', 'canonical')
IC_Ile_HelixAlpha_gminus_gminus = build_template('ILE', 'alpha-helix', 'g-/g-')
IC_Ile_HelixAlpha_gminus_trans = build_template('ILE', 'alpha-helix', 'g-/t')
IC_Ile_HelixAlpha_gminus_gplus = build_template('ILE', 'alpha-helix', 'g-/g+')
IC_Ile_HelixAlpha_trans_gminus = build_template('ILE', 'alpha-helix', 't/g-')
IC_Ile_HelixAlpha_trans_trans = build_template('ILE', 'alpha-helix', 't/t')
IC_Ile_HelixAlpha_trans_gplus = build_template('ILE', 'alpha-helix', 't/g+')
IC_Ile_HelixAlpha_gplus_gminus = build_template('ILE', 'alpha-helix', 'g+/g-')
IC_Ile_HelixAlpha_gplus_trans = build_template('ILE', 'alpha-helix', 'g+/t')
IC_Ile_HelixAlpha_gplus_gplus = build_template('ILE', 'alpha-helix', 'g+/g+')
IC_Ile_Helix310 = build_template('ILE', '3-10-helix', 'canonical')
IC_Ile_Helix310_gminus_gminus = build_template('ILE', '3-10-helix', 'g-/g-')
IC_Ile_Helix310_gminus_trans = build_template('ILE', '3-10-helix', 'g-/t')
IC_Ile_Helix310_gminus_gplus = build_template('ILE', '3-10-helix', 'g-/g+')
IC_Ile_Helix310_trans_gminus = build_template('ILE', '3-10-helix', 't/g-')
IC_Ile_Helix310_trans_trans = build_template('ILE', '3-10-helix', 't/t')
IC_Ile_Helix310_trans_gplus = build_template('ILE', '3-10-helix', 't/g+')
IC_Ile_Helix310_gplus_gminus = build_template('ILE', '3-10-helix', 'g+/g-')
IC_Ile_Helix310_gplus_trans = build_template('ILE', '3-10-helix', 'g+/t')
IC_Ile_Helix310_gplus_gplus = build_template('ILE', '3-10-helix', 'g+/g+')
IC_Ile_HelixPi = build_template('ILE', 'pi-helix', 'canonical')
IC_Ile_HelixPi_gminus_gminus = build_template('ILE', 'pi-helix', 'g-/g-')
IC_Ile_HelixPi_gminus_trans = build_template('ILE', 'pi-helix', 'g-/t')
IC_Ile_HelixPi_gminus_gplus = build_template('ILE', 'pi-helix', 'g-/g+')
IC_Ile_HelixPi_trans_gminus = build_template('ILE', 'pi-helix', 't/g-')
IC_Ile_HelixPi_trans_trans = build_template('ILE', 'pi-helix', 't/t')
IC_Ile_HelixPi_trans_gplus = build_template('ILE', 'pi-helix', 't/g+')
IC_Ile_HelixPi_gplus_gminus = build_template('ILE', 'pi-helix', 'g+/g-')
IC_Ile_HelixPi_gplus_trans = build_template('ILE', 'pi-helix', 'g+/t')
IC_Ile_HelixPi_gplus_gplus = build_template('ILE', 'pi-helix', 'g+/g+')
IC_Ile_HelixPPII = build_template('ILE', 'polyproline-II', 'canonical')
IC_Ile_HelixPPII_gminus_gminus = build_template('ILE', 'polyproline-II', 'g-/g-')
IC_Ile_HelixPPII_gminus_trans = build_template('ILE', 'polyproline-II', 'g-/t')
IC_Ile_HelixPPII_gminus_gplus = build_template('ILE', 'polyproline-II', 'g-/g+')
IC_Ile_HelixPPII_trans_gminus = build_template('ILE', 'polyproline-II', 't/g-')
IC_Ile_HelixPPII_trans_trans = build_template('ILE', 'polyproline-II', 't/t')
IC_Ile_HelixPPII_trans_gplus = build_template('ILE', 'polyproline-II', 't/g+')
IC_Ile_HelixPPII_gplus_gminus = build_template('ILE', 'polyproline-II', 'g+/g-')
IC_Ile_HelixPPII_gplus_trans = build_template('ILE', 'polyproline-II', 'g+/t')
IC_Ile_HelixPPII_gplus_gplus = build_template('ILE', 'polyproline-II', 'g+/g+')
IC_Ile_Strand = build_template('ILE', 'beta-strand', 'canonical')
IC_Ile_Strand_gminus_gminus = build_template('ILE', 'beta-strand', 'g-/g-')
IC_Ile_Strand_gminus_trans = build_template('ILE', 'beta-strand', 'g-/t')
IC_Ile_Strand_gminus_gplus = build_template('ILE', 'beta-strand', 'g-/g+')
IC_Ile_Strand_trans_gminus = build_template('ILE', 'beta-strand', 't/g-')
IC_Ile_Strand_trans_trans = build_template('ILE', 'beta-strand', 't/t')
IC_Ile_Strand_trans_gplus = build_template('ILE', 'beta-strand', 't/g+')
IC_Ile_Strand_gplus_gminus = build_template('ILE', 'beta-strand', 'g+/g-')
IC_Ile_Strand_gplus_trans = build_template('ILE', 'beta-strand', 'g+/t')
IC_Ile_Strand_gplus_gplus = build_template('ILE', 'beta-strand', 'g+/g+')
IC_Ile_StrandParallel = build_template('ILE', 'parallel-beta-strand', 'canonical')
IC_Ile_StrandParallel_gminus_gminus = build_template('ILE', 'parallel-beta-strand', 'g-/g-')
IC_Ile_StrandParallel_gminus_trans = build_template('ILE', 'parallel-beta-strand', 'g-/t')
IC_Ile_StrandParallel_gminus_gplus = build_template('ILE', 'parallel-beta-strand', 'g-/g+')
IC_Ile_StrandParallel_trans_gminus = build_template('ILE', 'parallel-beta-strand', 't/g-')
IC_Ile_StrandParallel_trans_trans = build_template('ILE', 'parallel-beta-strand', 't/t')
IC_Ile_StrandParallel_trans_gplus = build_template('ILE', 'parallel-beta-strand', 't/g+')
IC_Ile_StrandParallel_gplus_gminus = build_template('ILE', 'parallel-beta-strand', 'g+/g-')
IC_Ile_StrandParallel_gplus_trans = build_template('ILE', 'parallel-beta-strand', 'g+/t')
IC_Ile_StrandParallel_gplus_gplus = build_template('ILE', 'parallel-beta-strand', 'g+/g+')
IC_Ile_StrandAntiParallel = build_template('ILE', 'antiparallel-beta-strand', 'canonical')
IC_Ile_StrandAntiParallel_gminus_gminus = build_template('ILE', 'antiparallel-beta-strand', 'g-/g-')
IC_Ile_StrandAntiParallel_gminus_trans = build_template('ILE', 'antiparallel-beta-strand', 'g-/t')
IC_Ile_StrandAntiParallel_gminus_gplus = build_template('ILE', 'antiparallel-beta-strand', 'g-/g+')
IC_Ile_StrandAntiParallel_trans_gminus = build_template('ILE', 'antiparallel-beta-strand', 't/g-')
IC_Ile_StrandAntiParallel_trans_trans = build_template('ILE', 'antiparallel-beta-strand', 't/t')
IC_Ile_StrandAntiParallel_trans_gplus = build_template('ILE', 'antiparallel-beta-strand', 't/g+')
IC_Ile_StrandAntiParallel_gplus_gminus = build_template('ILE', 'antiparallel-beta-strand', 'g+/g-')
IC_Ile_StrandAntiParallel_gplus_trans = build_template('ILE', 'antiparallel-beta-strand', 'g+/t')
IC_Ile_StrandAntiParallel_gplus_gplus = build_template('ILE', 'antiparallel-beta-strand', 'g+/g+')
IC_Ile_Bridge = build_template('ILE', 'beta-bridge', 'canonical')
IC_Ile_Bridge_gminus_gminus = build_template('ILE', 'beta-bridge', 'g-/g-')
IC_Ile_Bridge_gminus_trans = build_template('ILE', 'beta-bridge', 'g-/t')
IC_Ile_Bridge_gminus_gplus = build_template('ILE', 'beta-bridge', 'g-/g+')
IC_Ile_Bridge_trans_gminus = build_template('ILE', 'beta-bridge', 't/g-')
IC_Ile_Bridge_trans_trans = build_template('ILE', 'beta-bridge', 't/t')
IC_Ile_Bridge_trans_gplus = build_template('ILE', 'beta-bridge', 't/g+')
IC_Ile_Bridge_gplus_gminus = build_template('ILE', 'beta-bridge', 'g+/g-')
IC_Ile_Bridge_gplus_trans = build_template('ILE', 'beta-bridge', 'g+/t')
IC_Ile_Bridge_gplus_gplus = build_template('ILE', 'beta-bridge', 'g+/g+')
IC_Ile_Turn = build_template('ILE', 'turn', 'canonical')
IC_Ile_Turn_gminus_gminus = build_template('ILE', 'turn', 'g-/g-')
IC_Ile_Turn_gminus_trans = build_template('ILE', 'turn', 'g-/t')
IC_Ile_Turn_gminus_gplus = build_template('ILE', 'turn', 'g-/g+')
IC_Ile_Turn_trans_gminus = build_template('ILE', 'turn', 't/g-')
IC_Ile_Turn_trans_trans = build_template('ILE', 'turn', 't/t')
IC_Ile_Turn_trans_gplus = build_template('ILE', 'turn', 't/g+')
IC_Ile_Turn_gplus_gminus = build_template('ILE', 'turn', 'g+/g-')
IC_Ile_Turn_gplus_trans = build_template('ILE', 'turn', 'g+/t')
IC_Ile_Turn_gplus_gplus = build_template('ILE', 'turn', 'g+/g+')
IC_Ile_Bend = build_template('ILE', 'bend', 'canonical')
IC_Ile_Bend_gminus_gminus = build_template('ILE', 'bend', 'g-/g-')
IC_Ile_Bend_gminus_trans = build_template('ILE', 'bend', 'g-/t')
IC_Ile_Bend_gminus_gplus = build_template('ILE', 'bend', 'g-/g+')
IC_Ile_Bend_trans_gminus = build_template('ILE', 'bend', 't/g-')
IC_Ile_Bend_trans_trans = build_template('ILE', 'bend', 't/t')
IC_Ile_Bend_trans_gplus = build_template('ILE', 'bend', 't/g+')
IC_Ile_Bend_gplus_gminus = build_template('ILE', 'bend', 'g+/g-')
IC_Ile_Bend_gplus_trans = build_template('ILE', 'bend', 'g+/t')
IC_Ile_Bend_gplus_gplus = build_template('ILE', 'bend', 'g+/g+')
IC_Ile_Coil = build_template('ILE', 'coil', 'canonical')
IC_Ile_Coil_gminus_gminus = build_template('ILE', 'coil', 'g-/g-')
IC_Ile_Coil_gminus_trans = build_template('ILE', 'coil', 'g-/t')
IC_Ile_Coil_gminus_gplus = build_template('ILE', 'coil', 'g-/g+')
IC_Ile_Coil_trans_gminus = build_template('ILE', 'coil', 't/g-')
IC_Ile_Coil_trans_trans = build_template('ILE', 'coil', 't/t')
IC_Ile_Coil_trans_gplus = build_template('ILE', 'coil', 't/g+')
IC_Ile_Coil_gplus_gminus = build_template('ILE', 'coil', 'g+/g-')
IC_Ile_Coil_gplus_trans = build_template('ILE', 'coil', 'g+/t')
IC_Ile_Coil_gplus_gplus = build_template('ILE', 'coil', 'g+/g+')
IC_Ile_CisPeptide = build_template('ILE', 'cis-peptide-bond', 'canonical')
IC_Ile_CisPeptide_gminus_gminus = build_template('ILE', 'cis-peptide-bond', 'g-/g-')
IC_Ile_CisPeptide_gminus_trans = build_template('ILE', 'cis-peptide-bond', 'g-/t')
IC_Ile_CisPeptide_gminus_gplus = build_template('ILE', 'cis-peptide-bond', 'g-/g+')
IC_Ile_CisPeptide_trans_gminus = build_template('ILE', 'cis-peptide-bond', 't/g-')
IC_Ile_CisPeptide_trans_trans = build_template('ILE', 'cis-peptide-bond', 't/t')
IC_Ile_CisPeptide_trans_gplus = build_template('ILE', 'cis-peptide-bond', 't/g+')
IC_Ile_CisPeptide_gplus_gminus = build_template('ILE', 'cis-peptide-bond', 'g+/g-')
IC_Ile_CisPeptide_gplus_trans = build_template('ILE', 'cis-peptide-bond', 'g+/t')
IC_Ile_CisPeptide_gplus_gplus = build_template('ILE', 'cis-peptide-bond', 'g+/g+')

__all__ = [
    'IC_Ile_Bend',
    'IC_Ile_Bend_gminus_gminus',
    'IC_Ile_Bend_gminus_gplus',
    'IC_Ile_Bend_gminus_trans',
    'IC_Ile_Bend_gplus_gminus',
    'IC_Ile_Bend_gplus_gplus',
    'IC_Ile_Bend_gplus_trans',
    'IC_Ile_Bend_trans_gminus',
    'IC_Ile_Bend_trans_gplus',
    'IC_Ile_Bend_trans_trans',
    'IC_Ile_Bridge',
    'IC_Ile_Bridge_gminus_gminus',
    'IC_Ile_Bridge_gminus_gplus',
    'IC_Ile_Bridge_gminus_trans',
    'IC_Ile_Bridge_gplus_gminus',
    'IC_Ile_Bridge_gplus_gplus',
    'IC_Ile_Bridge_gplus_trans',
    'IC_Ile_Bridge_trans_gminus',
    'IC_Ile_Bridge_trans_gplus',
    'IC_Ile_Bridge_trans_trans',
    'IC_Ile_CisPeptide',
    'IC_Ile_CisPeptide_gminus_gminus',
    'IC_Ile_CisPeptide_gminus_gplus',
    'IC_Ile_CisPeptide_gminus_trans',
    'IC_Ile_CisPeptide_gplus_gminus',
    'IC_Ile_CisPeptide_gplus_gplus',
    'IC_Ile_CisPeptide_gplus_trans',
    'IC_Ile_CisPeptide_trans_gminus',
    'IC_Ile_CisPeptide_trans_gplus',
    'IC_Ile_CisPeptide_trans_trans',
    'IC_Ile_Coil',
    'IC_Ile_Coil_gminus_gminus',
    'IC_Ile_Coil_gminus_gplus',
    'IC_Ile_Coil_gminus_trans',
    'IC_Ile_Coil_gplus_gminus',
    'IC_Ile_Coil_gplus_gplus',
    'IC_Ile_Coil_gplus_trans',
    'IC_Ile_Coil_trans_gminus',
    'IC_Ile_Coil_trans_gplus',
    'IC_Ile_Coil_trans_trans',
    'IC_Ile_Helix310',
    'IC_Ile_Helix310_gminus_gminus',
    'IC_Ile_Helix310_gminus_gplus',
    'IC_Ile_Helix310_gminus_trans',
    'IC_Ile_Helix310_gplus_gminus',
    'IC_Ile_Helix310_gplus_gplus',
    'IC_Ile_Helix310_gplus_trans',
    'IC_Ile_Helix310_trans_gminus',
    'IC_Ile_Helix310_trans_gplus',
    'IC_Ile_Helix310_trans_trans',
    'IC_Ile_HelixAlpha',
    'IC_Ile_HelixAlpha_gminus_gminus',
    'IC_Ile_HelixAlpha_gminus_gplus',
    'IC_Ile_HelixAlpha_gminus_trans',
    'IC_Ile_HelixAlpha_gplus_gminus',
    'IC_Ile_HelixAlpha_gplus_gplus',
    'IC_Ile_HelixAlpha_gplus_trans',
    'IC_Ile_HelixAlpha_trans_gminus',
    'IC_Ile_HelixAlpha_trans_gplus',
    'IC_Ile_HelixAlpha_trans_trans',
    'IC_Ile_HelixPPII',
    'IC_Ile_HelixPPII_gminus_gminus',
    'IC_Ile_HelixPPII_gminus_gplus',
    'IC_Ile_HelixPPII_gminus_trans',
    'IC_Ile_HelixPPII_gplus_gminus',
    'IC_Ile_HelixPPII_gplus_gplus',
    'IC_Ile_HelixPPII_gplus_trans',
    'IC_Ile_HelixPPII_trans_gminus',
    'IC_Ile_HelixPPII_trans_gplus',
    'IC_Ile_HelixPPII_trans_trans',
    'IC_Ile_HelixPi',
    'IC_Ile_HelixPi_gminus_gminus',
    'IC_Ile_HelixPi_gminus_gplus',
    'IC_Ile_HelixPi_gminus_trans',
    'IC_Ile_HelixPi_gplus_gminus',
    'IC_Ile_HelixPi_gplus_gplus',
    'IC_Ile_HelixPi_gplus_trans',
    'IC_Ile_HelixPi_trans_gminus',
    'IC_Ile_HelixPi_trans_gplus',
    'IC_Ile_HelixPi_trans_trans',
    'IC_Ile_Strand',
    'IC_Ile_StrandAntiParallel',
    'IC_Ile_StrandAntiParallel_gminus_gminus',
    'IC_Ile_StrandAntiParallel_gminus_gplus',
    'IC_Ile_StrandAntiParallel_gminus_trans',
    'IC_Ile_StrandAntiParallel_gplus_gminus',
    'IC_Ile_StrandAntiParallel_gplus_gplus',
    'IC_Ile_StrandAntiParallel_gplus_trans',
    'IC_Ile_StrandAntiParallel_trans_gminus',
    'IC_Ile_StrandAntiParallel_trans_gplus',
    'IC_Ile_StrandAntiParallel_trans_trans',
    'IC_Ile_StrandParallel',
    'IC_Ile_StrandParallel_gminus_gminus',
    'IC_Ile_StrandParallel_gminus_gplus',
    'IC_Ile_StrandParallel_gminus_trans',
    'IC_Ile_StrandParallel_gplus_gminus',
    'IC_Ile_StrandParallel_gplus_gplus',
    'IC_Ile_StrandParallel_gplus_trans',
    'IC_Ile_StrandParallel_trans_gminus',
    'IC_Ile_StrandParallel_trans_gplus',
    'IC_Ile_StrandParallel_trans_trans',
    'IC_Ile_Strand_gminus_gminus',
    'IC_Ile_Strand_gminus_gplus',
    'IC_Ile_Strand_gminus_trans',
    'IC_Ile_Strand_gplus_gminus',
    'IC_Ile_Strand_gplus_gplus',
    'IC_Ile_Strand_gplus_trans',
    'IC_Ile_Strand_trans_gminus',
    'IC_Ile_Strand_trans_gplus',
    'IC_Ile_Strand_trans_trans',
    'IC_Ile_Turn',
    'IC_Ile_Turn_gminus_gminus',
    'IC_Ile_Turn_gminus_gplus',
    'IC_Ile_Turn_gminus_trans',
    'IC_Ile_Turn_gplus_gminus',
    'IC_Ile_Turn_gplus_gplus',
    'IC_Ile_Turn_gplus_trans',
    'IC_Ile_Turn_trans_gminus',
    'IC_Ile_Turn_trans_gplus',
    'IC_Ile_Turn_trans_trans',
]
