# -*- coding: utf-8 -*-
"""Ideal {secondary-structure x rotamer} {names[aa]} {var} ``InternalCoord`` templates.

Each template is the ideal heavy-atom topology for {aa} anchored at
{{N, CA, C}}, with the side chain grown to ``rotamer`` and the
``phi``/``psi``/``omega`` of ``ss`` carried on the instance.  See
:mod:`biorazer.database.molecule.icoor.protein.template` and
:mod:`._builder` for the build rule and the per-conformer caveat.
"""

from biorazer.database.molecule.icoor.protein.template._builder import build_template

RESN = "MET"

# === GENERATED IC_* named templates - do not edit below ===
# Regenerate with: /opt/envs/BioRazer/bin/python scripts/generate_internal_coord_template_named.py

IC_Met_HelixAlpha = build_template('MET', 'alpha-helix', 'canonical')
IC_Met_HelixAlpha_gminus_gminus = build_template('MET', 'alpha-helix', 'g-/g-')
IC_Met_HelixAlpha_gminus_trans = build_template('MET', 'alpha-helix', 'g-/t')
IC_Met_HelixAlpha_gminus_gplus = build_template('MET', 'alpha-helix', 'g-/g+')
IC_Met_HelixAlpha_trans_gminus = build_template('MET', 'alpha-helix', 't/g-')
IC_Met_HelixAlpha_trans_trans = build_template('MET', 'alpha-helix', 't/t')
IC_Met_HelixAlpha_trans_gplus = build_template('MET', 'alpha-helix', 't/g+')
IC_Met_HelixAlpha_gplus_gminus = build_template('MET', 'alpha-helix', 'g+/g-')
IC_Met_HelixAlpha_gplus_trans = build_template('MET', 'alpha-helix', 'g+/t')
IC_Met_HelixAlpha_gplus_gplus = build_template('MET', 'alpha-helix', 'g+/g+')
IC_Met_Helix310 = build_template('MET', '3-10-helix', 'canonical')
IC_Met_Helix310_gminus_gminus = build_template('MET', '3-10-helix', 'g-/g-')
IC_Met_Helix310_gminus_trans = build_template('MET', '3-10-helix', 'g-/t')
IC_Met_Helix310_gminus_gplus = build_template('MET', '3-10-helix', 'g-/g+')
IC_Met_Helix310_trans_gminus = build_template('MET', '3-10-helix', 't/g-')
IC_Met_Helix310_trans_trans = build_template('MET', '3-10-helix', 't/t')
IC_Met_Helix310_trans_gplus = build_template('MET', '3-10-helix', 't/g+')
IC_Met_Helix310_gplus_gminus = build_template('MET', '3-10-helix', 'g+/g-')
IC_Met_Helix310_gplus_trans = build_template('MET', '3-10-helix', 'g+/t')
IC_Met_Helix310_gplus_gplus = build_template('MET', '3-10-helix', 'g+/g+')
IC_Met_HelixPi = build_template('MET', 'pi-helix', 'canonical')
IC_Met_HelixPi_gminus_gminus = build_template('MET', 'pi-helix', 'g-/g-')
IC_Met_HelixPi_gminus_trans = build_template('MET', 'pi-helix', 'g-/t')
IC_Met_HelixPi_gminus_gplus = build_template('MET', 'pi-helix', 'g-/g+')
IC_Met_HelixPi_trans_gminus = build_template('MET', 'pi-helix', 't/g-')
IC_Met_HelixPi_trans_trans = build_template('MET', 'pi-helix', 't/t')
IC_Met_HelixPi_trans_gplus = build_template('MET', 'pi-helix', 't/g+')
IC_Met_HelixPi_gplus_gminus = build_template('MET', 'pi-helix', 'g+/g-')
IC_Met_HelixPi_gplus_trans = build_template('MET', 'pi-helix', 'g+/t')
IC_Met_HelixPi_gplus_gplus = build_template('MET', 'pi-helix', 'g+/g+')
IC_Met_HelixPPII = build_template('MET', 'polyproline-II', 'canonical')
IC_Met_HelixPPII_gminus_gminus = build_template('MET', 'polyproline-II', 'g-/g-')
IC_Met_HelixPPII_gminus_trans = build_template('MET', 'polyproline-II', 'g-/t')
IC_Met_HelixPPII_gminus_gplus = build_template('MET', 'polyproline-II', 'g-/g+')
IC_Met_HelixPPII_trans_gminus = build_template('MET', 'polyproline-II', 't/g-')
IC_Met_HelixPPII_trans_trans = build_template('MET', 'polyproline-II', 't/t')
IC_Met_HelixPPII_trans_gplus = build_template('MET', 'polyproline-II', 't/g+')
IC_Met_HelixPPII_gplus_gminus = build_template('MET', 'polyproline-II', 'g+/g-')
IC_Met_HelixPPII_gplus_trans = build_template('MET', 'polyproline-II', 'g+/t')
IC_Met_HelixPPII_gplus_gplus = build_template('MET', 'polyproline-II', 'g+/g+')
IC_Met_Strand = build_template('MET', 'beta-strand', 'canonical')
IC_Met_Strand_gminus_gminus = build_template('MET', 'beta-strand', 'g-/g-')
IC_Met_Strand_gminus_trans = build_template('MET', 'beta-strand', 'g-/t')
IC_Met_Strand_gminus_gplus = build_template('MET', 'beta-strand', 'g-/g+')
IC_Met_Strand_trans_gminus = build_template('MET', 'beta-strand', 't/g-')
IC_Met_Strand_trans_trans = build_template('MET', 'beta-strand', 't/t')
IC_Met_Strand_trans_gplus = build_template('MET', 'beta-strand', 't/g+')
IC_Met_Strand_gplus_gminus = build_template('MET', 'beta-strand', 'g+/g-')
IC_Met_Strand_gplus_trans = build_template('MET', 'beta-strand', 'g+/t')
IC_Met_Strand_gplus_gplus = build_template('MET', 'beta-strand', 'g+/g+')
IC_Met_StrandParallel = build_template('MET', 'parallel-beta-strand', 'canonical')
IC_Met_StrandParallel_gminus_gminus = build_template('MET', 'parallel-beta-strand', 'g-/g-')
IC_Met_StrandParallel_gminus_trans = build_template('MET', 'parallel-beta-strand', 'g-/t')
IC_Met_StrandParallel_gminus_gplus = build_template('MET', 'parallel-beta-strand', 'g-/g+')
IC_Met_StrandParallel_trans_gminus = build_template('MET', 'parallel-beta-strand', 't/g-')
IC_Met_StrandParallel_trans_trans = build_template('MET', 'parallel-beta-strand', 't/t')
IC_Met_StrandParallel_trans_gplus = build_template('MET', 'parallel-beta-strand', 't/g+')
IC_Met_StrandParallel_gplus_gminus = build_template('MET', 'parallel-beta-strand', 'g+/g-')
IC_Met_StrandParallel_gplus_trans = build_template('MET', 'parallel-beta-strand', 'g+/t')
IC_Met_StrandParallel_gplus_gplus = build_template('MET', 'parallel-beta-strand', 'g+/g+')
IC_Met_StrandAntiParallel = build_template('MET', 'antiparallel-beta-strand', 'canonical')
IC_Met_StrandAntiParallel_gminus_gminus = build_template('MET', 'antiparallel-beta-strand', 'g-/g-')
IC_Met_StrandAntiParallel_gminus_trans = build_template('MET', 'antiparallel-beta-strand', 'g-/t')
IC_Met_StrandAntiParallel_gminus_gplus = build_template('MET', 'antiparallel-beta-strand', 'g-/g+')
IC_Met_StrandAntiParallel_trans_gminus = build_template('MET', 'antiparallel-beta-strand', 't/g-')
IC_Met_StrandAntiParallel_trans_trans = build_template('MET', 'antiparallel-beta-strand', 't/t')
IC_Met_StrandAntiParallel_trans_gplus = build_template('MET', 'antiparallel-beta-strand', 't/g+')
IC_Met_StrandAntiParallel_gplus_gminus = build_template('MET', 'antiparallel-beta-strand', 'g+/g-')
IC_Met_StrandAntiParallel_gplus_trans = build_template('MET', 'antiparallel-beta-strand', 'g+/t')
IC_Met_StrandAntiParallel_gplus_gplus = build_template('MET', 'antiparallel-beta-strand', 'g+/g+')
IC_Met_Bridge = build_template('MET', 'beta-bridge', 'canonical')
IC_Met_Bridge_gminus_gminus = build_template('MET', 'beta-bridge', 'g-/g-')
IC_Met_Bridge_gminus_trans = build_template('MET', 'beta-bridge', 'g-/t')
IC_Met_Bridge_gminus_gplus = build_template('MET', 'beta-bridge', 'g-/g+')
IC_Met_Bridge_trans_gminus = build_template('MET', 'beta-bridge', 't/g-')
IC_Met_Bridge_trans_trans = build_template('MET', 'beta-bridge', 't/t')
IC_Met_Bridge_trans_gplus = build_template('MET', 'beta-bridge', 't/g+')
IC_Met_Bridge_gplus_gminus = build_template('MET', 'beta-bridge', 'g+/g-')
IC_Met_Bridge_gplus_trans = build_template('MET', 'beta-bridge', 'g+/t')
IC_Met_Bridge_gplus_gplus = build_template('MET', 'beta-bridge', 'g+/g+')
IC_Met_Turn = build_template('MET', 'turn', 'canonical')
IC_Met_Turn_gminus_gminus = build_template('MET', 'turn', 'g-/g-')
IC_Met_Turn_gminus_trans = build_template('MET', 'turn', 'g-/t')
IC_Met_Turn_gminus_gplus = build_template('MET', 'turn', 'g-/g+')
IC_Met_Turn_trans_gminus = build_template('MET', 'turn', 't/g-')
IC_Met_Turn_trans_trans = build_template('MET', 'turn', 't/t')
IC_Met_Turn_trans_gplus = build_template('MET', 'turn', 't/g+')
IC_Met_Turn_gplus_gminus = build_template('MET', 'turn', 'g+/g-')
IC_Met_Turn_gplus_trans = build_template('MET', 'turn', 'g+/t')
IC_Met_Turn_gplus_gplus = build_template('MET', 'turn', 'g+/g+')
IC_Met_Bend = build_template('MET', 'bend', 'canonical')
IC_Met_Bend_gminus_gminus = build_template('MET', 'bend', 'g-/g-')
IC_Met_Bend_gminus_trans = build_template('MET', 'bend', 'g-/t')
IC_Met_Bend_gminus_gplus = build_template('MET', 'bend', 'g-/g+')
IC_Met_Bend_trans_gminus = build_template('MET', 'bend', 't/g-')
IC_Met_Bend_trans_trans = build_template('MET', 'bend', 't/t')
IC_Met_Bend_trans_gplus = build_template('MET', 'bend', 't/g+')
IC_Met_Bend_gplus_gminus = build_template('MET', 'bend', 'g+/g-')
IC_Met_Bend_gplus_trans = build_template('MET', 'bend', 'g+/t')
IC_Met_Bend_gplus_gplus = build_template('MET', 'bend', 'g+/g+')
IC_Met_Coil = build_template('MET', 'coil', 'canonical')
IC_Met_Coil_gminus_gminus = build_template('MET', 'coil', 'g-/g-')
IC_Met_Coil_gminus_trans = build_template('MET', 'coil', 'g-/t')
IC_Met_Coil_gminus_gplus = build_template('MET', 'coil', 'g-/g+')
IC_Met_Coil_trans_gminus = build_template('MET', 'coil', 't/g-')
IC_Met_Coil_trans_trans = build_template('MET', 'coil', 't/t')
IC_Met_Coil_trans_gplus = build_template('MET', 'coil', 't/g+')
IC_Met_Coil_gplus_gminus = build_template('MET', 'coil', 'g+/g-')
IC_Met_Coil_gplus_trans = build_template('MET', 'coil', 'g+/t')
IC_Met_Coil_gplus_gplus = build_template('MET', 'coil', 'g+/g+')
IC_Met_CisPeptide = build_template('MET', 'cis-peptide-bond', 'canonical')
IC_Met_CisPeptide_gminus_gminus = build_template('MET', 'cis-peptide-bond', 'g-/g-')
IC_Met_CisPeptide_gminus_trans = build_template('MET', 'cis-peptide-bond', 'g-/t')
IC_Met_CisPeptide_gminus_gplus = build_template('MET', 'cis-peptide-bond', 'g-/g+')
IC_Met_CisPeptide_trans_gminus = build_template('MET', 'cis-peptide-bond', 't/g-')
IC_Met_CisPeptide_trans_trans = build_template('MET', 'cis-peptide-bond', 't/t')
IC_Met_CisPeptide_trans_gplus = build_template('MET', 'cis-peptide-bond', 't/g+')
IC_Met_CisPeptide_gplus_gminus = build_template('MET', 'cis-peptide-bond', 'g+/g-')
IC_Met_CisPeptide_gplus_trans = build_template('MET', 'cis-peptide-bond', 'g+/t')
IC_Met_CisPeptide_gplus_gplus = build_template('MET', 'cis-peptide-bond', 'g+/g+')

__all__ = [
    'IC_Met_Bend',
    'IC_Met_Bend_gminus_gminus',
    'IC_Met_Bend_gminus_gplus',
    'IC_Met_Bend_gminus_trans',
    'IC_Met_Bend_gplus_gminus',
    'IC_Met_Bend_gplus_gplus',
    'IC_Met_Bend_gplus_trans',
    'IC_Met_Bend_trans_gminus',
    'IC_Met_Bend_trans_gplus',
    'IC_Met_Bend_trans_trans',
    'IC_Met_Bridge',
    'IC_Met_Bridge_gminus_gminus',
    'IC_Met_Bridge_gminus_gplus',
    'IC_Met_Bridge_gminus_trans',
    'IC_Met_Bridge_gplus_gminus',
    'IC_Met_Bridge_gplus_gplus',
    'IC_Met_Bridge_gplus_trans',
    'IC_Met_Bridge_trans_gminus',
    'IC_Met_Bridge_trans_gplus',
    'IC_Met_Bridge_trans_trans',
    'IC_Met_CisPeptide',
    'IC_Met_CisPeptide_gminus_gminus',
    'IC_Met_CisPeptide_gminus_gplus',
    'IC_Met_CisPeptide_gminus_trans',
    'IC_Met_CisPeptide_gplus_gminus',
    'IC_Met_CisPeptide_gplus_gplus',
    'IC_Met_CisPeptide_gplus_trans',
    'IC_Met_CisPeptide_trans_gminus',
    'IC_Met_CisPeptide_trans_gplus',
    'IC_Met_CisPeptide_trans_trans',
    'IC_Met_Coil',
    'IC_Met_Coil_gminus_gminus',
    'IC_Met_Coil_gminus_gplus',
    'IC_Met_Coil_gminus_trans',
    'IC_Met_Coil_gplus_gminus',
    'IC_Met_Coil_gplus_gplus',
    'IC_Met_Coil_gplus_trans',
    'IC_Met_Coil_trans_gminus',
    'IC_Met_Coil_trans_gplus',
    'IC_Met_Coil_trans_trans',
    'IC_Met_Helix310',
    'IC_Met_Helix310_gminus_gminus',
    'IC_Met_Helix310_gminus_gplus',
    'IC_Met_Helix310_gminus_trans',
    'IC_Met_Helix310_gplus_gminus',
    'IC_Met_Helix310_gplus_gplus',
    'IC_Met_Helix310_gplus_trans',
    'IC_Met_Helix310_trans_gminus',
    'IC_Met_Helix310_trans_gplus',
    'IC_Met_Helix310_trans_trans',
    'IC_Met_HelixAlpha',
    'IC_Met_HelixAlpha_gminus_gminus',
    'IC_Met_HelixAlpha_gminus_gplus',
    'IC_Met_HelixAlpha_gminus_trans',
    'IC_Met_HelixAlpha_gplus_gminus',
    'IC_Met_HelixAlpha_gplus_gplus',
    'IC_Met_HelixAlpha_gplus_trans',
    'IC_Met_HelixAlpha_trans_gminus',
    'IC_Met_HelixAlpha_trans_gplus',
    'IC_Met_HelixAlpha_trans_trans',
    'IC_Met_HelixPPII',
    'IC_Met_HelixPPII_gminus_gminus',
    'IC_Met_HelixPPII_gminus_gplus',
    'IC_Met_HelixPPII_gminus_trans',
    'IC_Met_HelixPPII_gplus_gminus',
    'IC_Met_HelixPPII_gplus_gplus',
    'IC_Met_HelixPPII_gplus_trans',
    'IC_Met_HelixPPII_trans_gminus',
    'IC_Met_HelixPPII_trans_gplus',
    'IC_Met_HelixPPII_trans_trans',
    'IC_Met_HelixPi',
    'IC_Met_HelixPi_gminus_gminus',
    'IC_Met_HelixPi_gminus_gplus',
    'IC_Met_HelixPi_gminus_trans',
    'IC_Met_HelixPi_gplus_gminus',
    'IC_Met_HelixPi_gplus_gplus',
    'IC_Met_HelixPi_gplus_trans',
    'IC_Met_HelixPi_trans_gminus',
    'IC_Met_HelixPi_trans_gplus',
    'IC_Met_HelixPi_trans_trans',
    'IC_Met_Strand',
    'IC_Met_StrandAntiParallel',
    'IC_Met_StrandAntiParallel_gminus_gminus',
    'IC_Met_StrandAntiParallel_gminus_gplus',
    'IC_Met_StrandAntiParallel_gminus_trans',
    'IC_Met_StrandAntiParallel_gplus_gminus',
    'IC_Met_StrandAntiParallel_gplus_gplus',
    'IC_Met_StrandAntiParallel_gplus_trans',
    'IC_Met_StrandAntiParallel_trans_gminus',
    'IC_Met_StrandAntiParallel_trans_gplus',
    'IC_Met_StrandAntiParallel_trans_trans',
    'IC_Met_StrandParallel',
    'IC_Met_StrandParallel_gminus_gminus',
    'IC_Met_StrandParallel_gminus_gplus',
    'IC_Met_StrandParallel_gminus_trans',
    'IC_Met_StrandParallel_gplus_gminus',
    'IC_Met_StrandParallel_gplus_gplus',
    'IC_Met_StrandParallel_gplus_trans',
    'IC_Met_StrandParallel_trans_gminus',
    'IC_Met_StrandParallel_trans_gplus',
    'IC_Met_StrandParallel_trans_trans',
    'IC_Met_Strand_gminus_gminus',
    'IC_Met_Strand_gminus_gplus',
    'IC_Met_Strand_gminus_trans',
    'IC_Met_Strand_gplus_gminus',
    'IC_Met_Strand_gplus_gplus',
    'IC_Met_Strand_gplus_trans',
    'IC_Met_Strand_trans_gminus',
    'IC_Met_Strand_trans_gplus',
    'IC_Met_Strand_trans_trans',
    'IC_Met_Turn',
    'IC_Met_Turn_gminus_gminus',
    'IC_Met_Turn_gminus_gplus',
    'IC_Met_Turn_gminus_trans',
    'IC_Met_Turn_gplus_gminus',
    'IC_Met_Turn_gplus_gplus',
    'IC_Met_Turn_gplus_trans',
    'IC_Met_Turn_trans_gminus',
    'IC_Met_Turn_trans_gplus',
    'IC_Met_Turn_trans_trans',
]
