# -*- coding: utf-8 -*-
"""Ideal {secondary-structure x rotamer} {names[aa]} {var} ``InternalCoord`` templates.

Each template is the ideal heavy-atom topology for {aa} anchored at
{{N, CA, C}}, with the side chain grown to ``rotamer`` and the
``phi``/``psi``/``omega`` of ``ss`` carried on the instance.  See
:mod:`biorazer.database.molecule.icoor.protein.template` and
:mod:`._builder` for the build rule and the per-conformer caveat.
"""

from biorazer.database.molecule.icoor.protein.template._builder import build_template

RESN = "HIS"

# === GENERATED IC_* named templates - do not edit below ===
# Regenerate with: /opt/envs/BioRazer/bin/python scripts/generate_internal_coord_template_named.py

IC_His_HelixAlpha = build_template('HIS', 'alpha-helix', 'canonical')
IC_His_HelixAlpha_gminus_gminus = build_template('HIS', 'alpha-helix', 'g-/g-')
IC_His_HelixAlpha_gminus_trans = build_template('HIS', 'alpha-helix', 'g-/t')
IC_His_HelixAlpha_gminus_gplus = build_template('HIS', 'alpha-helix', 'g-/g+')
IC_His_HelixAlpha_trans_gminus = build_template('HIS', 'alpha-helix', 't/g-')
IC_His_HelixAlpha_trans_trans = build_template('HIS', 'alpha-helix', 't/t')
IC_His_HelixAlpha_trans_gplus = build_template('HIS', 'alpha-helix', 't/g+')
IC_His_HelixAlpha_gplus_gminus = build_template('HIS', 'alpha-helix', 'g+/g-')
IC_His_HelixAlpha_gplus_trans = build_template('HIS', 'alpha-helix', 'g+/t')
IC_His_HelixAlpha_gplus_gplus = build_template('HIS', 'alpha-helix', 'g+/g+')
IC_His_Helix310 = build_template('HIS', '3-10-helix', 'canonical')
IC_His_Helix310_gminus_gminus = build_template('HIS', '3-10-helix', 'g-/g-')
IC_His_Helix310_gminus_trans = build_template('HIS', '3-10-helix', 'g-/t')
IC_His_Helix310_gminus_gplus = build_template('HIS', '3-10-helix', 'g-/g+')
IC_His_Helix310_trans_gminus = build_template('HIS', '3-10-helix', 't/g-')
IC_His_Helix310_trans_trans = build_template('HIS', '3-10-helix', 't/t')
IC_His_Helix310_trans_gplus = build_template('HIS', '3-10-helix', 't/g+')
IC_His_Helix310_gplus_gminus = build_template('HIS', '3-10-helix', 'g+/g-')
IC_His_Helix310_gplus_trans = build_template('HIS', '3-10-helix', 'g+/t')
IC_His_Helix310_gplus_gplus = build_template('HIS', '3-10-helix', 'g+/g+')
IC_His_HelixPi = build_template('HIS', 'pi-helix', 'canonical')
IC_His_HelixPi_gminus_gminus = build_template('HIS', 'pi-helix', 'g-/g-')
IC_His_HelixPi_gminus_trans = build_template('HIS', 'pi-helix', 'g-/t')
IC_His_HelixPi_gminus_gplus = build_template('HIS', 'pi-helix', 'g-/g+')
IC_His_HelixPi_trans_gminus = build_template('HIS', 'pi-helix', 't/g-')
IC_His_HelixPi_trans_trans = build_template('HIS', 'pi-helix', 't/t')
IC_His_HelixPi_trans_gplus = build_template('HIS', 'pi-helix', 't/g+')
IC_His_HelixPi_gplus_gminus = build_template('HIS', 'pi-helix', 'g+/g-')
IC_His_HelixPi_gplus_trans = build_template('HIS', 'pi-helix', 'g+/t')
IC_His_HelixPi_gplus_gplus = build_template('HIS', 'pi-helix', 'g+/g+')
IC_His_HelixPPII = build_template('HIS', 'polyproline-II', 'canonical')
IC_His_HelixPPII_gminus_gminus = build_template('HIS', 'polyproline-II', 'g-/g-')
IC_His_HelixPPII_gminus_trans = build_template('HIS', 'polyproline-II', 'g-/t')
IC_His_HelixPPII_gminus_gplus = build_template('HIS', 'polyproline-II', 'g-/g+')
IC_His_HelixPPII_trans_gminus = build_template('HIS', 'polyproline-II', 't/g-')
IC_His_HelixPPII_trans_trans = build_template('HIS', 'polyproline-II', 't/t')
IC_His_HelixPPII_trans_gplus = build_template('HIS', 'polyproline-II', 't/g+')
IC_His_HelixPPII_gplus_gminus = build_template('HIS', 'polyproline-II', 'g+/g-')
IC_His_HelixPPII_gplus_trans = build_template('HIS', 'polyproline-II', 'g+/t')
IC_His_HelixPPII_gplus_gplus = build_template('HIS', 'polyproline-II', 'g+/g+')
IC_His_Strand = build_template('HIS', 'beta-strand', 'canonical')
IC_His_Strand_gminus_gminus = build_template('HIS', 'beta-strand', 'g-/g-')
IC_His_Strand_gminus_trans = build_template('HIS', 'beta-strand', 'g-/t')
IC_His_Strand_gminus_gplus = build_template('HIS', 'beta-strand', 'g-/g+')
IC_His_Strand_trans_gminus = build_template('HIS', 'beta-strand', 't/g-')
IC_His_Strand_trans_trans = build_template('HIS', 'beta-strand', 't/t')
IC_His_Strand_trans_gplus = build_template('HIS', 'beta-strand', 't/g+')
IC_His_Strand_gplus_gminus = build_template('HIS', 'beta-strand', 'g+/g-')
IC_His_Strand_gplus_trans = build_template('HIS', 'beta-strand', 'g+/t')
IC_His_Strand_gplus_gplus = build_template('HIS', 'beta-strand', 'g+/g+')
IC_His_StrandParallel = build_template('HIS', 'parallel-beta-strand', 'canonical')
IC_His_StrandParallel_gminus_gminus = build_template('HIS', 'parallel-beta-strand', 'g-/g-')
IC_His_StrandParallel_gminus_trans = build_template('HIS', 'parallel-beta-strand', 'g-/t')
IC_His_StrandParallel_gminus_gplus = build_template('HIS', 'parallel-beta-strand', 'g-/g+')
IC_His_StrandParallel_trans_gminus = build_template('HIS', 'parallel-beta-strand', 't/g-')
IC_His_StrandParallel_trans_trans = build_template('HIS', 'parallel-beta-strand', 't/t')
IC_His_StrandParallel_trans_gplus = build_template('HIS', 'parallel-beta-strand', 't/g+')
IC_His_StrandParallel_gplus_gminus = build_template('HIS', 'parallel-beta-strand', 'g+/g-')
IC_His_StrandParallel_gplus_trans = build_template('HIS', 'parallel-beta-strand', 'g+/t')
IC_His_StrandParallel_gplus_gplus = build_template('HIS', 'parallel-beta-strand', 'g+/g+')
IC_His_StrandAntiParallel = build_template('HIS', 'antiparallel-beta-strand', 'canonical')
IC_His_StrandAntiParallel_gminus_gminus = build_template('HIS', 'antiparallel-beta-strand', 'g-/g-')
IC_His_StrandAntiParallel_gminus_trans = build_template('HIS', 'antiparallel-beta-strand', 'g-/t')
IC_His_StrandAntiParallel_gminus_gplus = build_template('HIS', 'antiparallel-beta-strand', 'g-/g+')
IC_His_StrandAntiParallel_trans_gminus = build_template('HIS', 'antiparallel-beta-strand', 't/g-')
IC_His_StrandAntiParallel_trans_trans = build_template('HIS', 'antiparallel-beta-strand', 't/t')
IC_His_StrandAntiParallel_trans_gplus = build_template('HIS', 'antiparallel-beta-strand', 't/g+')
IC_His_StrandAntiParallel_gplus_gminus = build_template('HIS', 'antiparallel-beta-strand', 'g+/g-')
IC_His_StrandAntiParallel_gplus_trans = build_template('HIS', 'antiparallel-beta-strand', 'g+/t')
IC_His_StrandAntiParallel_gplus_gplus = build_template('HIS', 'antiparallel-beta-strand', 'g+/g+')
IC_His_Bridge = build_template('HIS', 'beta-bridge', 'canonical')
IC_His_Bridge_gminus_gminus = build_template('HIS', 'beta-bridge', 'g-/g-')
IC_His_Bridge_gminus_trans = build_template('HIS', 'beta-bridge', 'g-/t')
IC_His_Bridge_gminus_gplus = build_template('HIS', 'beta-bridge', 'g-/g+')
IC_His_Bridge_trans_gminus = build_template('HIS', 'beta-bridge', 't/g-')
IC_His_Bridge_trans_trans = build_template('HIS', 'beta-bridge', 't/t')
IC_His_Bridge_trans_gplus = build_template('HIS', 'beta-bridge', 't/g+')
IC_His_Bridge_gplus_gminus = build_template('HIS', 'beta-bridge', 'g+/g-')
IC_His_Bridge_gplus_trans = build_template('HIS', 'beta-bridge', 'g+/t')
IC_His_Bridge_gplus_gplus = build_template('HIS', 'beta-bridge', 'g+/g+')
IC_His_Turn = build_template('HIS', 'turn', 'canonical')
IC_His_Turn_gminus_gminus = build_template('HIS', 'turn', 'g-/g-')
IC_His_Turn_gminus_trans = build_template('HIS', 'turn', 'g-/t')
IC_His_Turn_gminus_gplus = build_template('HIS', 'turn', 'g-/g+')
IC_His_Turn_trans_gminus = build_template('HIS', 'turn', 't/g-')
IC_His_Turn_trans_trans = build_template('HIS', 'turn', 't/t')
IC_His_Turn_trans_gplus = build_template('HIS', 'turn', 't/g+')
IC_His_Turn_gplus_gminus = build_template('HIS', 'turn', 'g+/g-')
IC_His_Turn_gplus_trans = build_template('HIS', 'turn', 'g+/t')
IC_His_Turn_gplus_gplus = build_template('HIS', 'turn', 'g+/g+')
IC_His_Bend = build_template('HIS', 'bend', 'canonical')
IC_His_Bend_gminus_gminus = build_template('HIS', 'bend', 'g-/g-')
IC_His_Bend_gminus_trans = build_template('HIS', 'bend', 'g-/t')
IC_His_Bend_gminus_gplus = build_template('HIS', 'bend', 'g-/g+')
IC_His_Bend_trans_gminus = build_template('HIS', 'bend', 't/g-')
IC_His_Bend_trans_trans = build_template('HIS', 'bend', 't/t')
IC_His_Bend_trans_gplus = build_template('HIS', 'bend', 't/g+')
IC_His_Bend_gplus_gminus = build_template('HIS', 'bend', 'g+/g-')
IC_His_Bend_gplus_trans = build_template('HIS', 'bend', 'g+/t')
IC_His_Bend_gplus_gplus = build_template('HIS', 'bend', 'g+/g+')
IC_His_Coil = build_template('HIS', 'coil', 'canonical')
IC_His_Coil_gminus_gminus = build_template('HIS', 'coil', 'g-/g-')
IC_His_Coil_gminus_trans = build_template('HIS', 'coil', 'g-/t')
IC_His_Coil_gminus_gplus = build_template('HIS', 'coil', 'g-/g+')
IC_His_Coil_trans_gminus = build_template('HIS', 'coil', 't/g-')
IC_His_Coil_trans_trans = build_template('HIS', 'coil', 't/t')
IC_His_Coil_trans_gplus = build_template('HIS', 'coil', 't/g+')
IC_His_Coil_gplus_gminus = build_template('HIS', 'coil', 'g+/g-')
IC_His_Coil_gplus_trans = build_template('HIS', 'coil', 'g+/t')
IC_His_Coil_gplus_gplus = build_template('HIS', 'coil', 'g+/g+')
IC_His_CisPeptide = build_template('HIS', 'cis-peptide-bond', 'canonical')
IC_His_CisPeptide_gminus_gminus = build_template('HIS', 'cis-peptide-bond', 'g-/g-')
IC_His_CisPeptide_gminus_trans = build_template('HIS', 'cis-peptide-bond', 'g-/t')
IC_His_CisPeptide_gminus_gplus = build_template('HIS', 'cis-peptide-bond', 'g-/g+')
IC_His_CisPeptide_trans_gminus = build_template('HIS', 'cis-peptide-bond', 't/g-')
IC_His_CisPeptide_trans_trans = build_template('HIS', 'cis-peptide-bond', 't/t')
IC_His_CisPeptide_trans_gplus = build_template('HIS', 'cis-peptide-bond', 't/g+')
IC_His_CisPeptide_gplus_gminus = build_template('HIS', 'cis-peptide-bond', 'g+/g-')
IC_His_CisPeptide_gplus_trans = build_template('HIS', 'cis-peptide-bond', 'g+/t')
IC_His_CisPeptide_gplus_gplus = build_template('HIS', 'cis-peptide-bond', 'g+/g+')

__all__ = [
    'IC_His_Bend',
    'IC_His_Bend_gminus_gminus',
    'IC_His_Bend_gminus_gplus',
    'IC_His_Bend_gminus_trans',
    'IC_His_Bend_gplus_gminus',
    'IC_His_Bend_gplus_gplus',
    'IC_His_Bend_gplus_trans',
    'IC_His_Bend_trans_gminus',
    'IC_His_Bend_trans_gplus',
    'IC_His_Bend_trans_trans',
    'IC_His_Bridge',
    'IC_His_Bridge_gminus_gminus',
    'IC_His_Bridge_gminus_gplus',
    'IC_His_Bridge_gminus_trans',
    'IC_His_Bridge_gplus_gminus',
    'IC_His_Bridge_gplus_gplus',
    'IC_His_Bridge_gplus_trans',
    'IC_His_Bridge_trans_gminus',
    'IC_His_Bridge_trans_gplus',
    'IC_His_Bridge_trans_trans',
    'IC_His_CisPeptide',
    'IC_His_CisPeptide_gminus_gminus',
    'IC_His_CisPeptide_gminus_gplus',
    'IC_His_CisPeptide_gminus_trans',
    'IC_His_CisPeptide_gplus_gminus',
    'IC_His_CisPeptide_gplus_gplus',
    'IC_His_CisPeptide_gplus_trans',
    'IC_His_CisPeptide_trans_gminus',
    'IC_His_CisPeptide_trans_gplus',
    'IC_His_CisPeptide_trans_trans',
    'IC_His_Coil',
    'IC_His_Coil_gminus_gminus',
    'IC_His_Coil_gminus_gplus',
    'IC_His_Coil_gminus_trans',
    'IC_His_Coil_gplus_gminus',
    'IC_His_Coil_gplus_gplus',
    'IC_His_Coil_gplus_trans',
    'IC_His_Coil_trans_gminus',
    'IC_His_Coil_trans_gplus',
    'IC_His_Coil_trans_trans',
    'IC_His_Helix310',
    'IC_His_Helix310_gminus_gminus',
    'IC_His_Helix310_gminus_gplus',
    'IC_His_Helix310_gminus_trans',
    'IC_His_Helix310_gplus_gminus',
    'IC_His_Helix310_gplus_gplus',
    'IC_His_Helix310_gplus_trans',
    'IC_His_Helix310_trans_gminus',
    'IC_His_Helix310_trans_gplus',
    'IC_His_Helix310_trans_trans',
    'IC_His_HelixAlpha',
    'IC_His_HelixAlpha_gminus_gminus',
    'IC_His_HelixAlpha_gminus_gplus',
    'IC_His_HelixAlpha_gminus_trans',
    'IC_His_HelixAlpha_gplus_gminus',
    'IC_His_HelixAlpha_gplus_gplus',
    'IC_His_HelixAlpha_gplus_trans',
    'IC_His_HelixAlpha_trans_gminus',
    'IC_His_HelixAlpha_trans_gplus',
    'IC_His_HelixAlpha_trans_trans',
    'IC_His_HelixPPII',
    'IC_His_HelixPPII_gminus_gminus',
    'IC_His_HelixPPII_gminus_gplus',
    'IC_His_HelixPPII_gminus_trans',
    'IC_His_HelixPPII_gplus_gminus',
    'IC_His_HelixPPII_gplus_gplus',
    'IC_His_HelixPPII_gplus_trans',
    'IC_His_HelixPPII_trans_gminus',
    'IC_His_HelixPPII_trans_gplus',
    'IC_His_HelixPPII_trans_trans',
    'IC_His_HelixPi',
    'IC_His_HelixPi_gminus_gminus',
    'IC_His_HelixPi_gminus_gplus',
    'IC_His_HelixPi_gminus_trans',
    'IC_His_HelixPi_gplus_gminus',
    'IC_His_HelixPi_gplus_gplus',
    'IC_His_HelixPi_gplus_trans',
    'IC_His_HelixPi_trans_gminus',
    'IC_His_HelixPi_trans_gplus',
    'IC_His_HelixPi_trans_trans',
    'IC_His_Strand',
    'IC_His_StrandAntiParallel',
    'IC_His_StrandAntiParallel_gminus_gminus',
    'IC_His_StrandAntiParallel_gminus_gplus',
    'IC_His_StrandAntiParallel_gminus_trans',
    'IC_His_StrandAntiParallel_gplus_gminus',
    'IC_His_StrandAntiParallel_gplus_gplus',
    'IC_His_StrandAntiParallel_gplus_trans',
    'IC_His_StrandAntiParallel_trans_gminus',
    'IC_His_StrandAntiParallel_trans_gplus',
    'IC_His_StrandAntiParallel_trans_trans',
    'IC_His_StrandParallel',
    'IC_His_StrandParallel_gminus_gminus',
    'IC_His_StrandParallel_gminus_gplus',
    'IC_His_StrandParallel_gminus_trans',
    'IC_His_StrandParallel_gplus_gminus',
    'IC_His_StrandParallel_gplus_gplus',
    'IC_His_StrandParallel_gplus_trans',
    'IC_His_StrandParallel_trans_gminus',
    'IC_His_StrandParallel_trans_gplus',
    'IC_His_StrandParallel_trans_trans',
    'IC_His_Strand_gminus_gminus',
    'IC_His_Strand_gminus_gplus',
    'IC_His_Strand_gminus_trans',
    'IC_His_Strand_gplus_gminus',
    'IC_His_Strand_gplus_gplus',
    'IC_His_Strand_gplus_trans',
    'IC_His_Strand_trans_gminus',
    'IC_His_Strand_trans_gplus',
    'IC_His_Strand_trans_trans',
    'IC_His_Turn',
    'IC_His_Turn_gminus_gminus',
    'IC_His_Turn_gminus_gplus',
    'IC_His_Turn_gminus_trans',
    'IC_His_Turn_gplus_gminus',
    'IC_His_Turn_gplus_gplus',
    'IC_His_Turn_gplus_trans',
    'IC_His_Turn_trans_gminus',
    'IC_His_Turn_trans_gplus',
    'IC_His_Turn_trans_trans',
]
