# -*- coding: utf-8 -*-
"""Ideal {secondary-structure x rotamer} {names[aa]} {var} ``InternalCoord`` templates.

Each template is the ideal heavy-atom topology for {aa} anchored at
{{N, CA, C}}, with the side chain grown to ``rotamer`` and the
``phi``/``psi``/``omega`` of ``ss`` carried on the instance.  See
:mod:`biorazer.database.molecule.icoor.protein.template` and
:mod:`._builder` for the build rule and the per-conformer caveat.
"""

from biorazer.database.molecule.icoor.protein.template._builder import build_template

RESN = "GLU"

# === GENERATED IC_* named templates - do not edit below ===
# Regenerate with: /opt/envs/BioRazer/bin/python scripts/generate_internal_coord_template_named.py

IC_Glu_HelixAlpha = build_template('GLU', 'alpha-helix', 'canonical')
IC_Glu_HelixAlpha_gminus_gminus = build_template('GLU', 'alpha-helix', 'g-/g-')
IC_Glu_HelixAlpha_gminus_trans = build_template('GLU', 'alpha-helix', 'g-/t')
IC_Glu_HelixAlpha_gminus_gplus = build_template('GLU', 'alpha-helix', 'g-/g+')
IC_Glu_HelixAlpha_trans_gminus = build_template('GLU', 'alpha-helix', 't/g-')
IC_Glu_HelixAlpha_trans_trans = build_template('GLU', 'alpha-helix', 't/t')
IC_Glu_HelixAlpha_trans_gplus = build_template('GLU', 'alpha-helix', 't/g+')
IC_Glu_HelixAlpha_gplus_gminus = build_template('GLU', 'alpha-helix', 'g+/g-')
IC_Glu_HelixAlpha_gplus_trans = build_template('GLU', 'alpha-helix', 'g+/t')
IC_Glu_HelixAlpha_gplus_gplus = build_template('GLU', 'alpha-helix', 'g+/g+')
IC_Glu_Helix310 = build_template('GLU', '3-10-helix', 'canonical')
IC_Glu_Helix310_gminus_gminus = build_template('GLU', '3-10-helix', 'g-/g-')
IC_Glu_Helix310_gminus_trans = build_template('GLU', '3-10-helix', 'g-/t')
IC_Glu_Helix310_gminus_gplus = build_template('GLU', '3-10-helix', 'g-/g+')
IC_Glu_Helix310_trans_gminus = build_template('GLU', '3-10-helix', 't/g-')
IC_Glu_Helix310_trans_trans = build_template('GLU', '3-10-helix', 't/t')
IC_Glu_Helix310_trans_gplus = build_template('GLU', '3-10-helix', 't/g+')
IC_Glu_Helix310_gplus_gminus = build_template('GLU', '3-10-helix', 'g+/g-')
IC_Glu_Helix310_gplus_trans = build_template('GLU', '3-10-helix', 'g+/t')
IC_Glu_Helix310_gplus_gplus = build_template('GLU', '3-10-helix', 'g+/g+')
IC_Glu_HelixPi = build_template('GLU', 'pi-helix', 'canonical')
IC_Glu_HelixPi_gminus_gminus = build_template('GLU', 'pi-helix', 'g-/g-')
IC_Glu_HelixPi_gminus_trans = build_template('GLU', 'pi-helix', 'g-/t')
IC_Glu_HelixPi_gminus_gplus = build_template('GLU', 'pi-helix', 'g-/g+')
IC_Glu_HelixPi_trans_gminus = build_template('GLU', 'pi-helix', 't/g-')
IC_Glu_HelixPi_trans_trans = build_template('GLU', 'pi-helix', 't/t')
IC_Glu_HelixPi_trans_gplus = build_template('GLU', 'pi-helix', 't/g+')
IC_Glu_HelixPi_gplus_gminus = build_template('GLU', 'pi-helix', 'g+/g-')
IC_Glu_HelixPi_gplus_trans = build_template('GLU', 'pi-helix', 'g+/t')
IC_Glu_HelixPi_gplus_gplus = build_template('GLU', 'pi-helix', 'g+/g+')
IC_Glu_HelixPPII = build_template('GLU', 'polyproline-II', 'canonical')
IC_Glu_HelixPPII_gminus_gminus = build_template('GLU', 'polyproline-II', 'g-/g-')
IC_Glu_HelixPPII_gminus_trans = build_template('GLU', 'polyproline-II', 'g-/t')
IC_Glu_HelixPPII_gminus_gplus = build_template('GLU', 'polyproline-II', 'g-/g+')
IC_Glu_HelixPPII_trans_gminus = build_template('GLU', 'polyproline-II', 't/g-')
IC_Glu_HelixPPII_trans_trans = build_template('GLU', 'polyproline-II', 't/t')
IC_Glu_HelixPPII_trans_gplus = build_template('GLU', 'polyproline-II', 't/g+')
IC_Glu_HelixPPII_gplus_gminus = build_template('GLU', 'polyproline-II', 'g+/g-')
IC_Glu_HelixPPII_gplus_trans = build_template('GLU', 'polyproline-II', 'g+/t')
IC_Glu_HelixPPII_gplus_gplus = build_template('GLU', 'polyproline-II', 'g+/g+')
IC_Glu_Strand = build_template('GLU', 'beta-strand', 'canonical')
IC_Glu_Strand_gminus_gminus = build_template('GLU', 'beta-strand', 'g-/g-')
IC_Glu_Strand_gminus_trans = build_template('GLU', 'beta-strand', 'g-/t')
IC_Glu_Strand_gminus_gplus = build_template('GLU', 'beta-strand', 'g-/g+')
IC_Glu_Strand_trans_gminus = build_template('GLU', 'beta-strand', 't/g-')
IC_Glu_Strand_trans_trans = build_template('GLU', 'beta-strand', 't/t')
IC_Glu_Strand_trans_gplus = build_template('GLU', 'beta-strand', 't/g+')
IC_Glu_Strand_gplus_gminus = build_template('GLU', 'beta-strand', 'g+/g-')
IC_Glu_Strand_gplus_trans = build_template('GLU', 'beta-strand', 'g+/t')
IC_Glu_Strand_gplus_gplus = build_template('GLU', 'beta-strand', 'g+/g+')
IC_Glu_StrandParallel = build_template('GLU', 'parallel-beta-strand', 'canonical')
IC_Glu_StrandParallel_gminus_gminus = build_template('GLU', 'parallel-beta-strand', 'g-/g-')
IC_Glu_StrandParallel_gminus_trans = build_template('GLU', 'parallel-beta-strand', 'g-/t')
IC_Glu_StrandParallel_gminus_gplus = build_template('GLU', 'parallel-beta-strand', 'g-/g+')
IC_Glu_StrandParallel_trans_gminus = build_template('GLU', 'parallel-beta-strand', 't/g-')
IC_Glu_StrandParallel_trans_trans = build_template('GLU', 'parallel-beta-strand', 't/t')
IC_Glu_StrandParallel_trans_gplus = build_template('GLU', 'parallel-beta-strand', 't/g+')
IC_Glu_StrandParallel_gplus_gminus = build_template('GLU', 'parallel-beta-strand', 'g+/g-')
IC_Glu_StrandParallel_gplus_trans = build_template('GLU', 'parallel-beta-strand', 'g+/t')
IC_Glu_StrandParallel_gplus_gplus = build_template('GLU', 'parallel-beta-strand', 'g+/g+')
IC_Glu_StrandAntiParallel = build_template('GLU', 'antiparallel-beta-strand', 'canonical')
IC_Glu_StrandAntiParallel_gminus_gminus = build_template('GLU', 'antiparallel-beta-strand', 'g-/g-')
IC_Glu_StrandAntiParallel_gminus_trans = build_template('GLU', 'antiparallel-beta-strand', 'g-/t')
IC_Glu_StrandAntiParallel_gminus_gplus = build_template('GLU', 'antiparallel-beta-strand', 'g-/g+')
IC_Glu_StrandAntiParallel_trans_gminus = build_template('GLU', 'antiparallel-beta-strand', 't/g-')
IC_Glu_StrandAntiParallel_trans_trans = build_template('GLU', 'antiparallel-beta-strand', 't/t')
IC_Glu_StrandAntiParallel_trans_gplus = build_template('GLU', 'antiparallel-beta-strand', 't/g+')
IC_Glu_StrandAntiParallel_gplus_gminus = build_template('GLU', 'antiparallel-beta-strand', 'g+/g-')
IC_Glu_StrandAntiParallel_gplus_trans = build_template('GLU', 'antiparallel-beta-strand', 'g+/t')
IC_Glu_StrandAntiParallel_gplus_gplus = build_template('GLU', 'antiparallel-beta-strand', 'g+/g+')
IC_Glu_Bridge = build_template('GLU', 'beta-bridge', 'canonical')
IC_Glu_Bridge_gminus_gminus = build_template('GLU', 'beta-bridge', 'g-/g-')
IC_Glu_Bridge_gminus_trans = build_template('GLU', 'beta-bridge', 'g-/t')
IC_Glu_Bridge_gminus_gplus = build_template('GLU', 'beta-bridge', 'g-/g+')
IC_Glu_Bridge_trans_gminus = build_template('GLU', 'beta-bridge', 't/g-')
IC_Glu_Bridge_trans_trans = build_template('GLU', 'beta-bridge', 't/t')
IC_Glu_Bridge_trans_gplus = build_template('GLU', 'beta-bridge', 't/g+')
IC_Glu_Bridge_gplus_gminus = build_template('GLU', 'beta-bridge', 'g+/g-')
IC_Glu_Bridge_gplus_trans = build_template('GLU', 'beta-bridge', 'g+/t')
IC_Glu_Bridge_gplus_gplus = build_template('GLU', 'beta-bridge', 'g+/g+')
IC_Glu_Turn = build_template('GLU', 'turn', 'canonical')
IC_Glu_Turn_gminus_gminus = build_template('GLU', 'turn', 'g-/g-')
IC_Glu_Turn_gminus_trans = build_template('GLU', 'turn', 'g-/t')
IC_Glu_Turn_gminus_gplus = build_template('GLU', 'turn', 'g-/g+')
IC_Glu_Turn_trans_gminus = build_template('GLU', 'turn', 't/g-')
IC_Glu_Turn_trans_trans = build_template('GLU', 'turn', 't/t')
IC_Glu_Turn_trans_gplus = build_template('GLU', 'turn', 't/g+')
IC_Glu_Turn_gplus_gminus = build_template('GLU', 'turn', 'g+/g-')
IC_Glu_Turn_gplus_trans = build_template('GLU', 'turn', 'g+/t')
IC_Glu_Turn_gplus_gplus = build_template('GLU', 'turn', 'g+/g+')
IC_Glu_Bend = build_template('GLU', 'bend', 'canonical')
IC_Glu_Bend_gminus_gminus = build_template('GLU', 'bend', 'g-/g-')
IC_Glu_Bend_gminus_trans = build_template('GLU', 'bend', 'g-/t')
IC_Glu_Bend_gminus_gplus = build_template('GLU', 'bend', 'g-/g+')
IC_Glu_Bend_trans_gminus = build_template('GLU', 'bend', 't/g-')
IC_Glu_Bend_trans_trans = build_template('GLU', 'bend', 't/t')
IC_Glu_Bend_trans_gplus = build_template('GLU', 'bend', 't/g+')
IC_Glu_Bend_gplus_gminus = build_template('GLU', 'bend', 'g+/g-')
IC_Glu_Bend_gplus_trans = build_template('GLU', 'bend', 'g+/t')
IC_Glu_Bend_gplus_gplus = build_template('GLU', 'bend', 'g+/g+')
IC_Glu_Coil = build_template('GLU', 'coil', 'canonical')
IC_Glu_Coil_gminus_gminus = build_template('GLU', 'coil', 'g-/g-')
IC_Glu_Coil_gminus_trans = build_template('GLU', 'coil', 'g-/t')
IC_Glu_Coil_gminus_gplus = build_template('GLU', 'coil', 'g-/g+')
IC_Glu_Coil_trans_gminus = build_template('GLU', 'coil', 't/g-')
IC_Glu_Coil_trans_trans = build_template('GLU', 'coil', 't/t')
IC_Glu_Coil_trans_gplus = build_template('GLU', 'coil', 't/g+')
IC_Glu_Coil_gplus_gminus = build_template('GLU', 'coil', 'g+/g-')
IC_Glu_Coil_gplus_trans = build_template('GLU', 'coil', 'g+/t')
IC_Glu_Coil_gplus_gplus = build_template('GLU', 'coil', 'g+/g+')
IC_Glu_CisPeptide = build_template('GLU', 'cis-peptide-bond', 'canonical')
IC_Glu_CisPeptide_gminus_gminus = build_template('GLU', 'cis-peptide-bond', 'g-/g-')
IC_Glu_CisPeptide_gminus_trans = build_template('GLU', 'cis-peptide-bond', 'g-/t')
IC_Glu_CisPeptide_gminus_gplus = build_template('GLU', 'cis-peptide-bond', 'g-/g+')
IC_Glu_CisPeptide_trans_gminus = build_template('GLU', 'cis-peptide-bond', 't/g-')
IC_Glu_CisPeptide_trans_trans = build_template('GLU', 'cis-peptide-bond', 't/t')
IC_Glu_CisPeptide_trans_gplus = build_template('GLU', 'cis-peptide-bond', 't/g+')
IC_Glu_CisPeptide_gplus_gminus = build_template('GLU', 'cis-peptide-bond', 'g+/g-')
IC_Glu_CisPeptide_gplus_trans = build_template('GLU', 'cis-peptide-bond', 'g+/t')
IC_Glu_CisPeptide_gplus_gplus = build_template('GLU', 'cis-peptide-bond', 'g+/g+')

__all__ = [
    'IC_Glu_Bend',
    'IC_Glu_Bend_gminus_gminus',
    'IC_Glu_Bend_gminus_gplus',
    'IC_Glu_Bend_gminus_trans',
    'IC_Glu_Bend_gplus_gminus',
    'IC_Glu_Bend_gplus_gplus',
    'IC_Glu_Bend_gplus_trans',
    'IC_Glu_Bend_trans_gminus',
    'IC_Glu_Bend_trans_gplus',
    'IC_Glu_Bend_trans_trans',
    'IC_Glu_Bridge',
    'IC_Glu_Bridge_gminus_gminus',
    'IC_Glu_Bridge_gminus_gplus',
    'IC_Glu_Bridge_gminus_trans',
    'IC_Glu_Bridge_gplus_gminus',
    'IC_Glu_Bridge_gplus_gplus',
    'IC_Glu_Bridge_gplus_trans',
    'IC_Glu_Bridge_trans_gminus',
    'IC_Glu_Bridge_trans_gplus',
    'IC_Glu_Bridge_trans_trans',
    'IC_Glu_CisPeptide',
    'IC_Glu_CisPeptide_gminus_gminus',
    'IC_Glu_CisPeptide_gminus_gplus',
    'IC_Glu_CisPeptide_gminus_trans',
    'IC_Glu_CisPeptide_gplus_gminus',
    'IC_Glu_CisPeptide_gplus_gplus',
    'IC_Glu_CisPeptide_gplus_trans',
    'IC_Glu_CisPeptide_trans_gminus',
    'IC_Glu_CisPeptide_trans_gplus',
    'IC_Glu_CisPeptide_trans_trans',
    'IC_Glu_Coil',
    'IC_Glu_Coil_gminus_gminus',
    'IC_Glu_Coil_gminus_gplus',
    'IC_Glu_Coil_gminus_trans',
    'IC_Glu_Coil_gplus_gminus',
    'IC_Glu_Coil_gplus_gplus',
    'IC_Glu_Coil_gplus_trans',
    'IC_Glu_Coil_trans_gminus',
    'IC_Glu_Coil_trans_gplus',
    'IC_Glu_Coil_trans_trans',
    'IC_Glu_Helix310',
    'IC_Glu_Helix310_gminus_gminus',
    'IC_Glu_Helix310_gminus_gplus',
    'IC_Glu_Helix310_gminus_trans',
    'IC_Glu_Helix310_gplus_gminus',
    'IC_Glu_Helix310_gplus_gplus',
    'IC_Glu_Helix310_gplus_trans',
    'IC_Glu_Helix310_trans_gminus',
    'IC_Glu_Helix310_trans_gplus',
    'IC_Glu_Helix310_trans_trans',
    'IC_Glu_HelixAlpha',
    'IC_Glu_HelixAlpha_gminus_gminus',
    'IC_Glu_HelixAlpha_gminus_gplus',
    'IC_Glu_HelixAlpha_gminus_trans',
    'IC_Glu_HelixAlpha_gplus_gminus',
    'IC_Glu_HelixAlpha_gplus_gplus',
    'IC_Glu_HelixAlpha_gplus_trans',
    'IC_Glu_HelixAlpha_trans_gminus',
    'IC_Glu_HelixAlpha_trans_gplus',
    'IC_Glu_HelixAlpha_trans_trans',
    'IC_Glu_HelixPPII',
    'IC_Glu_HelixPPII_gminus_gminus',
    'IC_Glu_HelixPPII_gminus_gplus',
    'IC_Glu_HelixPPII_gminus_trans',
    'IC_Glu_HelixPPII_gplus_gminus',
    'IC_Glu_HelixPPII_gplus_gplus',
    'IC_Glu_HelixPPII_gplus_trans',
    'IC_Glu_HelixPPII_trans_gminus',
    'IC_Glu_HelixPPII_trans_gplus',
    'IC_Glu_HelixPPII_trans_trans',
    'IC_Glu_HelixPi',
    'IC_Glu_HelixPi_gminus_gminus',
    'IC_Glu_HelixPi_gminus_gplus',
    'IC_Glu_HelixPi_gminus_trans',
    'IC_Glu_HelixPi_gplus_gminus',
    'IC_Glu_HelixPi_gplus_gplus',
    'IC_Glu_HelixPi_gplus_trans',
    'IC_Glu_HelixPi_trans_gminus',
    'IC_Glu_HelixPi_trans_gplus',
    'IC_Glu_HelixPi_trans_trans',
    'IC_Glu_Strand',
    'IC_Glu_StrandAntiParallel',
    'IC_Glu_StrandAntiParallel_gminus_gminus',
    'IC_Glu_StrandAntiParallel_gminus_gplus',
    'IC_Glu_StrandAntiParallel_gminus_trans',
    'IC_Glu_StrandAntiParallel_gplus_gminus',
    'IC_Glu_StrandAntiParallel_gplus_gplus',
    'IC_Glu_StrandAntiParallel_gplus_trans',
    'IC_Glu_StrandAntiParallel_trans_gminus',
    'IC_Glu_StrandAntiParallel_trans_gplus',
    'IC_Glu_StrandAntiParallel_trans_trans',
    'IC_Glu_StrandParallel',
    'IC_Glu_StrandParallel_gminus_gminus',
    'IC_Glu_StrandParallel_gminus_gplus',
    'IC_Glu_StrandParallel_gminus_trans',
    'IC_Glu_StrandParallel_gplus_gminus',
    'IC_Glu_StrandParallel_gplus_gplus',
    'IC_Glu_StrandParallel_gplus_trans',
    'IC_Glu_StrandParallel_trans_gminus',
    'IC_Glu_StrandParallel_trans_gplus',
    'IC_Glu_StrandParallel_trans_trans',
    'IC_Glu_Strand_gminus_gminus',
    'IC_Glu_Strand_gminus_gplus',
    'IC_Glu_Strand_gminus_trans',
    'IC_Glu_Strand_gplus_gminus',
    'IC_Glu_Strand_gplus_gplus',
    'IC_Glu_Strand_gplus_trans',
    'IC_Glu_Strand_trans_gminus',
    'IC_Glu_Strand_trans_gplus',
    'IC_Glu_Strand_trans_trans',
    'IC_Glu_Turn',
    'IC_Glu_Turn_gminus_gminus',
    'IC_Glu_Turn_gminus_gplus',
    'IC_Glu_Turn_gminus_trans',
    'IC_Glu_Turn_gplus_gminus',
    'IC_Glu_Turn_gplus_gplus',
    'IC_Glu_Turn_gplus_trans',
    'IC_Glu_Turn_trans_gminus',
    'IC_Glu_Turn_trans_gplus',
    'IC_Glu_Turn_trans_trans',
]
