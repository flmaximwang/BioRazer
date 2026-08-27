# -*- coding: utf-8 -*-
"""Ideal {secondary-structure x rotamer} {names[aa]} {var} ``InternalCoord`` templates.

Each template is the ideal heavy-atom topology for {aa} anchored at
{{N, CA, C}}, with the side chain grown to ``rotamer`` and the
``phi``/``psi``/``omega`` of ``ss`` carried on the instance.  See
:mod:`biorazer.database.molecule.icoor.protein.template` and
:mod:`._builder` for the build rule and the per-conformer caveat.
"""

from biorazer.database.molecule.icoor.protein.template._builder import build_template

RESN = "LYS"

# === GENERATED IC_* named templates - do not edit below ===
# Regenerate with: /opt/envs/BioRazer/bin/python scripts/generate_internal_coord_template_named.py

IC_Lys_HelixAlpha = build_template('LYS', 'alpha-helix', 'canonical')
IC_Lys_HelixAlpha_gminus_gminus = build_template('LYS', 'alpha-helix', 'g-/g-')
IC_Lys_HelixAlpha_gminus_trans = build_template('LYS', 'alpha-helix', 'g-/t')
IC_Lys_HelixAlpha_gminus_gplus = build_template('LYS', 'alpha-helix', 'g-/g+')
IC_Lys_HelixAlpha_trans_gminus = build_template('LYS', 'alpha-helix', 't/g-')
IC_Lys_HelixAlpha_trans_trans = build_template('LYS', 'alpha-helix', 't/t')
IC_Lys_HelixAlpha_trans_gplus = build_template('LYS', 'alpha-helix', 't/g+')
IC_Lys_HelixAlpha_gplus_gminus = build_template('LYS', 'alpha-helix', 'g+/g-')
IC_Lys_HelixAlpha_gplus_trans = build_template('LYS', 'alpha-helix', 'g+/t')
IC_Lys_HelixAlpha_gplus_gplus = build_template('LYS', 'alpha-helix', 'g+/g+')
IC_Lys_Helix310 = build_template('LYS', '3-10-helix', 'canonical')
IC_Lys_Helix310_gminus_gminus = build_template('LYS', '3-10-helix', 'g-/g-')
IC_Lys_Helix310_gminus_trans = build_template('LYS', '3-10-helix', 'g-/t')
IC_Lys_Helix310_gminus_gplus = build_template('LYS', '3-10-helix', 'g-/g+')
IC_Lys_Helix310_trans_gminus = build_template('LYS', '3-10-helix', 't/g-')
IC_Lys_Helix310_trans_trans = build_template('LYS', '3-10-helix', 't/t')
IC_Lys_Helix310_trans_gplus = build_template('LYS', '3-10-helix', 't/g+')
IC_Lys_Helix310_gplus_gminus = build_template('LYS', '3-10-helix', 'g+/g-')
IC_Lys_Helix310_gplus_trans = build_template('LYS', '3-10-helix', 'g+/t')
IC_Lys_Helix310_gplus_gplus = build_template('LYS', '3-10-helix', 'g+/g+')
IC_Lys_HelixPi = build_template('LYS', 'pi-helix', 'canonical')
IC_Lys_HelixPi_gminus_gminus = build_template('LYS', 'pi-helix', 'g-/g-')
IC_Lys_HelixPi_gminus_trans = build_template('LYS', 'pi-helix', 'g-/t')
IC_Lys_HelixPi_gminus_gplus = build_template('LYS', 'pi-helix', 'g-/g+')
IC_Lys_HelixPi_trans_gminus = build_template('LYS', 'pi-helix', 't/g-')
IC_Lys_HelixPi_trans_trans = build_template('LYS', 'pi-helix', 't/t')
IC_Lys_HelixPi_trans_gplus = build_template('LYS', 'pi-helix', 't/g+')
IC_Lys_HelixPi_gplus_gminus = build_template('LYS', 'pi-helix', 'g+/g-')
IC_Lys_HelixPi_gplus_trans = build_template('LYS', 'pi-helix', 'g+/t')
IC_Lys_HelixPi_gplus_gplus = build_template('LYS', 'pi-helix', 'g+/g+')
IC_Lys_HelixPPII = build_template('LYS', 'polyproline-II', 'canonical')
IC_Lys_HelixPPII_gminus_gminus = build_template('LYS', 'polyproline-II', 'g-/g-')
IC_Lys_HelixPPII_gminus_trans = build_template('LYS', 'polyproline-II', 'g-/t')
IC_Lys_HelixPPII_gminus_gplus = build_template('LYS', 'polyproline-II', 'g-/g+')
IC_Lys_HelixPPII_trans_gminus = build_template('LYS', 'polyproline-II', 't/g-')
IC_Lys_HelixPPII_trans_trans = build_template('LYS', 'polyproline-II', 't/t')
IC_Lys_HelixPPII_trans_gplus = build_template('LYS', 'polyproline-II', 't/g+')
IC_Lys_HelixPPII_gplus_gminus = build_template('LYS', 'polyproline-II', 'g+/g-')
IC_Lys_HelixPPII_gplus_trans = build_template('LYS', 'polyproline-II', 'g+/t')
IC_Lys_HelixPPII_gplus_gplus = build_template('LYS', 'polyproline-II', 'g+/g+')
IC_Lys_Strand = build_template('LYS', 'beta-strand', 'canonical')
IC_Lys_Strand_gminus_gminus = build_template('LYS', 'beta-strand', 'g-/g-')
IC_Lys_Strand_gminus_trans = build_template('LYS', 'beta-strand', 'g-/t')
IC_Lys_Strand_gminus_gplus = build_template('LYS', 'beta-strand', 'g-/g+')
IC_Lys_Strand_trans_gminus = build_template('LYS', 'beta-strand', 't/g-')
IC_Lys_Strand_trans_trans = build_template('LYS', 'beta-strand', 't/t')
IC_Lys_Strand_trans_gplus = build_template('LYS', 'beta-strand', 't/g+')
IC_Lys_Strand_gplus_gminus = build_template('LYS', 'beta-strand', 'g+/g-')
IC_Lys_Strand_gplus_trans = build_template('LYS', 'beta-strand', 'g+/t')
IC_Lys_Strand_gplus_gplus = build_template('LYS', 'beta-strand', 'g+/g+')
IC_Lys_StrandParallel = build_template('LYS', 'parallel-beta-strand', 'canonical')
IC_Lys_StrandParallel_gminus_gminus = build_template('LYS', 'parallel-beta-strand', 'g-/g-')
IC_Lys_StrandParallel_gminus_trans = build_template('LYS', 'parallel-beta-strand', 'g-/t')
IC_Lys_StrandParallel_gminus_gplus = build_template('LYS', 'parallel-beta-strand', 'g-/g+')
IC_Lys_StrandParallel_trans_gminus = build_template('LYS', 'parallel-beta-strand', 't/g-')
IC_Lys_StrandParallel_trans_trans = build_template('LYS', 'parallel-beta-strand', 't/t')
IC_Lys_StrandParallel_trans_gplus = build_template('LYS', 'parallel-beta-strand', 't/g+')
IC_Lys_StrandParallel_gplus_gminus = build_template('LYS', 'parallel-beta-strand', 'g+/g-')
IC_Lys_StrandParallel_gplus_trans = build_template('LYS', 'parallel-beta-strand', 'g+/t')
IC_Lys_StrandParallel_gplus_gplus = build_template('LYS', 'parallel-beta-strand', 'g+/g+')
IC_Lys_StrandAntiParallel = build_template('LYS', 'antiparallel-beta-strand', 'canonical')
IC_Lys_StrandAntiParallel_gminus_gminus = build_template('LYS', 'antiparallel-beta-strand', 'g-/g-')
IC_Lys_StrandAntiParallel_gminus_trans = build_template('LYS', 'antiparallel-beta-strand', 'g-/t')
IC_Lys_StrandAntiParallel_gminus_gplus = build_template('LYS', 'antiparallel-beta-strand', 'g-/g+')
IC_Lys_StrandAntiParallel_trans_gminus = build_template('LYS', 'antiparallel-beta-strand', 't/g-')
IC_Lys_StrandAntiParallel_trans_trans = build_template('LYS', 'antiparallel-beta-strand', 't/t')
IC_Lys_StrandAntiParallel_trans_gplus = build_template('LYS', 'antiparallel-beta-strand', 't/g+')
IC_Lys_StrandAntiParallel_gplus_gminus = build_template('LYS', 'antiparallel-beta-strand', 'g+/g-')
IC_Lys_StrandAntiParallel_gplus_trans = build_template('LYS', 'antiparallel-beta-strand', 'g+/t')
IC_Lys_StrandAntiParallel_gplus_gplus = build_template('LYS', 'antiparallel-beta-strand', 'g+/g+')
IC_Lys_Bridge = build_template('LYS', 'beta-bridge', 'canonical')
IC_Lys_Bridge_gminus_gminus = build_template('LYS', 'beta-bridge', 'g-/g-')
IC_Lys_Bridge_gminus_trans = build_template('LYS', 'beta-bridge', 'g-/t')
IC_Lys_Bridge_gminus_gplus = build_template('LYS', 'beta-bridge', 'g-/g+')
IC_Lys_Bridge_trans_gminus = build_template('LYS', 'beta-bridge', 't/g-')
IC_Lys_Bridge_trans_trans = build_template('LYS', 'beta-bridge', 't/t')
IC_Lys_Bridge_trans_gplus = build_template('LYS', 'beta-bridge', 't/g+')
IC_Lys_Bridge_gplus_gminus = build_template('LYS', 'beta-bridge', 'g+/g-')
IC_Lys_Bridge_gplus_trans = build_template('LYS', 'beta-bridge', 'g+/t')
IC_Lys_Bridge_gplus_gplus = build_template('LYS', 'beta-bridge', 'g+/g+')
IC_Lys_Turn = build_template('LYS', 'turn', 'canonical')
IC_Lys_Turn_gminus_gminus = build_template('LYS', 'turn', 'g-/g-')
IC_Lys_Turn_gminus_trans = build_template('LYS', 'turn', 'g-/t')
IC_Lys_Turn_gminus_gplus = build_template('LYS', 'turn', 'g-/g+')
IC_Lys_Turn_trans_gminus = build_template('LYS', 'turn', 't/g-')
IC_Lys_Turn_trans_trans = build_template('LYS', 'turn', 't/t')
IC_Lys_Turn_trans_gplus = build_template('LYS', 'turn', 't/g+')
IC_Lys_Turn_gplus_gminus = build_template('LYS', 'turn', 'g+/g-')
IC_Lys_Turn_gplus_trans = build_template('LYS', 'turn', 'g+/t')
IC_Lys_Turn_gplus_gplus = build_template('LYS', 'turn', 'g+/g+')
IC_Lys_Bend = build_template('LYS', 'bend', 'canonical')
IC_Lys_Bend_gminus_gminus = build_template('LYS', 'bend', 'g-/g-')
IC_Lys_Bend_gminus_trans = build_template('LYS', 'bend', 'g-/t')
IC_Lys_Bend_gminus_gplus = build_template('LYS', 'bend', 'g-/g+')
IC_Lys_Bend_trans_gminus = build_template('LYS', 'bend', 't/g-')
IC_Lys_Bend_trans_trans = build_template('LYS', 'bend', 't/t')
IC_Lys_Bend_trans_gplus = build_template('LYS', 'bend', 't/g+')
IC_Lys_Bend_gplus_gminus = build_template('LYS', 'bend', 'g+/g-')
IC_Lys_Bend_gplus_trans = build_template('LYS', 'bend', 'g+/t')
IC_Lys_Bend_gplus_gplus = build_template('LYS', 'bend', 'g+/g+')
IC_Lys_Coil = build_template('LYS', 'coil', 'canonical')
IC_Lys_Coil_gminus_gminus = build_template('LYS', 'coil', 'g-/g-')
IC_Lys_Coil_gminus_trans = build_template('LYS', 'coil', 'g-/t')
IC_Lys_Coil_gminus_gplus = build_template('LYS', 'coil', 'g-/g+')
IC_Lys_Coil_trans_gminus = build_template('LYS', 'coil', 't/g-')
IC_Lys_Coil_trans_trans = build_template('LYS', 'coil', 't/t')
IC_Lys_Coil_trans_gplus = build_template('LYS', 'coil', 't/g+')
IC_Lys_Coil_gplus_gminus = build_template('LYS', 'coil', 'g+/g-')
IC_Lys_Coil_gplus_trans = build_template('LYS', 'coil', 'g+/t')
IC_Lys_Coil_gplus_gplus = build_template('LYS', 'coil', 'g+/g+')
IC_Lys_CisPeptide = build_template('LYS', 'cis-peptide-bond', 'canonical')
IC_Lys_CisPeptide_gminus_gminus = build_template('LYS', 'cis-peptide-bond', 'g-/g-')
IC_Lys_CisPeptide_gminus_trans = build_template('LYS', 'cis-peptide-bond', 'g-/t')
IC_Lys_CisPeptide_gminus_gplus = build_template('LYS', 'cis-peptide-bond', 'g-/g+')
IC_Lys_CisPeptide_trans_gminus = build_template('LYS', 'cis-peptide-bond', 't/g-')
IC_Lys_CisPeptide_trans_trans = build_template('LYS', 'cis-peptide-bond', 't/t')
IC_Lys_CisPeptide_trans_gplus = build_template('LYS', 'cis-peptide-bond', 't/g+')
IC_Lys_CisPeptide_gplus_gminus = build_template('LYS', 'cis-peptide-bond', 'g+/g-')
IC_Lys_CisPeptide_gplus_trans = build_template('LYS', 'cis-peptide-bond', 'g+/t')
IC_Lys_CisPeptide_gplus_gplus = build_template('LYS', 'cis-peptide-bond', 'g+/g+')

__all__ = [
    'IC_Lys_Bend',
    'IC_Lys_Bend_gminus_gminus',
    'IC_Lys_Bend_gminus_gplus',
    'IC_Lys_Bend_gminus_trans',
    'IC_Lys_Bend_gplus_gminus',
    'IC_Lys_Bend_gplus_gplus',
    'IC_Lys_Bend_gplus_trans',
    'IC_Lys_Bend_trans_gminus',
    'IC_Lys_Bend_trans_gplus',
    'IC_Lys_Bend_trans_trans',
    'IC_Lys_Bridge',
    'IC_Lys_Bridge_gminus_gminus',
    'IC_Lys_Bridge_gminus_gplus',
    'IC_Lys_Bridge_gminus_trans',
    'IC_Lys_Bridge_gplus_gminus',
    'IC_Lys_Bridge_gplus_gplus',
    'IC_Lys_Bridge_gplus_trans',
    'IC_Lys_Bridge_trans_gminus',
    'IC_Lys_Bridge_trans_gplus',
    'IC_Lys_Bridge_trans_trans',
    'IC_Lys_CisPeptide',
    'IC_Lys_CisPeptide_gminus_gminus',
    'IC_Lys_CisPeptide_gminus_gplus',
    'IC_Lys_CisPeptide_gminus_trans',
    'IC_Lys_CisPeptide_gplus_gminus',
    'IC_Lys_CisPeptide_gplus_gplus',
    'IC_Lys_CisPeptide_gplus_trans',
    'IC_Lys_CisPeptide_trans_gminus',
    'IC_Lys_CisPeptide_trans_gplus',
    'IC_Lys_CisPeptide_trans_trans',
    'IC_Lys_Coil',
    'IC_Lys_Coil_gminus_gminus',
    'IC_Lys_Coil_gminus_gplus',
    'IC_Lys_Coil_gminus_trans',
    'IC_Lys_Coil_gplus_gminus',
    'IC_Lys_Coil_gplus_gplus',
    'IC_Lys_Coil_gplus_trans',
    'IC_Lys_Coil_trans_gminus',
    'IC_Lys_Coil_trans_gplus',
    'IC_Lys_Coil_trans_trans',
    'IC_Lys_Helix310',
    'IC_Lys_Helix310_gminus_gminus',
    'IC_Lys_Helix310_gminus_gplus',
    'IC_Lys_Helix310_gminus_trans',
    'IC_Lys_Helix310_gplus_gminus',
    'IC_Lys_Helix310_gplus_gplus',
    'IC_Lys_Helix310_gplus_trans',
    'IC_Lys_Helix310_trans_gminus',
    'IC_Lys_Helix310_trans_gplus',
    'IC_Lys_Helix310_trans_trans',
    'IC_Lys_HelixAlpha',
    'IC_Lys_HelixAlpha_gminus_gminus',
    'IC_Lys_HelixAlpha_gminus_gplus',
    'IC_Lys_HelixAlpha_gminus_trans',
    'IC_Lys_HelixAlpha_gplus_gminus',
    'IC_Lys_HelixAlpha_gplus_gplus',
    'IC_Lys_HelixAlpha_gplus_trans',
    'IC_Lys_HelixAlpha_trans_gminus',
    'IC_Lys_HelixAlpha_trans_gplus',
    'IC_Lys_HelixAlpha_trans_trans',
    'IC_Lys_HelixPPII',
    'IC_Lys_HelixPPII_gminus_gminus',
    'IC_Lys_HelixPPII_gminus_gplus',
    'IC_Lys_HelixPPII_gminus_trans',
    'IC_Lys_HelixPPII_gplus_gminus',
    'IC_Lys_HelixPPII_gplus_gplus',
    'IC_Lys_HelixPPII_gplus_trans',
    'IC_Lys_HelixPPII_trans_gminus',
    'IC_Lys_HelixPPII_trans_gplus',
    'IC_Lys_HelixPPII_trans_trans',
    'IC_Lys_HelixPi',
    'IC_Lys_HelixPi_gminus_gminus',
    'IC_Lys_HelixPi_gminus_gplus',
    'IC_Lys_HelixPi_gminus_trans',
    'IC_Lys_HelixPi_gplus_gminus',
    'IC_Lys_HelixPi_gplus_gplus',
    'IC_Lys_HelixPi_gplus_trans',
    'IC_Lys_HelixPi_trans_gminus',
    'IC_Lys_HelixPi_trans_gplus',
    'IC_Lys_HelixPi_trans_trans',
    'IC_Lys_Strand',
    'IC_Lys_StrandAntiParallel',
    'IC_Lys_StrandAntiParallel_gminus_gminus',
    'IC_Lys_StrandAntiParallel_gminus_gplus',
    'IC_Lys_StrandAntiParallel_gminus_trans',
    'IC_Lys_StrandAntiParallel_gplus_gminus',
    'IC_Lys_StrandAntiParallel_gplus_gplus',
    'IC_Lys_StrandAntiParallel_gplus_trans',
    'IC_Lys_StrandAntiParallel_trans_gminus',
    'IC_Lys_StrandAntiParallel_trans_gplus',
    'IC_Lys_StrandAntiParallel_trans_trans',
    'IC_Lys_StrandParallel',
    'IC_Lys_StrandParallel_gminus_gminus',
    'IC_Lys_StrandParallel_gminus_gplus',
    'IC_Lys_StrandParallel_gminus_trans',
    'IC_Lys_StrandParallel_gplus_gminus',
    'IC_Lys_StrandParallel_gplus_gplus',
    'IC_Lys_StrandParallel_gplus_trans',
    'IC_Lys_StrandParallel_trans_gminus',
    'IC_Lys_StrandParallel_trans_gplus',
    'IC_Lys_StrandParallel_trans_trans',
    'IC_Lys_Strand_gminus_gminus',
    'IC_Lys_Strand_gminus_gplus',
    'IC_Lys_Strand_gminus_trans',
    'IC_Lys_Strand_gplus_gminus',
    'IC_Lys_Strand_gplus_gplus',
    'IC_Lys_Strand_gplus_trans',
    'IC_Lys_Strand_trans_gminus',
    'IC_Lys_Strand_trans_gplus',
    'IC_Lys_Strand_trans_trans',
    'IC_Lys_Turn',
    'IC_Lys_Turn_gminus_gminus',
    'IC_Lys_Turn_gminus_gplus',
    'IC_Lys_Turn_gminus_trans',
    'IC_Lys_Turn_gplus_gminus',
    'IC_Lys_Turn_gplus_gplus',
    'IC_Lys_Turn_gplus_trans',
    'IC_Lys_Turn_trans_gminus',
    'IC_Lys_Turn_trans_gplus',
    'IC_Lys_Turn_trans_trans',
]
