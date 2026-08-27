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

IC_Asn_HelixAlpha = build_template('ASN', 'alpha-helix', 'canonical')
IC_Asn_HelixAlpha_gminus_gminus = build_template('ASN', 'alpha-helix', 'g-/g-')
IC_Asn_HelixAlpha_gminus_trans = build_template('ASN', 'alpha-helix', 'g-/t')
IC_Asn_HelixAlpha_gminus_gplus = build_template('ASN', 'alpha-helix', 'g-/g+')
IC_Asn_HelixAlpha_trans_gminus = build_template('ASN', 'alpha-helix', 't/g-')
IC_Asn_HelixAlpha_trans_trans = build_template('ASN', 'alpha-helix', 't/t')
IC_Asn_HelixAlpha_trans_gplus = build_template('ASN', 'alpha-helix', 't/g+')
IC_Asn_HelixAlpha_gplus_gminus = build_template('ASN', 'alpha-helix', 'g+/g-')
IC_Asn_HelixAlpha_gplus_trans = build_template('ASN', 'alpha-helix', 'g+/t')
IC_Asn_HelixAlpha_gplus_gplus = build_template('ASN', 'alpha-helix', 'g+/g+')
IC_Asn_Helix310 = build_template('ASN', '3-10-helix', 'canonical')
IC_Asn_Helix310_gminus_gminus = build_template('ASN', '3-10-helix', 'g-/g-')
IC_Asn_Helix310_gminus_trans = build_template('ASN', '3-10-helix', 'g-/t')
IC_Asn_Helix310_gminus_gplus = build_template('ASN', '3-10-helix', 'g-/g+')
IC_Asn_Helix310_trans_gminus = build_template('ASN', '3-10-helix', 't/g-')
IC_Asn_Helix310_trans_trans = build_template('ASN', '3-10-helix', 't/t')
IC_Asn_Helix310_trans_gplus = build_template('ASN', '3-10-helix', 't/g+')
IC_Asn_Helix310_gplus_gminus = build_template('ASN', '3-10-helix', 'g+/g-')
IC_Asn_Helix310_gplus_trans = build_template('ASN', '3-10-helix', 'g+/t')
IC_Asn_Helix310_gplus_gplus = build_template('ASN', '3-10-helix', 'g+/g+')
IC_Asn_HelixPi = build_template('ASN', 'pi-helix', 'canonical')
IC_Asn_HelixPi_gminus_gminus = build_template('ASN', 'pi-helix', 'g-/g-')
IC_Asn_HelixPi_gminus_trans = build_template('ASN', 'pi-helix', 'g-/t')
IC_Asn_HelixPi_gminus_gplus = build_template('ASN', 'pi-helix', 'g-/g+')
IC_Asn_HelixPi_trans_gminus = build_template('ASN', 'pi-helix', 't/g-')
IC_Asn_HelixPi_trans_trans = build_template('ASN', 'pi-helix', 't/t')
IC_Asn_HelixPi_trans_gplus = build_template('ASN', 'pi-helix', 't/g+')
IC_Asn_HelixPi_gplus_gminus = build_template('ASN', 'pi-helix', 'g+/g-')
IC_Asn_HelixPi_gplus_trans = build_template('ASN', 'pi-helix', 'g+/t')
IC_Asn_HelixPi_gplus_gplus = build_template('ASN', 'pi-helix', 'g+/g+')
IC_Asn_HelixPPII = build_template('ASN', 'polyproline-II', 'canonical')
IC_Asn_HelixPPII_gminus_gminus = build_template('ASN', 'polyproline-II', 'g-/g-')
IC_Asn_HelixPPII_gminus_trans = build_template('ASN', 'polyproline-II', 'g-/t')
IC_Asn_HelixPPII_gminus_gplus = build_template('ASN', 'polyproline-II', 'g-/g+')
IC_Asn_HelixPPII_trans_gminus = build_template('ASN', 'polyproline-II', 't/g-')
IC_Asn_HelixPPII_trans_trans = build_template('ASN', 'polyproline-II', 't/t')
IC_Asn_HelixPPII_trans_gplus = build_template('ASN', 'polyproline-II', 't/g+')
IC_Asn_HelixPPII_gplus_gminus = build_template('ASN', 'polyproline-II', 'g+/g-')
IC_Asn_HelixPPII_gplus_trans = build_template('ASN', 'polyproline-II', 'g+/t')
IC_Asn_HelixPPII_gplus_gplus = build_template('ASN', 'polyproline-II', 'g+/g+')
IC_Asn_Strand = build_template('ASN', 'beta-strand', 'canonical')
IC_Asn_Strand_gminus_gminus = build_template('ASN', 'beta-strand', 'g-/g-')
IC_Asn_Strand_gminus_trans = build_template('ASN', 'beta-strand', 'g-/t')
IC_Asn_Strand_gminus_gplus = build_template('ASN', 'beta-strand', 'g-/g+')
IC_Asn_Strand_trans_gminus = build_template('ASN', 'beta-strand', 't/g-')
IC_Asn_Strand_trans_trans = build_template('ASN', 'beta-strand', 't/t')
IC_Asn_Strand_trans_gplus = build_template('ASN', 'beta-strand', 't/g+')
IC_Asn_Strand_gplus_gminus = build_template('ASN', 'beta-strand', 'g+/g-')
IC_Asn_Strand_gplus_trans = build_template('ASN', 'beta-strand', 'g+/t')
IC_Asn_Strand_gplus_gplus = build_template('ASN', 'beta-strand', 'g+/g+')
IC_Asn_StrandParallel = build_template('ASN', 'parallel-beta-strand', 'canonical')
IC_Asn_StrandParallel_gminus_gminus = build_template('ASN', 'parallel-beta-strand', 'g-/g-')
IC_Asn_StrandParallel_gminus_trans = build_template('ASN', 'parallel-beta-strand', 'g-/t')
IC_Asn_StrandParallel_gminus_gplus = build_template('ASN', 'parallel-beta-strand', 'g-/g+')
IC_Asn_StrandParallel_trans_gminus = build_template('ASN', 'parallel-beta-strand', 't/g-')
IC_Asn_StrandParallel_trans_trans = build_template('ASN', 'parallel-beta-strand', 't/t')
IC_Asn_StrandParallel_trans_gplus = build_template('ASN', 'parallel-beta-strand', 't/g+')
IC_Asn_StrandParallel_gplus_gminus = build_template('ASN', 'parallel-beta-strand', 'g+/g-')
IC_Asn_StrandParallel_gplus_trans = build_template('ASN', 'parallel-beta-strand', 'g+/t')
IC_Asn_StrandParallel_gplus_gplus = build_template('ASN', 'parallel-beta-strand', 'g+/g+')
IC_Asn_StrandAntiParallel = build_template('ASN', 'antiparallel-beta-strand', 'canonical')
IC_Asn_StrandAntiParallel_gminus_gminus = build_template('ASN', 'antiparallel-beta-strand', 'g-/g-')
IC_Asn_StrandAntiParallel_gminus_trans = build_template('ASN', 'antiparallel-beta-strand', 'g-/t')
IC_Asn_StrandAntiParallel_gminus_gplus = build_template('ASN', 'antiparallel-beta-strand', 'g-/g+')
IC_Asn_StrandAntiParallel_trans_gminus = build_template('ASN', 'antiparallel-beta-strand', 't/g-')
IC_Asn_StrandAntiParallel_trans_trans = build_template('ASN', 'antiparallel-beta-strand', 't/t')
IC_Asn_StrandAntiParallel_trans_gplus = build_template('ASN', 'antiparallel-beta-strand', 't/g+')
IC_Asn_StrandAntiParallel_gplus_gminus = build_template('ASN', 'antiparallel-beta-strand', 'g+/g-')
IC_Asn_StrandAntiParallel_gplus_trans = build_template('ASN', 'antiparallel-beta-strand', 'g+/t')
IC_Asn_StrandAntiParallel_gplus_gplus = build_template('ASN', 'antiparallel-beta-strand', 'g+/g+')
IC_Asn_Bridge = build_template('ASN', 'beta-bridge', 'canonical')
IC_Asn_Bridge_gminus_gminus = build_template('ASN', 'beta-bridge', 'g-/g-')
IC_Asn_Bridge_gminus_trans = build_template('ASN', 'beta-bridge', 'g-/t')
IC_Asn_Bridge_gminus_gplus = build_template('ASN', 'beta-bridge', 'g-/g+')
IC_Asn_Bridge_trans_gminus = build_template('ASN', 'beta-bridge', 't/g-')
IC_Asn_Bridge_trans_trans = build_template('ASN', 'beta-bridge', 't/t')
IC_Asn_Bridge_trans_gplus = build_template('ASN', 'beta-bridge', 't/g+')
IC_Asn_Bridge_gplus_gminus = build_template('ASN', 'beta-bridge', 'g+/g-')
IC_Asn_Bridge_gplus_trans = build_template('ASN', 'beta-bridge', 'g+/t')
IC_Asn_Bridge_gplus_gplus = build_template('ASN', 'beta-bridge', 'g+/g+')
IC_Asn_Turn = build_template('ASN', 'turn', 'canonical')
IC_Asn_Turn_gminus_gminus = build_template('ASN', 'turn', 'g-/g-')
IC_Asn_Turn_gminus_trans = build_template('ASN', 'turn', 'g-/t')
IC_Asn_Turn_gminus_gplus = build_template('ASN', 'turn', 'g-/g+')
IC_Asn_Turn_trans_gminus = build_template('ASN', 'turn', 't/g-')
IC_Asn_Turn_trans_trans = build_template('ASN', 'turn', 't/t')
IC_Asn_Turn_trans_gplus = build_template('ASN', 'turn', 't/g+')
IC_Asn_Turn_gplus_gminus = build_template('ASN', 'turn', 'g+/g-')
IC_Asn_Turn_gplus_trans = build_template('ASN', 'turn', 'g+/t')
IC_Asn_Turn_gplus_gplus = build_template('ASN', 'turn', 'g+/g+')
IC_Asn_Bend = build_template('ASN', 'bend', 'canonical')
IC_Asn_Bend_gminus_gminus = build_template('ASN', 'bend', 'g-/g-')
IC_Asn_Bend_gminus_trans = build_template('ASN', 'bend', 'g-/t')
IC_Asn_Bend_gminus_gplus = build_template('ASN', 'bend', 'g-/g+')
IC_Asn_Bend_trans_gminus = build_template('ASN', 'bend', 't/g-')
IC_Asn_Bend_trans_trans = build_template('ASN', 'bend', 't/t')
IC_Asn_Bend_trans_gplus = build_template('ASN', 'bend', 't/g+')
IC_Asn_Bend_gplus_gminus = build_template('ASN', 'bend', 'g+/g-')
IC_Asn_Bend_gplus_trans = build_template('ASN', 'bend', 'g+/t')
IC_Asn_Bend_gplus_gplus = build_template('ASN', 'bend', 'g+/g+')
IC_Asn_Coil = build_template('ASN', 'coil', 'canonical')
IC_Asn_Coil_gminus_gminus = build_template('ASN', 'coil', 'g-/g-')
IC_Asn_Coil_gminus_trans = build_template('ASN', 'coil', 'g-/t')
IC_Asn_Coil_gminus_gplus = build_template('ASN', 'coil', 'g-/g+')
IC_Asn_Coil_trans_gminus = build_template('ASN', 'coil', 't/g-')
IC_Asn_Coil_trans_trans = build_template('ASN', 'coil', 't/t')
IC_Asn_Coil_trans_gplus = build_template('ASN', 'coil', 't/g+')
IC_Asn_Coil_gplus_gminus = build_template('ASN', 'coil', 'g+/g-')
IC_Asn_Coil_gplus_trans = build_template('ASN', 'coil', 'g+/t')
IC_Asn_Coil_gplus_gplus = build_template('ASN', 'coil', 'g+/g+')
IC_Asn_CisPeptide = build_template('ASN', 'cis-peptide-bond', 'canonical')
IC_Asn_CisPeptide_gminus_gminus = build_template('ASN', 'cis-peptide-bond', 'g-/g-')
IC_Asn_CisPeptide_gminus_trans = build_template('ASN', 'cis-peptide-bond', 'g-/t')
IC_Asn_CisPeptide_gminus_gplus = build_template('ASN', 'cis-peptide-bond', 'g-/g+')
IC_Asn_CisPeptide_trans_gminus = build_template('ASN', 'cis-peptide-bond', 't/g-')
IC_Asn_CisPeptide_trans_trans = build_template('ASN', 'cis-peptide-bond', 't/t')
IC_Asn_CisPeptide_trans_gplus = build_template('ASN', 'cis-peptide-bond', 't/g+')
IC_Asn_CisPeptide_gplus_gminus = build_template('ASN', 'cis-peptide-bond', 'g+/g-')
IC_Asn_CisPeptide_gplus_trans = build_template('ASN', 'cis-peptide-bond', 'g+/t')
IC_Asn_CisPeptide_gplus_gplus = build_template('ASN', 'cis-peptide-bond', 'g+/g+')

__all__ = [
    'IC_Asn_Bend',
    'IC_Asn_Bend_gminus_gminus',
    'IC_Asn_Bend_gminus_gplus',
    'IC_Asn_Bend_gminus_trans',
    'IC_Asn_Bend_gplus_gminus',
    'IC_Asn_Bend_gplus_gplus',
    'IC_Asn_Bend_gplus_trans',
    'IC_Asn_Bend_trans_gminus',
    'IC_Asn_Bend_trans_gplus',
    'IC_Asn_Bend_trans_trans',
    'IC_Asn_Bridge',
    'IC_Asn_Bridge_gminus_gminus',
    'IC_Asn_Bridge_gminus_gplus',
    'IC_Asn_Bridge_gminus_trans',
    'IC_Asn_Bridge_gplus_gminus',
    'IC_Asn_Bridge_gplus_gplus',
    'IC_Asn_Bridge_gplus_trans',
    'IC_Asn_Bridge_trans_gminus',
    'IC_Asn_Bridge_trans_gplus',
    'IC_Asn_Bridge_trans_trans',
    'IC_Asn_CisPeptide',
    'IC_Asn_CisPeptide_gminus_gminus',
    'IC_Asn_CisPeptide_gminus_gplus',
    'IC_Asn_CisPeptide_gminus_trans',
    'IC_Asn_CisPeptide_gplus_gminus',
    'IC_Asn_CisPeptide_gplus_gplus',
    'IC_Asn_CisPeptide_gplus_trans',
    'IC_Asn_CisPeptide_trans_gminus',
    'IC_Asn_CisPeptide_trans_gplus',
    'IC_Asn_CisPeptide_trans_trans',
    'IC_Asn_Coil',
    'IC_Asn_Coil_gminus_gminus',
    'IC_Asn_Coil_gminus_gplus',
    'IC_Asn_Coil_gminus_trans',
    'IC_Asn_Coil_gplus_gminus',
    'IC_Asn_Coil_gplus_gplus',
    'IC_Asn_Coil_gplus_trans',
    'IC_Asn_Coil_trans_gminus',
    'IC_Asn_Coil_trans_gplus',
    'IC_Asn_Coil_trans_trans',
    'IC_Asn_Helix310',
    'IC_Asn_Helix310_gminus_gminus',
    'IC_Asn_Helix310_gminus_gplus',
    'IC_Asn_Helix310_gminus_trans',
    'IC_Asn_Helix310_gplus_gminus',
    'IC_Asn_Helix310_gplus_gplus',
    'IC_Asn_Helix310_gplus_trans',
    'IC_Asn_Helix310_trans_gminus',
    'IC_Asn_Helix310_trans_gplus',
    'IC_Asn_Helix310_trans_trans',
    'IC_Asn_HelixAlpha',
    'IC_Asn_HelixAlpha_gminus_gminus',
    'IC_Asn_HelixAlpha_gminus_gplus',
    'IC_Asn_HelixAlpha_gminus_trans',
    'IC_Asn_HelixAlpha_gplus_gminus',
    'IC_Asn_HelixAlpha_gplus_gplus',
    'IC_Asn_HelixAlpha_gplus_trans',
    'IC_Asn_HelixAlpha_trans_gminus',
    'IC_Asn_HelixAlpha_trans_gplus',
    'IC_Asn_HelixAlpha_trans_trans',
    'IC_Asn_HelixPPII',
    'IC_Asn_HelixPPII_gminus_gminus',
    'IC_Asn_HelixPPII_gminus_gplus',
    'IC_Asn_HelixPPII_gminus_trans',
    'IC_Asn_HelixPPII_gplus_gminus',
    'IC_Asn_HelixPPII_gplus_gplus',
    'IC_Asn_HelixPPII_gplus_trans',
    'IC_Asn_HelixPPII_trans_gminus',
    'IC_Asn_HelixPPII_trans_gplus',
    'IC_Asn_HelixPPII_trans_trans',
    'IC_Asn_HelixPi',
    'IC_Asn_HelixPi_gminus_gminus',
    'IC_Asn_HelixPi_gminus_gplus',
    'IC_Asn_HelixPi_gminus_trans',
    'IC_Asn_HelixPi_gplus_gminus',
    'IC_Asn_HelixPi_gplus_gplus',
    'IC_Asn_HelixPi_gplus_trans',
    'IC_Asn_HelixPi_trans_gminus',
    'IC_Asn_HelixPi_trans_gplus',
    'IC_Asn_HelixPi_trans_trans',
    'IC_Asn_Strand',
    'IC_Asn_StrandAntiParallel',
    'IC_Asn_StrandAntiParallel_gminus_gminus',
    'IC_Asn_StrandAntiParallel_gminus_gplus',
    'IC_Asn_StrandAntiParallel_gminus_trans',
    'IC_Asn_StrandAntiParallel_gplus_gminus',
    'IC_Asn_StrandAntiParallel_gplus_gplus',
    'IC_Asn_StrandAntiParallel_gplus_trans',
    'IC_Asn_StrandAntiParallel_trans_gminus',
    'IC_Asn_StrandAntiParallel_trans_gplus',
    'IC_Asn_StrandAntiParallel_trans_trans',
    'IC_Asn_StrandParallel',
    'IC_Asn_StrandParallel_gminus_gminus',
    'IC_Asn_StrandParallel_gminus_gplus',
    'IC_Asn_StrandParallel_gminus_trans',
    'IC_Asn_StrandParallel_gplus_gminus',
    'IC_Asn_StrandParallel_gplus_gplus',
    'IC_Asn_StrandParallel_gplus_trans',
    'IC_Asn_StrandParallel_trans_gminus',
    'IC_Asn_StrandParallel_trans_gplus',
    'IC_Asn_StrandParallel_trans_trans',
    'IC_Asn_Strand_gminus_gminus',
    'IC_Asn_Strand_gminus_gplus',
    'IC_Asn_Strand_gminus_trans',
    'IC_Asn_Strand_gplus_gminus',
    'IC_Asn_Strand_gplus_gplus',
    'IC_Asn_Strand_gplus_trans',
    'IC_Asn_Strand_trans_gminus',
    'IC_Asn_Strand_trans_gplus',
    'IC_Asn_Strand_trans_trans',
    'IC_Asn_Turn',
    'IC_Asn_Turn_gminus_gminus',
    'IC_Asn_Turn_gminus_gplus',
    'IC_Asn_Turn_gminus_trans',
    'IC_Asn_Turn_gplus_gminus',
    'IC_Asn_Turn_gplus_gplus',
    'IC_Asn_Turn_gplus_trans',
    'IC_Asn_Turn_trans_gminus',
    'IC_Asn_Turn_trans_gplus',
    'IC_Asn_Turn_trans_trans',
]
