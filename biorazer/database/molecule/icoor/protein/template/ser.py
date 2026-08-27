# -*- coding: utf-8 -*-
"""Ideal {secondary-structure x rotamer} {names[aa]} {var} ``InternalCoord`` templates.

Each template is the ideal heavy-atom topology for {aa} anchored at
{{N, CA, C}}, with the side chain grown to ``rotamer`` and the
``phi``/``psi``/``omega`` of ``ss`` carried on the instance.  See
:mod:`biorazer.database.molecule.icoor.protein.template` and
:mod:`._builder` for the build rule and the per-conformer caveat.
"""

from biorazer.database.molecule.icoor.protein.template._builder import build_template

RESN = "SER"

# === GENERATED IC_* named templates - do not edit below ===
# Regenerate with: /opt/envs/BioRazer/bin/python scripts/generate_internal_coord_template_named.py

IC_Ser_HelixAlpha = build_template('SER', 'alpha-helix', 'canonical')
IC_Ser_HelixAlpha_gminus = build_template('SER', 'alpha-helix', 'g-')
IC_Ser_HelixAlpha_trans = build_template('SER', 'alpha-helix', 't')
IC_Ser_HelixAlpha_gplus = build_template('SER', 'alpha-helix', 'g+')
IC_Ser_Helix310 = build_template('SER', '3-10-helix', 'canonical')
IC_Ser_Helix310_gminus = build_template('SER', '3-10-helix', 'g-')
IC_Ser_Helix310_trans = build_template('SER', '3-10-helix', 't')
IC_Ser_Helix310_gplus = build_template('SER', '3-10-helix', 'g+')
IC_Ser_HelixPi = build_template('SER', 'pi-helix', 'canonical')
IC_Ser_HelixPi_gminus = build_template('SER', 'pi-helix', 'g-')
IC_Ser_HelixPi_trans = build_template('SER', 'pi-helix', 't')
IC_Ser_HelixPi_gplus = build_template('SER', 'pi-helix', 'g+')
IC_Ser_HelixPPII = build_template('SER', 'polyproline-II', 'canonical')
IC_Ser_HelixPPII_gminus = build_template('SER', 'polyproline-II', 'g-')
IC_Ser_HelixPPII_trans = build_template('SER', 'polyproline-II', 't')
IC_Ser_HelixPPII_gplus = build_template('SER', 'polyproline-II', 'g+')
IC_Ser_Strand = build_template('SER', 'beta-strand', 'canonical')
IC_Ser_Strand_gminus = build_template('SER', 'beta-strand', 'g-')
IC_Ser_Strand_trans = build_template('SER', 'beta-strand', 't')
IC_Ser_Strand_gplus = build_template('SER', 'beta-strand', 'g+')
IC_Ser_StrandParallel = build_template('SER', 'parallel-beta-strand', 'canonical')
IC_Ser_StrandParallel_gminus = build_template('SER', 'parallel-beta-strand', 'g-')
IC_Ser_StrandParallel_trans = build_template('SER', 'parallel-beta-strand', 't')
IC_Ser_StrandParallel_gplus = build_template('SER', 'parallel-beta-strand', 'g+')
IC_Ser_StrandAntiParallel = build_template('SER', 'antiparallel-beta-strand', 'canonical')
IC_Ser_StrandAntiParallel_gminus = build_template('SER', 'antiparallel-beta-strand', 'g-')
IC_Ser_StrandAntiParallel_trans = build_template('SER', 'antiparallel-beta-strand', 't')
IC_Ser_StrandAntiParallel_gplus = build_template('SER', 'antiparallel-beta-strand', 'g+')
IC_Ser_Bridge = build_template('SER', 'beta-bridge', 'canonical')
IC_Ser_Bridge_gminus = build_template('SER', 'beta-bridge', 'g-')
IC_Ser_Bridge_trans = build_template('SER', 'beta-bridge', 't')
IC_Ser_Bridge_gplus = build_template('SER', 'beta-bridge', 'g+')
IC_Ser_Turn = build_template('SER', 'turn', 'canonical')
IC_Ser_Turn_gminus = build_template('SER', 'turn', 'g-')
IC_Ser_Turn_trans = build_template('SER', 'turn', 't')
IC_Ser_Turn_gplus = build_template('SER', 'turn', 'g+')
IC_Ser_Bend = build_template('SER', 'bend', 'canonical')
IC_Ser_Bend_gminus = build_template('SER', 'bend', 'g-')
IC_Ser_Bend_trans = build_template('SER', 'bend', 't')
IC_Ser_Bend_gplus = build_template('SER', 'bend', 'g+')
IC_Ser_Coil = build_template('SER', 'coil', 'canonical')
IC_Ser_Coil_gminus = build_template('SER', 'coil', 'g-')
IC_Ser_Coil_trans = build_template('SER', 'coil', 't')
IC_Ser_Coil_gplus = build_template('SER', 'coil', 'g+')
IC_Ser_CisPeptide = build_template('SER', 'cis-peptide-bond', 'canonical')
IC_Ser_CisPeptide_gminus = build_template('SER', 'cis-peptide-bond', 'g-')
IC_Ser_CisPeptide_trans = build_template('SER', 'cis-peptide-bond', 't')
IC_Ser_CisPeptide_gplus = build_template('SER', 'cis-peptide-bond', 'g+')

__all__ = [
    'IC_Ser_Bend',
    'IC_Ser_Bend_gminus',
    'IC_Ser_Bend_gplus',
    'IC_Ser_Bend_trans',
    'IC_Ser_Bridge',
    'IC_Ser_Bridge_gminus',
    'IC_Ser_Bridge_gplus',
    'IC_Ser_Bridge_trans',
    'IC_Ser_CisPeptide',
    'IC_Ser_CisPeptide_gminus',
    'IC_Ser_CisPeptide_gplus',
    'IC_Ser_CisPeptide_trans',
    'IC_Ser_Coil',
    'IC_Ser_Coil_gminus',
    'IC_Ser_Coil_gplus',
    'IC_Ser_Coil_trans',
    'IC_Ser_Helix310',
    'IC_Ser_Helix310_gminus',
    'IC_Ser_Helix310_gplus',
    'IC_Ser_Helix310_trans',
    'IC_Ser_HelixAlpha',
    'IC_Ser_HelixAlpha_gminus',
    'IC_Ser_HelixAlpha_gplus',
    'IC_Ser_HelixAlpha_trans',
    'IC_Ser_HelixPPII',
    'IC_Ser_HelixPPII_gminus',
    'IC_Ser_HelixPPII_gplus',
    'IC_Ser_HelixPPII_trans',
    'IC_Ser_HelixPi',
    'IC_Ser_HelixPi_gminus',
    'IC_Ser_HelixPi_gplus',
    'IC_Ser_HelixPi_trans',
    'IC_Ser_Strand',
    'IC_Ser_StrandAntiParallel',
    'IC_Ser_StrandAntiParallel_gminus',
    'IC_Ser_StrandAntiParallel_gplus',
    'IC_Ser_StrandAntiParallel_trans',
    'IC_Ser_StrandParallel',
    'IC_Ser_StrandParallel_gminus',
    'IC_Ser_StrandParallel_gplus',
    'IC_Ser_StrandParallel_trans',
    'IC_Ser_Strand_gminus',
    'IC_Ser_Strand_gplus',
    'IC_Ser_Strand_trans',
    'IC_Ser_Turn',
    'IC_Ser_Turn_gminus',
    'IC_Ser_Turn_gplus',
    'IC_Ser_Turn_trans',
]
