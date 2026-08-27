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

IC_Cys_HelixAlpha = build_template('CYS', 'alpha-helix', 'canonical')
IC_Cys_HelixAlpha_gminus = build_template('CYS', 'alpha-helix', 'g-')
IC_Cys_HelixAlpha_trans = build_template('CYS', 'alpha-helix', 't')
IC_Cys_HelixAlpha_gplus = build_template('CYS', 'alpha-helix', 'g+')
IC_Cys_Helix310 = build_template('CYS', '3-10-helix', 'canonical')
IC_Cys_Helix310_gminus = build_template('CYS', '3-10-helix', 'g-')
IC_Cys_Helix310_trans = build_template('CYS', '3-10-helix', 't')
IC_Cys_Helix310_gplus = build_template('CYS', '3-10-helix', 'g+')
IC_Cys_HelixPi = build_template('CYS', 'pi-helix', 'canonical')
IC_Cys_HelixPi_gminus = build_template('CYS', 'pi-helix', 'g-')
IC_Cys_HelixPi_trans = build_template('CYS', 'pi-helix', 't')
IC_Cys_HelixPi_gplus = build_template('CYS', 'pi-helix', 'g+')
IC_Cys_HelixPPII = build_template('CYS', 'polyproline-II', 'canonical')
IC_Cys_HelixPPII_gminus = build_template('CYS', 'polyproline-II', 'g-')
IC_Cys_HelixPPII_trans = build_template('CYS', 'polyproline-II', 't')
IC_Cys_HelixPPII_gplus = build_template('CYS', 'polyproline-II', 'g+')
IC_Cys_Strand = build_template('CYS', 'beta-strand', 'canonical')
IC_Cys_Strand_gminus = build_template('CYS', 'beta-strand', 'g-')
IC_Cys_Strand_trans = build_template('CYS', 'beta-strand', 't')
IC_Cys_Strand_gplus = build_template('CYS', 'beta-strand', 'g+')
IC_Cys_StrandParallel = build_template('CYS', 'parallel-beta-strand', 'canonical')
IC_Cys_StrandParallel_gminus = build_template('CYS', 'parallel-beta-strand', 'g-')
IC_Cys_StrandParallel_trans = build_template('CYS', 'parallel-beta-strand', 't')
IC_Cys_StrandParallel_gplus = build_template('CYS', 'parallel-beta-strand', 'g+')
IC_Cys_StrandAntiParallel = build_template('CYS', 'antiparallel-beta-strand', 'canonical')
IC_Cys_StrandAntiParallel_gminus = build_template('CYS', 'antiparallel-beta-strand', 'g-')
IC_Cys_StrandAntiParallel_trans = build_template('CYS', 'antiparallel-beta-strand', 't')
IC_Cys_StrandAntiParallel_gplus = build_template('CYS', 'antiparallel-beta-strand', 'g+')
IC_Cys_Bridge = build_template('CYS', 'beta-bridge', 'canonical')
IC_Cys_Bridge_gminus = build_template('CYS', 'beta-bridge', 'g-')
IC_Cys_Bridge_trans = build_template('CYS', 'beta-bridge', 't')
IC_Cys_Bridge_gplus = build_template('CYS', 'beta-bridge', 'g+')
IC_Cys_Turn = build_template('CYS', 'turn', 'canonical')
IC_Cys_Turn_gminus = build_template('CYS', 'turn', 'g-')
IC_Cys_Turn_trans = build_template('CYS', 'turn', 't')
IC_Cys_Turn_gplus = build_template('CYS', 'turn', 'g+')
IC_Cys_Bend = build_template('CYS', 'bend', 'canonical')
IC_Cys_Bend_gminus = build_template('CYS', 'bend', 'g-')
IC_Cys_Bend_trans = build_template('CYS', 'bend', 't')
IC_Cys_Bend_gplus = build_template('CYS', 'bend', 'g+')
IC_Cys_Coil = build_template('CYS', 'coil', 'canonical')
IC_Cys_Coil_gminus = build_template('CYS', 'coil', 'g-')
IC_Cys_Coil_trans = build_template('CYS', 'coil', 't')
IC_Cys_Coil_gplus = build_template('CYS', 'coil', 'g+')
IC_Cys_CisPeptide = build_template('CYS', 'cis-peptide-bond', 'canonical')
IC_Cys_CisPeptide_gminus = build_template('CYS', 'cis-peptide-bond', 'g-')
IC_Cys_CisPeptide_trans = build_template('CYS', 'cis-peptide-bond', 't')
IC_Cys_CisPeptide_gplus = build_template('CYS', 'cis-peptide-bond', 'g+')

__all__ = [
    'IC_Cys_Bend',
    'IC_Cys_Bend_gminus',
    'IC_Cys_Bend_gplus',
    'IC_Cys_Bend_trans',
    'IC_Cys_Bridge',
    'IC_Cys_Bridge_gminus',
    'IC_Cys_Bridge_gplus',
    'IC_Cys_Bridge_trans',
    'IC_Cys_CisPeptide',
    'IC_Cys_CisPeptide_gminus',
    'IC_Cys_CisPeptide_gplus',
    'IC_Cys_CisPeptide_trans',
    'IC_Cys_Coil',
    'IC_Cys_Coil_gminus',
    'IC_Cys_Coil_gplus',
    'IC_Cys_Coil_trans',
    'IC_Cys_Helix310',
    'IC_Cys_Helix310_gminus',
    'IC_Cys_Helix310_gplus',
    'IC_Cys_Helix310_trans',
    'IC_Cys_HelixAlpha',
    'IC_Cys_HelixAlpha_gminus',
    'IC_Cys_HelixAlpha_gplus',
    'IC_Cys_HelixAlpha_trans',
    'IC_Cys_HelixPPII',
    'IC_Cys_HelixPPII_gminus',
    'IC_Cys_HelixPPII_gplus',
    'IC_Cys_HelixPPII_trans',
    'IC_Cys_HelixPi',
    'IC_Cys_HelixPi_gminus',
    'IC_Cys_HelixPi_gplus',
    'IC_Cys_HelixPi_trans',
    'IC_Cys_Strand',
    'IC_Cys_StrandAntiParallel',
    'IC_Cys_StrandAntiParallel_gminus',
    'IC_Cys_StrandAntiParallel_gplus',
    'IC_Cys_StrandAntiParallel_trans',
    'IC_Cys_StrandParallel',
    'IC_Cys_StrandParallel_gminus',
    'IC_Cys_StrandParallel_gplus',
    'IC_Cys_StrandParallel_trans',
    'IC_Cys_Strand_gminus',
    'IC_Cys_Strand_gplus',
    'IC_Cys_Strand_trans',
    'IC_Cys_Turn',
    'IC_Cys_Turn_gminus',
    'IC_Cys_Turn_gplus',
    'IC_Cys_Turn_trans',
]
