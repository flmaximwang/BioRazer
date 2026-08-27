# -*- coding: utf-8 -*-
"""Ideal torsion (dihedral) angles (°) and torsion definitions.

* :mod:`.generic` -- molecule-agnostic torsion tables (empty placeholder;
  for future non-protein ligands).
* :mod:`.protein` -- protein torsion data, split into
  :mod:`.protein.by_ss` (backbone phi/psi/omega by secondary structure,
  DSSP-classified, beta-turn dihedrals, the cis/trans omega constants,
  official main-chain torsion definitions) and :mod:`.protein.by_residue`
  (side-chain chi definitions / canonical IC-frame dihedrals / rotamer bin
  framework, Dunbrack).

Every numeric entry is ``{mean, std, lb, up, source}`` (``np.nan`` when
unknown).
"""
