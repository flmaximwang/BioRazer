# -*- coding: utf-8 -*-
"""Name-token rules for the named ``IC_<Res>_<SS>[_<rotamer>]`` template
constants (shared by the package ``__init__`` and the generator script).

Secondary-structure class keys contain hyphens (``alpha-helix``,
``cis-peptide-bond`` ...) and rotamer keys contain ``-`` and ``/``
(``g-``, ``g-/t`` ...), none of which is a legal Python identifier -- so each
gets an identifier-safe token here.  Canonical templates are named
``IC_<Res>_<SS>``; rotamer variants append ``_<rotamer-token>``.
"""

#: Secondary-structure class -> identifier-safe token (per user naming,
#: ``alpha-helix`` is the fully-spelled ``HelixAlpha``).
SS_TOKEN = {
    "alpha-helix": "HelixAlpha",
    "3-10-helix": "Helix310",
    "pi-helix": "HelixPi",
    "polyproline-II": "HelixPPII",
    "beta-strand": "Strand",
    "parallel-beta-strand": "StrandParallel",
    "antiparallel-beta-strand": "StrandAntiParallel",
    "beta-bridge": "Bridge",
    "turn": "Turn",
    "bend": "Bend",
    "coil": "Coil",
    "cis-peptide-bond": "CisPeptide",
}

#: Rotamer key -> identifier-safe token (Dunbrack bin centers).
ROTAMER_TOKEN = {"g-": "gminus", "t": "trans", "g+": "gplus"}


def rotamer_token(rotamer):
    """Rotamer key (``'g-/t'``) -> identifier token (``'gminus_trans'``)."""
    return "_".join(ROTAMER_TOKEN[p] for p in rotamer.split("/"))


def template_token(resn, ss, rotamer="canonical"):
    """Identifier-safe ``IC_<Res>_<SS>[_<rotamer>]`` token for one template.

    ``resn`` is the 3-letter residue (upper-case), ``ss`` an
    :data:`SS_TOKEN` key, ``rotamer`` a rotamer key (``"canonical"`` for the
    no-suffix representative).
    """
    tok = f"IC_{resn.title()}_{SS_TOKEN[ss]}"
    if rotamer != "canonical":
        tok += "_" + rotamer_token(rotamer)
    return tok