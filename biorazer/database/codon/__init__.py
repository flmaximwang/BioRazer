"""Codon data tables and codon-usage frequency lookups.

``tables`` holds static codon-frequency data (e.g. the E. coli
``E_COLI_T_SNAPGENE_712`` / ``E_COLI_RT_SNAPGENE_712`` maps); ``usage`` wraps
the python-codon-tables package for per-species usage frequencies. This
``__init__`` re-exports every public name so
``from biorazer.database.codon import ...`` works for both.
"""

from .tables import E_COLI_RT_SNAPGENE_712, E_COLI_T_SNAPGENE_712
from .usage import (
    AVAILABLE_CODON_TABLES,
    get_codon_usage_table,
    get_codon_usage_table_by_aa,
)

__all__ = [
    "E_COLI_T_SNAPGENE_712",
    "E_COLI_RT_SNAPGENE_712",
    "AVAILABLE_CODON_TABLES",
    "get_codon_usage_table",
    "get_codon_usage_table_by_aa",
]
