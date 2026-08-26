"""Amino-acid classification (physicochemical type -> 1-letter codes).

Split out of ``biorazer.database.amino_acid`` (originally ``aa_types.py``).
"""

TYPES2AA = {
    "hydrophobic": "AFILMPVWY",
    "polar": "CDEGHKNRST",
    "charged": "DEHKR",
    "positive": "HKR",
    "negative": "DE",
    "aromatic": "FHWY",
    "aliphatic": "AILMPV",
    "small": "ACDGNPSTV",
    "large": "EFHIKLRWY",
    "all": "ACDEFGHIKLMNPQRSTVWY",
}
