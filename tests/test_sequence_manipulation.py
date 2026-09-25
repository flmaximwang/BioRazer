"""``biorazer.sequence.manipulation.util`` and the Ensembl module that uses it."""
import importlib

from biorazer.sequence.manipulation.util import get_mRNA_from_transcript


def test_exon_coordinates_are_end_inclusive():
    transcript = "AAACCCGGGTTT"

    # (0, 2) -> "AAA", (6, 8) -> "GGG"
    assert get_mRNA_from_transcript(transcript, [(0, 2), (6, 8)]) == "AAAGGG"
    # a single-base exon keeps its last base: [1:2], not [1:1]
    assert get_mRNA_from_transcript(transcript, [(1, 1)]) == "A"


def test_ensembl_access_module_imports():
    """``access.database.ensembl.utils`` imported a module that never existed.

    It reached for ``biorazer.access.genome_analyzer.utils.nt_utils`` (and an
    ``ensembl_rest`` package that is not a dependency), so the module could not
    be imported at all -- the helper now lives in
    :mod:`biorazer.sequence.manipulation.util`.
    """
    module = importlib.import_module("biorazer.access.database.ensembl.utils")

    assert module.get_mRNA_from_transcript is get_mRNA_from_transcript
