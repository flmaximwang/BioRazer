"""Utilities for manipulating sequence strings directly.

Operations here work on plain sequence strings (or biotite sequence objects),
outside any file-format or object-bridge concern: use
:mod:`biorazer.sequence.io` to get sequences in and out of files, and
:mod:`biorazer.sequence.bridge` to move between sequence and structure
representations.
"""


def get_mRNA_from_transcript(transcript_seq: str, exons: list[tuple[int, int]]) -> str:
    """
    Join the exons of a transcript sequence into the mature mRNA sequence.

    Parameters
    ----------
    transcript_seq : str
        The full transcript sequence the exon coordinates refer to.
    exons : list of tuple[int, int]
        0-based ``(start, end)`` positions of the exons, in the order they
        should be joined. ``end`` is **inclusive**, i.e. the Ensembl
        convention -- an Ensembl exon (1-based, counted from the transcript
        start) becomes ``(exon_start - transcript_start,
        exon_end - transcript_start)``.
    """
    return "".join(transcript_seq[start : end + 1] for start, end in exons)
