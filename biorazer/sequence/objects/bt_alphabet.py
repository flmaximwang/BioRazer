"""biotite's alphabet objects, re-exported for the object layer.

An ``Alphabet`` is the symbol set a sequence is defined over -- the letters,
the gaps, and any ambiguity codes -- so a sequence object is meaningless
without it: ``ProteinSequence.alphabet`` **is** the standard 20-letter
``LetterAlphabet``, and a profile carries the alphabet its columns count.

``LetterAlphabet`` is the concrete alphabet over 1-character symbols that
biorazer builds everywhere (the amino-acid families in
:mod:`biorazer.database.alphabet.protein`, the alphabet-with-gap variants in
the ColabFold MSA pipeline).  ``AlphabetError`` is the exception biotite raises
when a letter is not in the alphabet -- raised by ``ProteinSequence(...)`` for
an invalid residue, so it is imported from here wherever such a construction is
caught.

The alphabets biorazer itself defines (``COMMON_AA_ALPHABET`` and friends) are
*data*, not classes, and stay in :mod:`biorazer.database.alphabet`.
"""

from biotite.sequence import AlphabetError, LetterAlphabet

__all__ = [
    "LetterAlphabet",
    "AlphabetError",
]
