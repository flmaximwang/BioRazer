"""Generic IO helpers shared by the optional-dependency converters.

These deal with the two concrete output targets a ``Converter`` may write to:
a ``str``/``Path`` file path or an ``io.StringIO``.
"""

import io
from pathlib import Path


def _written_text(output_io):
    """Return the written text when the target is an ``io.StringIO``, else None."""
    if isinstance(output_io, io.StringIO):
        return output_io.getvalue()
    return None


def _io_target(output_io):
    """Bio.PDB PDBIO/MMCIFIO accept a str filename or a file object (not a Path)."""
    return str(output_io) if isinstance(output_io, Path) else output_io
