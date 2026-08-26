"""Helpers for the PyRosetta Pose converters.

PyRosetta is NOT a core dependency of biorazer (it is not on PyPI), so all
imports here are LAZY: ``_import_pyrosetta`` is only called when a Pose
converter (Pdb_Pose / Cif_Pose / Pose_Pdb / Pose_Cif) is actually used, and
the informative ImportError is raised then, not at module import time.
"""

import os
import tempfile
from pathlib import Path

_PYROSETTA_INITIALIZED = False


def _import_pyrosetta():
    """Lazily import PyRosetta, raising an informative ImportError if missing."""
    try:
        import pyrosetta
    except ImportError as e:
        raise ImportError(
            "PyRosetta is required for the Pose converters (Pdb_Pose / Cif_Pose / "
            "Pose_Pdb / Pose_Cif). PyRosetta is not on PyPI; install it, e.g. via "
            "conda (-c https://conda.rosettacommons.org -c conda-forge pyrosetta) or "
            "from a pyrosetta wheel, then `pip install .[pyrosetta]`."
        ) from e
    return pyrosetta


def _ensure_pyrosetta_init():
    """Call ``pyrosetta.init`` once per process (it is process-global state)."""
    global _PYROSETTA_INITIALIZED
    if _PYROSETTA_INITIALIZED:
        return
    _import_pyrosetta().init("-mute all", silent=True)
    _PYROSETTA_INITIALIZED = True


def _pose_from_io(input_io, suffix):
    """Read a PyRosetta Pose from a str/Path (auto-detected by extension) or StringIO."""
    pyrosetta = _import_pyrosetta()
    _ensure_pyrosetta_init()
    if isinstance(input_io, (str, Path)):
        return pyrosetta.pose_from_file(str(input_io))
    # io.StringIO: PyRosetta needs a real file, so stage the text with a
    # format suffix so pose_from_file can auto-detect the format.
    with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as fh:
        fh.write(input_io.getvalue().encode("utf-8"))
        name = fh.name
    try:
        return pyrosetta.pose_from_file(name)
    finally:
        os.remove(name)


def _dump_pose(pose, output_io, suffix):
    """Write a PyRosetta Pose to a str/Path or StringIO, in PDB or mmCIF format."""
    if isinstance(output_io, (str, Path)):
        if suffix == ".pdb":
            pose.dump_pdb(str(output_io))
        else:
            pose.dump_cif(str(output_io))
        return None
    # io.StringIO: dump to a temp file, then read the text back.
    with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as fh:
        name = fh.name
    try:
        if suffix == ".pdb":
            pose.dump_pdb(name)
        else:
            pose.dump_cif(name)
        with open(name, "r") as fh:
            output_io.write(fh.read())
    finally:
        os.remove(name)
    return output_io.getvalue()
