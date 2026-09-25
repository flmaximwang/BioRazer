"""PyRosetta's ``Pose``, enumerated for the object layer.

``rosetta.core.pose.Pose`` is the fourth structure object biorazer works with
(Rosetta's full-fledged pose: residues with their rotamers, fold tree, energies
and score function).  The Pose converters in
:mod:`biorazer.structure.io.protein` (``Pdb_Pose`` / ``Cif_Pose`` /
``Pose_Pdb`` / ``Pose_Cif``) return and accept instances of exactly this class.

PyRosetta is an **optional** dependency -- it is not on PyPI, so it is declared
as the ``pyrosetta`` extra in ``pyproject.toml`` -- therefore the class cannot
be imported at module level: doing so would make ``import
biorazer.structure.objects`` fail on every install without PyRosetta.  It is
enumerated as a *lazy accessor* instead: :func:`pose_class` imports and returns
the type on demand.

Nothing in biorazer annotates a Pose today (the converters get theirs from
``pyrosetta.pose_from_file``), so this accessor has no in-library caller yet.
It exists so that the object is *named* in one place instead of living only
inside a lazy import in the I/O layer.
"""

__all__ = [
    "pose_class",
]


def pose_class():
    """The PyRosetta ``Pose`` type, imported lazily.

    Returns the class object, e.g. for annotations::

        Pose = biorazer.structure.objects.pose_class()

    Raises
    ------
    ImportError
        If PyRosetta is not installed.  It is optional (not on PyPI); install
        it first, e.g. with ``conda -c
        https://conda.rosettacommons.org -c conda-forge pyrosetta``.
    """
    import pyrosetta
    return pyrosetta.Pose
