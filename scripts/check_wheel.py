"""Fail if the installed biorazer wheel is missing its compiled extension.

cibuildwheel runs this as each build leg's test command (``CIBW_TEST_COMMAND``),
so a wheel that was assembled without the extension fails the build instead of
being published.  0.9.5 is the reason this exists: its ``win_amd64`` wheels
shipped ``__init__.py`` and the ``.pyx`` and no ``.pyd``, and nothing in the
pipeline looked inside a wheel -- the test job tests the source tree on Linux,
and the repair step only handles binaries that are already there.

The check is deliberately about the *installed* package: the checkout that
built the wheel is still on disk while the tests run, so importing from there
would prove nothing.
"""

from __future__ import annotations

import importlib
import os
import pathlib
import sys

EXTENSION = "biorazer.structure.util.geometry.sphere.fibonacci_sphere_grid"
COMPILED_SUFFIXES = {".so", ".pyd", ".dylib"}


def _fail(message: str) -> int:
    print(f"FAIL: {message}", file=sys.stderr)
    return 1


def main() -> int:
    # The checkout that built the wheel is still on disk, and the venv may even
    # live inside it (the CI job creates .verify/ in the project directory), so
    # "is the path under the project?" is the wrong question.  What matters is
    # that biorazer did *not* come from this checkout's source package, and that
    # it sits in a site-packages directory.
    checkout_package = pathlib.Path(__file__).resolve().parents[1] / "biorazer"
    sys.path[:] = [p for p in sys.path if p not in ("", os.getcwd())]

    try:
        import biorazer
    except ImportError as exc:
        return _fail(
            f"cannot import biorazer ({exc}); a development environment points "
            "back at the checkout -- this check only means something against an "
            "installed wheel"
        )

    location_str = getattr(biorazer, "__file__", None)
    if not location_str:
        return _fail("biorazer has no __file__ (is it a namespace package?)")
    location = pathlib.Path(location_str).resolve()

    if location == checkout_package or checkout_package in location.parents:
        return _fail(
            f"biorazer was imported from this checkout ({location}), not from an "
            "installed wheel"
        )
    if not any(p.name in {"site-packages", "dist-packages"} for p in location.parents):
        return _fail(f"{location} is not inside a site-packages directory")

    try:
        module = importlib.import_module(EXTENSION)
    except ImportError as exc:
        return _fail(
            f"cannot import {EXTENSION} from the installed package ({exc}); "
            "the wheel was assembled without its extension"
        )
    extension_file = pathlib.Path(module.__file__).resolve()
    if extension_file.suffix not in COMPILED_SUFFIXES:
        return _fail(
            f"{EXTENSION} is not a compiled module ({extension_file}); the "
            "wheel was assembled without its extension"
        )
    if not extension_file.is_file():
        return _fail(f"{extension_file} does not exist")

    # And it has to be the only one: a wheel that also carries another
    # interpreter's artifact (a cp312 wheel holding a cp311 .so) is what the
    # build used to produce, and on Windows a stale .pyd is one build away
    # from being the only thing in the wheel.
    siblings = sorted(
        p.name for p in extension_file.parent.iterdir() if p.suffix in COMPILED_SUFFIXES
    )
    if siblings != [extension_file.name]:
        return _fail(
            f"expected exactly one compiled artifact next to {extension_file.name}, "
            f"found {siblings}"
        )

    # It has to work, not merely exist.
    import numpy as np

    points = np.asarray(module.fibonacci_sphere_grid(np.zeros(3), 1.0, 6))
    radius = np.linalg.norm(points - points.mean(axis=0), axis=1)
    if points.shape != (6, 3) or not (abs(radius - 1.0) < 0.2).all():
        return _fail(f"fibonacci_sphere_grid returned {points!r} (radii {radius})")

    print(
        f"wheel ok: biorazer {getattr(biorazer, '__version__', '?')} at {location}\n"
        f"          extension: {extension_file.name}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
