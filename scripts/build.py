from __future__ import annotations

import os, shutil, sys
from pathlib import Path
from Cython.Build import cythonize
from setuptools import Extension, Distribution
from setuptools.command.build_ext import build_ext
import numpy as np


def compile_settings(platform: str) -> tuple[list[str], list[str]]:
    """(extra_compile_args, libraries) for ``platform`` (``sys.platform``).

    MSVC (Windows) takes ``/O2`` and ``/fp:fast``, not ``-O3 -ffast-math``, and
    has no separate ``libm`` -- the math functions live in the CRT.  Linking
    ``m`` there fails, so the library is only requested where it exists.
    """
    if platform == "win32":
        return ["/O2", "/fp:fast"], []
    return ["-O3", "-ffast-math"], ["m"]


COMPILE_ARGS, LIBRARIES = compile_settings(sys.platform)
LINK_ARGS: list[str] = []
INCLUDE_DIRS: list[str] = []



def build():
    cython_build_dir = Path("build") / "cython"

    extensions = [
        Extension(
            "biorazer.structure.util.geometry.sphere.fibonacci_sphere_grid",
            sources=[
                "biorazer/structure/util/geometry/sphere/fibonacci_sphere_grid.pyx"
            ],
            extra_compile_args=COMPILE_ARGS,
            extra_link_args=LINK_ARGS,
            include_dirs=INCLUDE_DIRS + [np.get_include()],
            libraries=LIBRARIES,
        ),
    ]
    ext_modules = cythonize(
        extensions,
        compiler_directives={"language_level": "3"},
        build_dir=str(cython_build_dir),
    )

    distribution = Distribution({"name": "biorazer", "ext_modules": ext_modules})
    cmd = build_ext(distribution)
    cmd.ensure_finalized()
    cmd.run()

    # Copy built extensions back to the project directory
    for output in cmd.get_outputs():
        output_path = Path(output)
        relative_extension_path = output_path.relative_to(cmd.build_lib)

        # Drop compiled artifacts of this extension that belong to another
        # interpreter.  They would be picked up by the `*.so` / `*.pyd`
        # includes in pyproject.toml and shipped inside a wheel for the wrong
        # Python -- building cp311 and then cp312 in the same checkout (what
        # cibuildwheel does) put a cp311 .so into the cp312 wheel.
        for stale in relative_extension_path.parent.glob(
            f"{relative_extension_path.name.split('.')[0]}.*"
        ):
            if stale.name != relative_extension_path.name and stale.suffix in {
                ".so",
                ".pyd",
                ".dylib",
            }:
                stale.unlink()

        shutil.copyfile(output_path, relative_extension_path)
        mode = os.stat(relative_extension_path).st_mode
        mode |= (mode & 0o444) >> 2  # Copy R bits to X
        os.chmod(relative_extension_path, mode)


if __name__ == "__main__":
    build()
