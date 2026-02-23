# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

from pathlib import Path
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import subprocess
import os
import sys


class BuildExtWithLibrary(build_ext):
    """Custom build_ext that compiles the C library first."""

    def run(self):
        # Force MinGW compiler on Windows (MSVC cannot build this project)
        if sys.platform == "win32":
            self.compiler = "mingw32"
        # Build the C library using its Makefile
        self._build_c_library()
        # Then build the Python extension
        super().run()

    def _build_c_library(self):
        """Build libalea using make."""
        print("Building libalea C library...")

        # Get absolute path to vendor directory for make
        root = Path(__file__).parent.resolve()
        vendor = root / "csrc/libalea"

        # Build command: always RELEASE=1, optionally USE_OPENMP=1
        make_args = ["make", "lib", "RELEASE=1"]

        # Default to portable builds (safe for pip installs);
        # users who want -march=native can set PORTABLE=0
        if os.environ.get("PORTABLE", "1").strip() != "0":
            make_args.append("PORTABLE=1")

        if os.environ.get("USE_OPENMP", "").strip() == "1":
            make_args.append("USE_OPENMP=1")
            print("OpenMP support enabled")

        # Run make lib (builds core + modules + full library)
        result = subprocess.run(
            make_args,
            cwd=str(vendor),
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            print("STDOUT:", result.stdout)
            print("STDERR:", result.stderr)
            raise RuntimeError(f"Failed to build C library: {result.stderr}")

        print("C library built successfully")


def get_extension():
    """Create extension with paths resolved at setup time."""
    root = Path(__file__).parent.resolve()
    lib_path = str(root / "csrc" / "libalea" / "bin" / "libalea_full.a")

    compile_args = ["-std=c11", "-O2"]

    # macOS uses Apple ld64 which doesn't support --whole-archive
    if sys.platform == "darwin":
        link_args = [f"-Wl,-force_load,{lib_path}"]
    else:
        link_args = [
            "-Wl,--whole-archive",
            lib_path,
            "-Wl,--no-whole-archive",
        ]

    if os.environ.get("USE_OPENMP", "").strip() == "1":
        compile_args.append("-fopenmp")
        link_args.append("-fopenmp")

    if sys.platform == "win32":
        # -lm is not needed on Windows (math functions are in the standard library)
        libraries = []
        # Statically link the MinGW runtime so the .pyd is self-contained.
        # Python 3.8+ ignores PATH for DLL search, so without this the
        # extension would fail to import with missing libgcc_s_seh-1.dll.
        link_args.append("-static")
    else:
        libraries = ["m"]

    return Extension(
        name="_alea",
        sources=["src/aleathor/_binding/aleathor_binding.c"],
        include_dirs=["csrc/libalea/include"],
        library_dirs=["csrc/libalea/bin"],
        # Use the combined full library (core + all modules).
        # --whole-archive / -force_load is needed so the linker pulls in ALL
        # objects, overriding weak stubs with strong module implementations.
        libraries=libraries,
        extra_compile_args=compile_args,
        extra_link_args=link_args,
    )


setup(
    ext_modules=[get_extension()],
    cmdclass={"build_ext": BuildExtWithLibrary},
)
