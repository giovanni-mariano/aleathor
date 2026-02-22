# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

from pathlib import Path
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import subprocess
import os


class BuildExtWithLibrary(build_ext):
    """Custom build_ext that compiles the C library first."""

    def run(self):
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

        # Run make lib (builds core + modules + full library)
        result = subprocess.run(
            ["make", "lib", "RELEASE=1"],
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
    return Extension(
        name="_alea",
        sources=["src/aleathor/_binding/aleathor_binding.c"],
        include_dirs=["csrc/libalea/include"],
        library_dirs=["csrc/libalea/bin"],
        # Use the combined full library (core + all modules).
        # --whole-archive is needed so the linker pulls in ALL objects,
        # overriding weak stubs with strong module implementations.
        libraries=["m"],
        extra_compile_args=["-std=c11", "-O2"],
        extra_link_args=[
            "-Wl,--whole-archive",
            "csrc/libalea/bin/libalea_full.a",
            "-Wl,--no-whole-archive",
        ],
    )


setup(
    ext_modules=[get_extension()],
    cmdclass={"build_ext": BuildExtWithLibrary},
)
