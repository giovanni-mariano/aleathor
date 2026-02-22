# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""Pytest fixtures for aleathor tests."""

import pytest
import math


@pytest.fixture
def simple_model():
    """Create a simple model with a sphere in a box."""
    import aleathor as ath

    model = ath.Model("Test Model")

    # Inner sphere (fuel)
    sphere = ath.Sphere(0, 0, 0, radius=5.0)

    # Outer box (moderator)
    box = ath.Box(-10, 10, -10, 10, -10, 10)

    # Add cells
    model.add_cell(
        region=-sphere,
        material=1,
        density=10.0,
        name="fuel"
    )
    model.add_cell(
        region=-box & +sphere,
        material=2,
        density=1.0,
        name="moderator"
    )

    return model


@pytest.fixture
def multi_cell_model():
    """Create a model with multiple concentric spheres."""
    import aleathor as ath

    model = ath.Model("Multi-cell Model")

    # Concentric spheres
    s1 = ath.Sphere(0, 0, 0, radius=2.0)
    s2 = ath.Sphere(0, 0, 0, radius=4.0)
    s3 = ath.Sphere(0, 0, 0, radius=6.0)
    box = ath.Box(-10, 10, -10, 10, -10, 10)

    model.add_cell(region=-s1, material=1, density=10.0, name="core")
    model.add_cell(region=-s2 & +s1, material=2, density=8.0, name="inner_shell")
    model.add_cell(region=-s3 & +s2, material=3, density=6.0, name="outer_shell")
    model.add_cell(region=-box & +s3, material=4, density=1.0, name="exterior")

    return model


@pytest.fixture
def cylinder_model():
    """Create a model with a cylinder for testing non-spherical surfaces."""
    import aleathor as ath

    model = ath.Model("Cylinder Model")

    cyl = ath.CylinderZ(0, 0, radius=3.0)
    box = ath.Box(-10, 10, -10, 10, -5, 5)

    model.add_cell(region=-cyl & -box, material=1, density=5.0, name="cylinder")
    model.add_cell(region=-box & +cyl, material=2, density=1.0, name="exterior")

    return model


@pytest.fixture
def bounds_xy():
    """Standard XY bounds for testing."""
    return (-10.0, 10.0, -10.0, 10.0)


@pytest.fixture
def bounds_3d():
    """Standard 3D bounds for testing."""
    return (-10.0, 10.0, -10.0, 10.0, -10.0, 10.0)
