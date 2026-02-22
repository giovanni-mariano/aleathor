# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""Tests for analytical slice curves API."""

import pytest
import math


class TestSliceCurvesZ:
    """Tests for get_slice_curves_z (XY plane slices)."""

    def test_returns_dict_with_curves(self, simple_model, bounds_xy):
        """Slice curves result should contain 'curves' key."""
        result = simple_model.get_slice_curves_z(0, bounds_xy)

        assert isinstance(result, dict)
        assert 'curves' in result
        assert isinstance(result['curves'], list)

    def test_sphere_slice_has_circle(self, simple_model, bounds_xy):
        """Slicing a sphere at z=0 should produce a circle curve."""
        result = simple_model.get_slice_curves_z(0, bounds_xy)

        # Should have at least one circle (the sphere)
        circles = [c for c in result['curves'] if c.get('type') == 'circle']
        assert len(circles) >= 1

    def test_sphere_circle_has_correct_radius(self, simple_model, bounds_xy):
        """Circle from sphere slice should have radius 5 at z=0."""
        result = simple_model.get_slice_curves_z(0, bounds_xy)

        circles = [c for c in result['curves'] if c.get('type') == 'circle']
        # Find the sphere circle (radius 5)
        sphere_circles = [c for c in circles if abs(c['radius'] - 5.0) < 0.01]
        assert len(sphere_circles) == 1

    def test_slice_above_sphere_has_no_sphere_curve(self, simple_model, bounds_xy):
        """Slicing above the sphere (z=6) should not show the sphere."""
        result = simple_model.get_slice_curves_z(6, bounds_xy)

        # Sphere radius is 5, so at z=6 there's no intersection
        circles = [c for c in result['curves'] if c.get('type') == 'circle']
        sphere_circles = [c for c in circles if c['radius'] < 5.1]
        assert len(sphere_circles) == 0

    def test_slice_near_top_has_small_circle(self, simple_model, bounds_xy):
        """Slicing near top of sphere should produce smaller circle."""
        result = simple_model.get_slice_curves_z(4, bounds_xy)

        # At z=4, sphere of radius 5 gives circle radius sqrt(25-16) = 3
        circles = [c for c in result['curves'] if c.get('type') == 'circle']
        small_circles = [c for c in circles if abs(c['radius'] - 3.0) < 0.1]
        assert len(small_circles) == 1

    def test_bounds_info_included(self, simple_model, bounds_xy):
        """Result should include bounds information."""
        result = simple_model.get_slice_curves_z(0, bounds_xy)

        assert 'x_min' in result or 'u_min' in result
        assert 'x_max' in result or 'u_max' in result


class TestSliceCurvesY:
    """Tests for get_slice_curves_y (XZ plane slices)."""

    def test_returns_curves(self, simple_model):
        """XZ slice should return curves."""
        bounds = (-10, 10, -10, 10)  # x_min, x_max, z_min, z_max
        result = simple_model.get_slice_curves_y(0, bounds)

        assert isinstance(result, dict)
        assert 'curves' in result

    def test_sphere_produces_circle(self, simple_model):
        """Sphere sliced in XZ plane should produce circle."""
        bounds = (-10, 10, -10, 10)
        result = simple_model.get_slice_curves_y(0, bounds)

        circles = [c for c in result['curves'] if c.get('type') == 'circle']
        assert len(circles) >= 1


class TestSliceCurvesX:
    """Tests for get_slice_curves_x (YZ plane slices)."""

    def test_returns_curves(self, simple_model):
        """YZ slice should return curves."""
        bounds = (-10, 10, -10, 10)  # y_min, y_max, z_min, z_max
        result = simple_model.get_slice_curves_x(0, bounds)

        assert isinstance(result, dict)
        assert 'curves' in result


class TestSliceCurvesArbitrary:
    """Tests for get_slice_curves (arbitrary plane slices)."""

    def test_arbitrary_plane_returns_curves(self, simple_model):
        """Arbitrary plane slice should return curves."""
        origin = (0, 0, 0)
        normal = (1, 1, 0)  # 45 degree diagonal
        up = (0, 0, 1)
        bounds = (-10, 10, -10, 10)

        result = simple_model.get_slice_curves(origin, normal, up, bounds)

        assert isinstance(result, dict)
        assert 'curves' in result

    def test_diagonal_slice_through_sphere(self, simple_model):
        """Diagonal slice through sphere center should produce circle."""
        origin = (0, 0, 0)
        normal = (1, 0, 0)  # Same as X slice
        up = (0, 0, 1)
        bounds = (-10, 10, -10, 10)

        result = simple_model.get_slice_curves(origin, normal, up, bounds)

        circles = [c for c in result['curves'] if c.get('type') == 'circle']
        assert len(circles) >= 1

    def test_tilted_plane(self, simple_model):
        """Tilted plane should still produce valid curves."""
        origin = (0, 0, 0)
        tilt = math.radians(30)
        normal = (0, math.sin(tilt), math.cos(tilt))
        up = (0, math.cos(tilt), -math.sin(tilt))
        bounds = (-10, 10, -10, 10)

        result = simple_model.get_slice_curves(origin, normal, up, bounds)

        assert isinstance(result, dict)
        assert 'curves' in result


class TestCylinderSlices:
    """Tests for slicing cylinders."""

    def test_z_slice_cylinder_produces_circle(self, cylinder_model, bounds_xy):
        """Z slice of vertical cylinder should produce circle."""
        result = cylinder_model.get_slice_curves_z(0, bounds_xy)

        circles = [c for c in result['curves'] if c.get('type') == 'circle']
        # Should have cylinder circle with radius 3
        cyl_circles = [c for c in circles if abs(c['radius'] - 3.0) < 0.1]
        assert len(cyl_circles) == 1

    def test_y_slice_cylinder_produces_lines(self, cylinder_model):
        """Y slice of vertical cylinder should produce lines, not circles."""
        bounds = (-10, 10, -5, 5)
        result = cylinder_model.get_slice_curves_y(0, bounds)

        # XZ slice of Z-cylinder produces two vertical lines (as parallel_lines type)
        lines = [c for c in result['curves'] if c.get('type') in ('line', 'line_segment', 'parallel_lines')]
        assert len(lines) >= 1  # parallel_lines represents two lines in one curve


class TestMultipleSurfaces:
    """Tests for models with multiple surfaces."""

    def test_concentric_spheres_produce_multiple_circles(self, multi_cell_model, bounds_xy):
        """Concentric spheres should produce multiple circles."""
        result = multi_cell_model.get_slice_curves_z(0, bounds_xy)

        circles = [c for c in result['curves'] if c.get('type') == 'circle']
        # Should have 3 circles (radii 2, 4, 6)
        assert len(circles) >= 3

    def test_circles_have_correct_radii(self, multi_cell_model, bounds_xy):
        """Circles should have correct radii."""
        result = multi_cell_model.get_slice_curves_z(0, bounds_xy)

        circles = [c for c in result['curves'] if c.get('type') == 'circle']
        radii = sorted([c['radius'] for c in circles])

        # Check for radii 2, 4, 6 (with tolerance)
        expected = [2.0, 4.0, 6.0]
        for exp in expected:
            assert any(abs(r - exp) < 0.1 for r in radii), f"Missing circle with radius {exp}"
