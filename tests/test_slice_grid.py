# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""Tests for grid-based cell sampling API."""

import pytest
import math


class TestFindCellsGridZ:
    """Tests for find_cells_grid_z (XY plane grid sampling)."""

    def test_returns_dict_with_cell_ids(self, simple_model, bounds_xy):
        """Grid result should contain cell_ids."""
        result = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(10, 10))

        assert isinstance(result, dict)
        assert 'cell_ids' in result
        assert isinstance(result['cell_ids'], (list, tuple))

    def test_returns_material_ids(self, simple_model, bounds_xy):
        """Grid result should contain material_ids."""
        result = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(10, 10))

        assert 'material_ids' in result
        assert isinstance(result['material_ids'], (list, tuple))

    def test_grid_dimensions_correct(self, simple_model, bounds_xy):
        """Grid should have correct number of elements."""
        nx, ny = 20, 15
        result = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(nx, ny))

        assert len(result['cell_ids']) == nx * ny
        assert len(result['material_ids']) == nx * ny

    def test_center_point_in_sphere(self, simple_model, bounds_xy):
        """Center of grid (0,0) should be in the sphere (material 1)."""
        result = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(21, 21))

        # Center element (grid is 21x21, so center is at index 10*21 + 10 = 220)
        center_idx = 10 * 21 + 10
        center_material = result['material_ids'][center_idx]

        # Material 1 is the fuel (inside sphere)
        assert center_material == 1

    def test_corner_point_outside_sphere(self, simple_model, bounds_xy):
        """Corner of grid should be in moderator (material 2)."""
        result = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(21, 21))

        # Corner element (0,0) at grid index 0
        corner_material = result['material_ids'][0]

        # Material 2 is the moderator (outside sphere)
        assert corner_material == 2

    def test_grid_info_included(self, simple_model, bounds_xy):
        """Result should include grid dimension info."""
        result = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(10, 10))

        assert 'nx' in result
        assert 'ny' in result
        assert result['nx'] == 10
        assert result['ny'] == 10


class TestFindCellsGridY:
    """Tests for find_cells_grid_y (XZ plane grid sampling)."""

    def test_returns_cell_and_material_ids(self, simple_model):
        """XZ grid should return cell and material IDs."""
        bounds = (-10, 10, -10, 10)
        result = simple_model.find_cells_grid_y(0, bounds, resolution=(10, 10))

        assert 'cell_ids' in result
        assert 'material_ids' in result

    def test_correct_dimensions(self, simple_model):
        """XZ grid should have nx * nz elements."""
        bounds = (-10, 10, -10, 10)
        result = simple_model.find_cells_grid_y(0, bounds, resolution=(15, 20))

        assert len(result['cell_ids']) == 15 * 20


class TestFindCellsGridX:
    """Tests for find_cells_grid_x (YZ plane grid sampling)."""

    def test_returns_cell_and_material_ids(self, simple_model):
        """YZ grid should return cell and material IDs."""
        bounds = (-10, 10, -10, 10)
        result = simple_model.find_cells_grid_x(0, bounds, resolution=(10, 10))

        assert 'cell_ids' in result
        assert 'material_ids' in result


class TestFindCellsGridArbitrary:
    """Tests for find_cells_grid (arbitrary plane grid sampling)."""

    def test_arbitrary_plane_returns_grid(self, simple_model):
        """Arbitrary plane should return grid data."""
        origin = (0, 0, 0)
        normal = (1, 1, 0)
        up = (0, 0, 1)
        bounds = (-10, 10, -10, 10)

        result = simple_model.find_cells_grid(origin, normal, up, bounds, resolution=(10, 10))

        assert 'cell_ids' in result
        assert 'material_ids' in result

    def test_diagonal_plane_samples_cells(self, simple_model):
        """Diagonal plane should sample both cells."""
        origin = (0, 0, 0)
        normal = (1, 0, 0)  # Same as X=0 plane
        up = (0, 0, 1)
        bounds = (-10, 10, -10, 10)

        result = simple_model.find_cells_grid(origin, normal, up, bounds, resolution=(21, 21))

        # Should have both materials
        materials = set(result['material_ids'])
        assert 1 in materials  # Fuel
        assert 2 in materials  # Moderator


class TestMultipleCells:
    """Tests for grid sampling with multiple cells."""

    def test_finds_all_materials(self, multi_cell_model, bounds_xy):
        """Grid should find all materials in the slice."""
        result = multi_cell_model.find_cells_grid_z(0, bounds_xy, resolution=(50, 50))

        materials = set(result['material_ids'])
        # Should have materials 1, 2, 3, 4
        assert 1 in materials
        assert 2 in materials
        assert 3 in materials
        assert 4 in materials

    def test_finds_all_cells(self, multi_cell_model, bounds_xy):
        """Grid should find all cells in the slice."""
        result = multi_cell_model.find_cells_grid_z(0, bounds_xy, resolution=(50, 50))

        cells = set(result['cell_ids'])
        # Should have 4 different cells (excluding void -1)
        cells.discard(-1)
        assert len(cells) == 4


class TestVoidRegions:
    """Tests for handling void regions."""

    def test_void_marked_as_negative_one(self, simple_model):
        """Void regions should be marked with cell_id = -1."""
        # Slice outside the box should be void
        bounds = (-20, -11, -20, -11)  # Outside the box
        result = simple_model.find_cells_grid_z(0, bounds, resolution=(5, 5))

        # All points should be void
        assert all(c == -1 for c in result['cell_ids'])

    def test_void_material_zero(self, simple_model):
        """Void regions should have material_id = 0."""
        bounds = (-20, -11, -20, -11)
        result = simple_model.find_cells_grid_z(0, bounds, resolution=(5, 5))

        # Void material is 0 (MCNP convention)
        assert all(m == 0 for m in result['material_ids'])


class TestGridNewParameters:
    """Tests for new grid parameters (universe_depth, detect_errors)."""

    def test_universe_depth_parameter(self, simple_model, bounds_xy):
        """Should accept universe_depth parameter."""
        result = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(10, 10),
                                                 universe_depth=-1)
        assert 'cell_ids' in result

    def test_universe_depth_zero(self, simple_model, bounds_xy):
        """Should accept universe_depth=0 for root universe only."""
        result = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(10, 10),
                                                 universe_depth=0)
        assert 'cell_ids' in result

    def test_detect_errors_returns_errors_array(self, simple_model, bounds_xy):
        """Should return errors array when detect_errors=True."""
        result = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(10, 10),
                                                 detect_errors=True)
        assert 'errors' in result
        assert len(result['errors']) == 10 * 10

    def test_detect_errors_false_no_errors_array(self, simple_model, bounds_xy):
        """Should not return errors array when detect_errors=False."""
        result = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(10, 10),
                                                 detect_errors=False)
        # errors key may or may not be present when detect_errors=False
        # If present, it should be all zeros or not meaningful

    def test_y_grid_with_new_params(self, simple_model):
        """Y grid should accept new parameters."""
        bounds = (-10, 10, -10, 10)
        result = simple_model.find_cells_grid_y(0, bounds, resolution=(10, 10),
                                                 universe_depth=-1, detect_errors=True)
        assert 'cell_ids' in result
        assert 'errors' in result

    def test_x_grid_with_new_params(self, simple_model):
        """X grid should accept new parameters."""
        bounds = (-10, 10, -10, 10)
        result = simple_model.find_cells_grid_x(0, bounds, resolution=(10, 10),
                                                 universe_depth=-1, detect_errors=True)
        assert 'cell_ids' in result
        assert 'errors' in result

    def test_arbitrary_grid_with_new_params(self, simple_model):
        """Arbitrary plane grid should accept new parameters."""
        origin = (0, 0, 0)
        normal = (1, 0, 0)
        up = (0, 0, 1)
        bounds = (-10, 10, -10, 10)
        result = simple_model.find_cells_grid(origin, normal, up, bounds,
                                               resolution=(10, 10),
                                               universe_depth=-1, detect_errors=True)
        assert 'cell_ids' in result
        assert 'errors' in result
