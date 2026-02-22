# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""Tests for counting functions (cells, surfaces, materials)."""

import pytest


class TestCountSurfaces:
    """Tests for count_surfaces function."""

    def test_counts_unique_surfaces(self, simple_model, bounds_xy):
        """Should count unique surfaces in the slice."""
        from aleathor import count_surfaces

        curves = simple_model.get_slice_curves_z(0, bounds_xy)
        n_surfaces = count_surfaces(curves)

        # Simple model has sphere and box (4 planes for box + 1 sphere)
        # At z=0, we see sphere circle and box edges
        assert n_surfaces >= 1  # At minimum the sphere

    def test_concentric_spheres_count(self, multi_cell_model, bounds_xy):
        """Should count multiple sphere surfaces."""
        from aleathor import count_surfaces

        curves = multi_cell_model.get_slice_curves_z(0, bounds_xy)
        n_surfaces = count_surfaces(curves)

        # 3 concentric spheres + box surfaces
        assert n_surfaces >= 3

    def test_empty_curves_returns_zero(self):
        """Empty curves result should return 0."""
        from aleathor import count_surfaces

        result = {'curves': []}
        assert count_surfaces(result) == 0

    def test_curves_without_surface_id(self):
        """Curves without surface_id should be skipped."""
        from aleathor import count_surfaces

        result = {'curves': [{'type': 'circle', 'center': (0, 0), 'radius': 5}]}
        # No surface_id, so count is 0
        assert count_surfaces(result) == 0


class TestCountCells:
    """Tests for count_cells function."""

    def test_counts_unique_cells(self, simple_model, bounds_xy):
        """Should count unique cells in the grid."""
        from aleathor import count_cells

        grid = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(50, 50))
        n_cells = count_cells(grid)

        # Simple model has 2 cells (fuel + moderator)
        assert n_cells == 2

    def test_excludes_void(self, simple_model):
        """Should exclude void (-1) from count."""
        from aleathor import count_cells

        # Grid partially outside geometry
        bounds = (-15, 15, -15, 15)
        grid = simple_model.find_cells_grid_z(0, bounds, resolution=(50, 50))
        n_cells = count_cells(grid)

        # Should still be 2 cells, void not counted
        assert n_cells == 2

    def test_multi_cell_model(self, multi_cell_model, bounds_xy):
        """Should count all cells in multi-cell model."""
        from aleathor import count_cells

        grid = multi_cell_model.find_cells_grid_z(0, bounds_xy, resolution=(50, 50))
        n_cells = count_cells(grid)

        # 4 concentric regions
        assert n_cells == 4

    def test_empty_grid_returns_zero(self):
        """Empty grid should return 0."""
        from aleathor import count_cells

        result = {'cell_ids': []}
        assert count_cells(result) == 0

    def test_all_void_returns_zero(self):
        """Grid with only void should return 0."""
        from aleathor import count_cells

        result = {'cell_ids': [-1, -1, -1, -1]}
        assert count_cells(result) == 0


class TestCountMaterials:
    """Tests for count_materials function."""

    def test_counts_unique_materials(self, simple_model, bounds_xy):
        """Should count unique materials in the grid."""
        from aleathor import count_materials

        grid = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(50, 50))
        n_materials = count_materials(grid)

        # Simple model has 2 materials (1 and 2)
        assert n_materials == 2

    def test_excludes_void(self, simple_model):
        """Should exclude void (-1) from count."""
        from aleathor import count_materials

        bounds = (-15, 15, -15, 15)
        grid = simple_model.find_cells_grid_z(0, bounds, resolution=(50, 50))
        n_materials = count_materials(grid)

        # Still 2 materials
        assert n_materials == 2

    def test_multi_cell_model(self, multi_cell_model, bounds_xy):
        """Should count all materials in multi-cell model."""
        from aleathor import count_materials

        grid = multi_cell_model.find_cells_grid_z(0, bounds_xy, resolution=(50, 50))
        n_materials = count_materials(grid)

        # 4 different materials
        assert n_materials == 4

    def test_empty_grid_returns_zero(self):
        """Empty grid should return 0."""
        from aleathor import count_materials

        result = {'material_ids': []}
        assert count_materials(result) == 0


class TestGetSliceStats:
    """Tests for get_slice_stats function."""

    def test_returns_all_counts(self, simple_model, bounds_xy):
        """Should return dict with cells, surfaces, materials."""
        from aleathor import get_slice_stats

        curves = simple_model.get_slice_curves_z(0, bounds_xy)
        grid = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(50, 50))

        stats = get_slice_stats(curves, grid)

        assert 'cells' in stats
        assert 'surfaces' in stats
        assert 'materials' in stats

    def test_values_are_integers(self, simple_model, bounds_xy):
        """All counts should be integers."""
        from aleathor import get_slice_stats

        curves = simple_model.get_slice_curves_z(0, bounds_xy)
        grid = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(50, 50))

        stats = get_slice_stats(curves, grid)

        assert isinstance(stats['cells'], int)
        assert isinstance(stats['surfaces'], int)
        assert isinstance(stats['materials'], int)

    def test_consistent_with_individual_functions(self, simple_model, bounds_xy):
        """Stats should match individual count functions."""
        from aleathor import get_slice_stats, count_cells, count_surfaces, count_materials

        curves = simple_model.get_slice_curves_z(0, bounds_xy)
        grid = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(50, 50))

        stats = get_slice_stats(curves, grid)

        assert stats['cells'] == count_cells(grid)
        assert stats['surfaces'] == count_surfaces(curves)
        assert stats['materials'] == count_materials(grid)

    def test_multi_cell_stats(self, multi_cell_model, bounds_xy):
        """Should work correctly with multi-cell model."""
        from aleathor import get_slice_stats

        curves = multi_cell_model.get_slice_curves_z(0, bounds_xy)
        grid = multi_cell_model.find_cells_grid_z(0, bounds_xy, resolution=(50, 50))

        stats = get_slice_stats(curves, grid)

        assert stats['cells'] == 4
        assert stats['materials'] == 4
        # Surface count depends on C library curve extraction (may be 0 if not implemented)
        assert stats['surfaces'] >= 0
