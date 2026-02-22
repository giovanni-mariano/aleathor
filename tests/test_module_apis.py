# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""Tests for module-level function APIs (slicing helpers, plotting).

These test that the slicing DRY helper and plotting module-level functions
work correctly.
"""

import pytest


# =========================================================================
# slicing module
# =========================================================================

class TestSlicingImport:
    """Test that slicing functions are importable."""

    def test_import_module(self):
        from aleathor import slicing
        assert slicing is not None

    def test_all_functions_exist(self):
        from aleathor import slicing
        for name in [
            'find_label_positions', 'find_surface_label_positions',
            'check_grid_overlaps', '_extract_slice_params',
        ]:
            assert hasattr(slicing, name), f"Missing: slicing.{name}"


class TestSlicingFindLabelPositions:
    """Test slicing.find_label_positions."""

    def test_returns_list(self, simple_model, bounds_xy):
        from aleathor.slicing import find_label_positions
        grid = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(50, 50))
        labels = find_label_positions(simple_model, grid, min_pixels=10)
        assert isinstance(labels, list)

    def test_matches_method(self, simple_model, bounds_xy):
        from aleathor.slicing import find_label_positions
        grid = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(50, 50))
        r1 = find_label_positions(simple_model, grid, min_pixels=10)
        r2 = simple_model.find_label_positions(grid, min_pixels=10)
        # Compare deterministic fields (pixel positions may vary between calls)
        ids1 = sorted(d['id'] for d in r1)
        ids2 = sorted(d['id'] for d in r2)
        assert ids1 == ids2
        counts1 = sorted(d['pixel_count'] for d in r1)
        counts2 = sorted(d['pixel_count'] for d in r2)
        assert counts1 == counts2

    def test_by_material(self, simple_model, bounds_xy):
        from aleathor.slicing import find_label_positions
        grid = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(50, 50))
        labels = find_label_positions(simple_model, grid, by_material=True)
        assert isinstance(labels, list)


class TestSlicingFindSurfaceLabelPositions:
    """Test slicing.find_surface_label_positions."""

    def test_returns_list(self, simple_model, bounds_xy):
        from aleathor.slicing import find_surface_label_positions
        grid = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(50, 50))
        labels = find_surface_label_positions(simple_model, grid)
        assert isinstance(labels, list)

    def test_matches_method(self, simple_model, bounds_xy):
        from aleathor.slicing import find_surface_label_positions
        grid = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(50, 50))
        r1 = find_surface_label_positions(simple_model, grid, margin=20)
        r2 = simple_model.find_surface_label_positions(grid, margin=20)
        assert r1 == r2


class TestSlicingCheckGridOverlaps:
    """Test slicing.check_grid_overlaps."""

    def test_returns_list(self, simple_model, bounds_xy):
        from aleathor.slicing import check_grid_overlaps
        grid = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(10, 10),
                                               detect_errors=True)
        errors = check_grid_overlaps(simple_model, grid)
        assert isinstance(errors, list)
        assert len(errors) == 100

    def test_matches_method(self, simple_model, bounds_xy):
        from aleathor.slicing import check_grid_overlaps
        grid = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(10, 10),
                                               detect_errors=True)
        r1 = check_grid_overlaps(simple_model, grid)
        r2 = simple_model.check_grid_overlaps(grid)
        assert r1 == r2


class TestExtractSliceParams:
    """Test slicing._extract_slice_params DRY helper."""

    def test_z_grid(self, simple_model, bounds_xy):
        from aleathor.slicing import _extract_slice_params
        grid = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(10, 15))
        nu, nv, origin, normal, up, u_min, u_max, v_min, v_max = _extract_slice_params(grid)
        assert nu == 10
        assert nv == 15
        assert normal == (0, 0, 1)

    def test_y_grid(self, simple_model):
        from aleathor.slicing import _extract_slice_params
        bounds = (-10, 10, -10, 10)
        grid = simple_model.find_cells_grid_y(0, bounds, resolution=(12, 8))
        nu, nv, origin, normal, up, u_min, u_max, v_min, v_max = _extract_slice_params(grid)
        assert nu == 12
        assert nv == 8
        assert normal == (0, 1, 0)

    def test_x_grid(self, simple_model):
        from aleathor.slicing import _extract_slice_params
        bounds = (-10, 10, -10, 10)
        grid = simple_model.find_cells_grid_x(0, bounds, resolution=(7, 9))
        nu, nv, origin, normal, up, u_min, u_max, v_min, v_max = _extract_slice_params(grid)
        assert nu == 7
        assert nv == 9
        assert normal == (1, 0, 0)

    def test_arbitrary_grid(self, simple_model):
        from aleathor.slicing import _extract_slice_params
        origin = (0, 0, 0)
        normal = (1, 0, 0)
        up = (0, 0, 1)
        bounds = (-10, 10, -10, 10)
        grid = simple_model.find_cells_grid(origin, normal, up, bounds,
                                             resolution=(11, 13))
        nu, nv = _extract_slice_params(grid)[:2]
        assert nu == 11
        assert nv == 13

    def test_invalid_grid_raises(self):
        from aleathor.slicing import _extract_slice_params
        with pytest.raises(ValueError, match="Cannot determine"):
            _extract_slice_params({'foo': 'bar'})


# =========================================================================
# plotting module-level functions
# =========================================================================

pytest.importorskip("matplotlib")


class TestPlottingModulePlot:
    """Test plotting.plot as module-level function."""

    def test_z_slice(self, simple_model, bounds_xy):
        import matplotlib.pyplot as plt
        from aleathor.plotting import plot
        ax = plot(simple_model, z=0, bounds=bounds_xy)
        assert ax is not None
        plt.close('all')

    def test_y_slice(self, simple_model):
        import matplotlib.pyplot as plt
        from aleathor.plotting import plot
        ax = plot(simple_model, y=0, bounds=(-10, 10, -10, 10))
        assert ax is not None
        plt.close('all')

    def test_x_slice(self, simple_model):
        import matplotlib.pyplot as plt
        from aleathor.plotting import plot
        ax = plot(simple_model, x=0, bounds=(-10, 10, -10, 10))
        assert ax is not None
        plt.close('all')

    def test_default_z_zero(self, simple_model, bounds_xy):
        import matplotlib.pyplot as plt
        from aleathor.plotting import plot
        ax = plot(simple_model, bounds=bounds_xy)
        assert ax is not None
        plt.close('all')

    def test_by_material(self, simple_model, bounds_xy):
        import matplotlib.pyplot as plt
        from aleathor.plotting import plot
        ax = plot(simple_model, z=0, bounds=bounds_xy, by_material=True)
        assert ax is not None
        plt.close('all')

    def test_save_to_file(self, simple_model, bounds_xy, tmp_path):
        import matplotlib.pyplot as plt
        from aleathor.plotting import plot
        filepath = tmp_path / "test_mod_plot.png"
        plot(simple_model, z=0, bounds=bounds_xy, save=str(filepath))
        assert filepath.exists()
        plt.close('all')

    def test_matches_method(self, simple_model, bounds_xy):
        """Module function and method should both succeed."""
        import matplotlib.pyplot as plt
        from aleathor.plotting import plot
        ax1 = plot(simple_model, z=0, bounds=bounds_xy)
        assert ax1 is not None
        plt.close('all')
        ax2 = simple_model.plot(z=0, bounds=bounds_xy)
        assert ax2 is not None
        plt.close('all')


class TestPlottingModulePlotViews:
    """Test plotting.plot_views as module-level function."""

    def test_returns_figure(self, simple_model, bounds_3d):
        import matplotlib.pyplot as plt
        from aleathor.plotting import plot_views
        fig = plot_views(simple_model, bounds=bounds_3d)
        assert fig is not None
        assert len(fig.axes) == 3
        plt.close('all')

    def test_save_to_file(self, simple_model, bounds_3d, tmp_path):
        import matplotlib.pyplot as plt
        from aleathor.plotting import plot_views
        filepath = tmp_path / "test_mod_views.png"
        plot_views(simple_model, bounds=bounds_3d, save=str(filepath))
        assert filepath.exists()
        plt.close('all')

    def test_matches_method(self, simple_model, bounds_3d):
        import matplotlib.pyplot as plt
        from aleathor.plotting import plot_views
        fig1 = plot_views(simple_model, bounds=bounds_3d)
        assert len(fig1.axes) == 3
        plt.close('all')
        fig2 = simple_model.plot_views(bounds=bounds_3d)
        assert len(fig2.axes) == 3
        plt.close('all')
