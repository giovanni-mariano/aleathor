# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""Tests for plotting functions."""

import pytest
import math

# Skip all tests if matplotlib not available
pytest.importorskip("matplotlib")


class TestPlotSliceCurves:
    """Tests for plot_slice_curves function."""

    def test_returns_axes(self, simple_model, bounds_xy):
        """Should return matplotlib Axes object."""
        import matplotlib.pyplot as plt
        from aleathor import plot_slice_curves

        curves = simple_model.get_slice_curves_z(0, bounds_xy)
        ax = plot_slice_curves(curves)

        assert ax is not None
        plt.close('all')

    def test_accepts_existing_axes(self, simple_model, bounds_xy):
        """Should work with provided axes."""
        import matplotlib.pyplot as plt
        from aleathor import plot_slice_curves

        fig, ax = plt.subplots()
        curves = simple_model.get_slice_curves_z(0, bounds_xy)
        result_ax = plot_slice_curves(curves, ax=ax)

        assert result_ax is ax
        plt.close('all')

    def test_sets_title(self, simple_model, bounds_xy):
        """Should set title when provided."""
        import matplotlib.pyplot as plt
        from aleathor import plot_slice_curves

        curves = simple_model.get_slice_curves_z(0, bounds_xy)
        ax = plot_slice_curves(curves, title="Test Title")

        assert ax.get_title() == "Test Title"
        plt.close('all')

    def test_sets_labels(self, simple_model, bounds_xy):
        """Should set axis labels when provided."""
        import matplotlib.pyplot as plt
        from aleathor import plot_slice_curves

        curves = simple_model.get_slice_curves_z(0, bounds_xy)
        ax = plot_slice_curves(curves, xlabel="X", ylabel="Y")

        assert ax.get_xlabel() == "X"
        assert ax.get_ylabel() == "Y"
        plt.close('all')

    def test_custom_color(self, simple_model, bounds_xy):
        """Should accept custom color parameter."""
        import matplotlib.pyplot as plt
        from aleathor import plot_slice_curves

        curves = simple_model.get_slice_curves_z(0, bounds_xy)
        # Should not raise
        ax = plot_slice_curves(curves, color='red')
        plt.close('all')

    def test_custom_linewidth(self, simple_model, bounds_xy):
        """Should accept custom linewidth parameter."""
        import matplotlib.pyplot as plt
        from aleathor import plot_slice_curves

        curves = simple_model.get_slice_curves_z(0, bounds_xy)
        ax = plot_slice_curves(curves, linewidth=2.0)
        plt.close('all')


class TestPlotSliceFilled:
    """Tests for plot_slice_filled function."""

    def test_returns_axes(self, simple_model, bounds_xy):
        """Should return matplotlib Axes object."""
        import matplotlib.pyplot as plt
        from aleathor import plot_slice_filled

        grid = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(20, 20))
        ax = plot_slice_filled(grid)

        assert ax is not None
        plt.close('all')

    def test_show_contours_option(self, simple_model, bounds_xy):
        """Should draw grid-based contours when show_contours=True."""
        import matplotlib.pyplot as plt
        from aleathor import plot_slice_filled

        grid = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(20, 20))
        ax = plot_slice_filled(grid, show_contours=True)

        plt.close('all')

    def test_hide_contours_option(self, simple_model, bounds_xy):
        """Should hide contours when show_contours=False."""
        import matplotlib.pyplot as plt
        from aleathor import plot_slice_filled

        grid = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(20, 20))
        ax = plot_slice_filled(grid, show_contours=False)

        plt.close('all')

    def test_show_fill_option(self, simple_model, bounds_xy):
        """Should support show_fill parameter for contour-only mode."""
        import matplotlib.pyplot as plt
        from aleathor import plot_slice_filled

        grid = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(20, 20))
        ax = plot_slice_filled(grid, show_fill=False, show_contours=True)

        plt.close('all')

    def test_by_material_option(self, simple_model, bounds_xy):
        """Should support coloring by material."""
        import matplotlib.pyplot as plt
        from aleathor import plot_slice_filled

        grid = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(20, 20))
        ax = plot_slice_filled(grid, by_material=True)

        plt.close('all')

    def test_colorbar_option(self, simple_model, bounds_xy):
        """Should support showing colorbar."""
        import matplotlib.pyplot as plt
        from aleathor import plot_slice_filled

        grid = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(20, 20))
        ax = plot_slice_filled(grid, show_colorbar=True)

        plt.close('all')

    def test_custom_cmap(self, simple_model, bounds_xy):
        """Should accept custom colormap."""
        import matplotlib.pyplot as plt
        from aleathor import plot_slice_filled

        grid = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(20, 20))
        ax = plot_slice_filled(grid, cmap='viridis')

        plt.close('all')


class TestPlotCurve:
    """Tests for plot_curve function."""

    def test_plots_circle(self):
        """Should plot circle curve."""
        import matplotlib.pyplot as plt
        from aleathor import plot_curve

        fig, ax = plt.subplots()
        curve = {'type': 'circle', 'center': (0, 0), 'radius': 5}
        plot_curve(curve, ax, color='blue')

        plt.close('all')

    def test_plots_line(self):
        """Should plot line curve."""
        import matplotlib.pyplot as plt
        from aleathor import plot_curve

        fig, ax = plt.subplots()
        curve = {'type': 'line', 'point': (0, 0), 'direction': (1, 0)}
        plot_curve(curve, ax, color='blue')

        plt.close('all')

    def test_plots_ellipse(self):
        """Should plot ellipse curve."""
        import matplotlib.pyplot as plt
        from aleathor import plot_curve

        fig, ax = plt.subplots()
        curve = {'type': 'ellipse', 'center': (0, 0), 'semi_a': 5, 'semi_b': 3}
        plot_curve(curve, ax, color='blue')

        plt.close('all')


class TestPlotRayPath:
    """Tests for plot_ray_path function."""

    def test_returns_axes(self, simple_model):
        """Should return matplotlib Axes object."""
        import matplotlib.pyplot as plt
        from aleathor import plot_ray_path

        # Get ray trace segments
        trace = simple_model.trace(origin=(-10, 0, 0), direction=(1, 0, 0))
        segments = [
            {
                't_enter': seg.t_enter,
                't_exit': seg.t_exit,
                'cell_id': seg.cell.id if seg.cell else -1,
                'material_id': seg.material
            }
            for seg in trace
        ]

        ax = plot_ray_path(segments, t_max=20)

        assert ax is not None
        plt.close('all')


class TestContourBy:
    """Tests for contour_by parameter in plot_slice_filled and Model.plot."""

    def test_contour_by_cell(self, simple_model, bounds_xy):
        """Should accept contour_by='cell' without error."""
        import matplotlib.pyplot as plt
        from aleathor import plot_slice_filled

        grid = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(20, 20))
        ax = plot_slice_filled(grid, contour_by='cell')
        assert ax is not None
        plt.close('all')

    def test_contour_by_material(self, simple_model, bounds_xy):
        """Should accept contour_by='material' without error."""
        import matplotlib.pyplot as plt
        from aleathor import plot_slice_filled

        grid = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(20, 20))
        ax = plot_slice_filled(grid, contour_by='material')
        assert ax is not None
        plt.close('all')

    def test_auto_default_material_when_by_material(self, simple_model, bounds_xy):
        """contour_by=None should auto-select 'material' when by_material=True."""
        import matplotlib.pyplot as plt
        import numpy as np
        from aleathor import plot_slice_filled

        grid = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(20, 20))
        # by_material=True with default contour_by should not raise
        ax = plot_slice_filled(grid, by_material=True)
        assert ax is not None
        plt.close('all')

    def test_explicit_override_cell_with_by_material(self, simple_model, bounds_xy):
        """Should allow contour_by='cell' even when by_material=True."""
        import matplotlib.pyplot as plt
        from aleathor import plot_slice_filled

        grid = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(20, 20))
        ax = plot_slice_filled(grid, by_material=True, contour_by='cell')
        assert ax is not None
        plt.close('all')

    def test_invalid_contour_by_raises(self, simple_model, bounds_xy):
        """Should raise ValueError for invalid contour_by value."""
        from aleathor import plot_slice_filled

        grid = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(20, 20))
        with pytest.raises(ValueError, match="contour_by"):
            plot_slice_filled(grid, contour_by='invalid')

    def test_overlay_defaults_to_material_contours(self, simple_model, bounds_xy):
        """contour_by=None should auto-select 'material' when overlay is set."""
        import matplotlib.pyplot as plt
        import numpy as np
        from aleathor import plot_slice_filled

        grid = simple_model.find_cells_grid_z(0, bounds_xy, resolution=(20, 20))
        overlay = np.random.rand(20, 20)
        ax = plot_slice_filled(grid, overlay=overlay)
        assert ax is not None
        plt.close('all')

    def test_model_plot_contour_by(self, simple_model, bounds_xy):
        """Model.plot() should accept and forward contour_by parameter."""
        import matplotlib.pyplot as plt

        ax = simple_model.plot(z=0, bounds=bounds_xy, contour_by='material')
        assert ax is not None
        plt.close('all')


class TestModelPlot:
    """Tests for Model.plot() method."""

    def test_z_slice(self, simple_model, bounds_xy):
        """Should plot XY slice at z=constant."""
        import matplotlib.pyplot as plt

        ax = simple_model.plot(z=0, bounds=bounds_xy)
        assert ax is not None
        plt.close('all')

    def test_y_slice(self, simple_model):
        """Should plot XZ slice at y=constant."""
        import matplotlib.pyplot as plt

        bounds = (-10, 10, -10, 10)
        ax = simple_model.plot(y=0, bounds=bounds)
        assert ax is not None
        plt.close('all')

    def test_x_slice(self, simple_model):
        """Should plot YZ slice at x=constant."""
        import matplotlib.pyplot as plt

        bounds = (-10, 10, -10, 10)
        ax = simple_model.plot(x=0, bounds=bounds)
        assert ax is not None
        plt.close('all')

    def test_arbitrary_plane(self, simple_model):
        """Should plot arbitrary plane slice."""
        import matplotlib.pyplot as plt

        ax = simple_model.plot(
            origin=(0, 0, 0),
            normal=(1, 1, 0),
            up=(0, 0, 1),
            bounds=(-10, 10, -10, 10)
        )
        assert ax is not None
        plt.close('all')

    def test_filled_false(self, simple_model, bounds_xy):
        """Should plot curves only when filled=False."""
        import matplotlib.pyplot as plt

        ax = simple_model.plot(z=0, bounds=bounds_xy, filled=False)
        assert ax is not None
        plt.close('all')

    def test_by_material(self, simple_model, bounds_xy):
        """Should color by material when by_material=True."""
        import matplotlib.pyplot as plt

        ax = simple_model.plot(z=0, bounds=bounds_xy, by_material=True)
        assert ax is not None
        plt.close('all')

    def test_default_z_zero(self, simple_model, bounds_xy):
        """Should default to z=0 when no plane specified."""
        import matplotlib.pyplot as plt

        ax = simple_model.plot(bounds=bounds_xy)
        assert ax is not None
        plt.close('all')

    def test_save_to_file(self, simple_model, bounds_xy, tmp_path):
        """Should save plot to file."""
        import matplotlib.pyplot as plt

        filepath = tmp_path / "test_plot.png"
        simple_model.plot(z=0, bounds=bounds_xy, save=str(filepath))

        assert filepath.exists()
        plt.close('all')

    def test_show_contours_parameter(self, simple_model, bounds_xy):
        """Should support show_contours parameter."""
        import matplotlib.pyplot as plt

        ax = simple_model.plot(z=0, bounds=bounds_xy, show_contours=False)
        assert ax is not None
        plt.close('all')

    def test_detect_errors_parameter(self, simple_model, bounds_xy):
        """Should support detect_errors parameter."""
        import matplotlib.pyplot as plt

        ax = simple_model.plot(z=0, bounds=bounds_xy, detect_errors=True)
        assert ax is not None
        plt.close('all')

    def test_universe_depth_parameter(self, simple_model, bounds_xy):
        """Should support universe_depth parameter."""
        import matplotlib.pyplot as plt

        ax = simple_model.plot(z=0, bounds=bounds_xy, universe_depth=0)
        assert ax is not None
        plt.close('all')


class TestModelPlotViews:
    """Tests for Model.plot_views() method."""

    def test_returns_figure(self, simple_model, bounds_3d):
        """Should return matplotlib Figure."""
        import matplotlib.pyplot as plt

        fig = simple_model.plot_views(bounds=bounds_3d)
        assert fig is not None
        plt.close('all')

    def test_creates_three_subplots(self, simple_model, bounds_3d):
        """Should create figure with 3 subplots."""
        import matplotlib.pyplot as plt

        fig = simple_model.plot_views(bounds=bounds_3d)
        assert len(fig.axes) == 3
        plt.close('all')

    def test_save_to_file(self, simple_model, bounds_3d, tmp_path):
        """Should save views to file."""
        import matplotlib.pyplot as plt

        filepath = tmp_path / "test_views.png"
        simple_model.plot_views(bounds=bounds_3d, save=str(filepath))

        assert filepath.exists()
        plt.close('all')
