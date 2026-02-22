# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""
Matplotlib-based plotting for CSG geometry visualization.

Includes:
- Analytical curve plotting (plot_slice_curves, plot_curve)
- Grid-based filled plotting (plot_slice_filled)
- High-level plot/plot_views orchestration functions
- Cell/surface/material counting from analytical results

The ``plot`` and ``plot_views`` functions are also available as methods on
the Model class for convenience.
"""

from __future__ import annotations
from typing import Optional, Tuple, Dict, Any, Set, TYPE_CHECKING

if TYPE_CHECKING:
    from .model import Model

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None

try:
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    from matplotlib.axes import Axes
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


def plot_ray_path(segments: list,
                  ax: Optional['Axes'] = None,
                  show_materials: bool = True,
                  t_max: float = None) -> 'Axes':
    """Plot ray path through geometry as 1D bar chart.

    Args:
        segments: List of segment dicts from raycast
        ax: Matplotlib axes
        show_materials: Color by material
        t_max: Maximum t value to show (clips infinite segments)

    Returns:
        The matplotlib Axes object
    """
    if not HAS_MATPLOTLIB:
        raise ImportError("matplotlib is required for plotting")

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 2))

    # Color map for materials
    colors = plt.cm.tab10.colors

    for seg in segments:
        t_enter = seg['t_enter']
        t_exit = seg['t_exit']
        cell_id = seg['cell_id']
        mat_id = seg['material_id']

        # Clip infinite segments
        if t_max and t_exit > t_max:
            t_exit = t_max
        if t_exit > 1e30:
            continue  # Skip infinite void at end

        width = t_exit - t_enter

        if cell_id < 0:
            color = 'lightgray'
            label = 'void'
        else:
            color = colors[mat_id % len(colors)]
            label = f'mat {mat_id}'

        ax.barh(0, width, left=t_enter, height=0.8, color=color,
                edgecolor='black', linewidth=0.5)

        # Add cell label
        if cell_id >= 0 and width > 0.1:
            ax.text(t_enter + width/2, 0, str(cell_id),
                   ha='center', va='center', fontsize=8)

    ax.set_ylim(-0.5, 0.5)
    ax.set_yticks([])
    ax.set_xlabel('Distance along ray')
    ax.set_title('Ray path through geometry')

    return ax


def plot_curve(curve: Dict[str, Any], ax: 'Axes', **kwargs) -> None:
    """Plot a single curve from slice_curves result.

    Args:
        curve: Curve dict from get_slice_curves_*
        ax: Matplotlib axes
        **kwargs: Additional arguments passed to matplotlib (color, linewidth, etc.)
    """
    if not HAS_MATPLOTLIB:
        raise ImportError("matplotlib is required for plotting")

    curve_type = curve.get('type', 'none')

    if curve_type == 'circle':
        center = curve['center']
        radius = curve['radius']
        circle = plt.Circle(center, radius, fill=False, **kwargs)
        ax.add_patch(circle)

    elif curve_type == 'arc':
        from matplotlib.patches import Arc
        center = curve['center']
        radius = curve['radius']
        theta1 = np.degrees(curve.get('theta_start', 0))
        theta2 = np.degrees(curve.get('theta_end', 2 * np.pi))
        arc = Arc(center, 2*radius, 2*radius, angle=0,
                  theta1=theta1, theta2=theta2, **kwargs)
        ax.add_patch(arc)

    elif curve_type == 'ellipse':
        from matplotlib.patches import Ellipse
        center = curve['center']
        width = 2 * curve['semi_a']
        height = 2 * curve['semi_b']
        angle = np.degrees(curve.get('angle', 0))
        ellipse = Ellipse(center, width, height, angle=angle, fill=False, **kwargs)
        ax.add_patch(ellipse)

    elif curve_type in ('line', 'line_segment'):
        point = curve['point']
        direction = curve['direction']
        # For infinite lines, clip to axes limits
        t_min = curve.get('t_min', -1000)
        t_max = curve.get('t_max', 1000)
        p1 = (point[0] + t_min * direction[0], point[1] + t_min * direction[1])
        p2 = (point[0] + t_max * direction[0], point[1] + t_max * direction[1])
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], **kwargs)

    elif curve_type == 'polygon':
        vertices = curve['vertices']
        if curve.get('closed', False):
            vertices = vertices + [vertices[0]]  # Close the polygon
        xs = [v[0] for v in vertices]
        ys = [v[1] for v in vertices]
        ax.plot(xs, ys, **kwargs)

    elif curve_type == 'parallel_lines':
        # Two parallel lines from point1 and point2 in the given direction
        point1 = curve['point1']
        point2 = curve['point2']
        direction = curve['direction']
        t_min = curve.get('t_min', -1000)
        t_max = curve.get('t_max', 1000)
        # Line 1
        p1_start = (point1[0] + t_min * direction[0], point1[1] + t_min * direction[1])
        p1_end = (point1[0] + t_max * direction[0], point1[1] + t_max * direction[1])
        ax.plot([p1_start[0], p1_end[0]], [p1_start[1], p1_end[1]], **kwargs)
        # Line 2
        p2_start = (point2[0] + t_min * direction[0], point2[1] + t_min * direction[1])
        p2_end = (point2[0] + t_max * direction[0], point2[1] + t_max * direction[1])
        ax.plot([p2_start[0], p2_end[0]], [p2_start[1], p2_end[1]], **kwargs)


def plot_slice_curves(curves_result: Dict[str, Any],
                      ax: Optional['Axes'] = None,
                      color: str = 'black',
                      linewidth: float = 0.5,
                      title: Optional[str] = None,
                      xlabel: Optional[str] = None,
                      ylabel: Optional[str] = None,
                      xlim: Optional[Tuple[float, float]] = None,
                      ylim: Optional[Tuple[float, float]] = None) -> 'Axes':
    """Plot slice curves (vector graphics).

    Args:
        curves_result: Result from sys.get_slice_curves_z/y/x
        ax: Matplotlib axes (creates new if None)
        color: Curve color
        linewidth: Line width
        title, xlabel, ylabel: Labels
        xlim, ylim: Optional axis limits (defaults to viewport from curves_result)

    Returns:
        Matplotlib Axes
    """
    if not HAS_MATPLOTLIB:
        raise ImportError("matplotlib is required for plotting")
    if not HAS_NUMPY:
        raise ImportError("numpy is required for plotting")

    if ax is None:
        fig, ax = plt.subplots()

    for curve in curves_result.get('curves', []):
        plot_curve(curve, ax, color=color, linewidth=linewidth)

    # Set limits - use provided or from result
    # Determine which plane this is based on available keys:
    # Z slice: x_min/max, y_min/max
    # Y slice: x_min/max, z_min/max
    # X slice: y_min/max, z_min/max
    if xlim:
        ax.set_xlim(xlim)
    else:
        # Horizontal axis: x for Z/Y slices, y for X slice
        if 'x_min' in curves_result:
            ax.set_xlim(curves_result['x_min'], curves_result['x_max'])
        elif 'y_min' in curves_result:
            ax.set_xlim(curves_result['y_min'], curves_result['y_max'])
        else:
            ax.set_xlim(curves_result.get('u_min', -10), curves_result.get('u_max', 10))

    if ylim:
        ax.set_ylim(ylim)
    else:
        # Vertical axis: y for Z slice, z for Y/X slices
        if 'y_min' in curves_result and 'x_min' in curves_result:
            # Z slice (has both x and y)
            ax.set_ylim(curves_result['y_min'], curves_result['y_max'])
        elif 'z_min' in curves_result:
            # Y or X slice (has z)
            ax.set_ylim(curves_result['z_min'], curves_result['z_max'])
        else:
            ax.set_ylim(curves_result.get('v_min', -10), curves_result.get('v_max', 10))

    ax.set_aspect('equal')

    if title:
        ax.set_title(title)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)

    return ax


def plot_slice_filled(grid_result: Dict[str, Any],
                      ax: Optional['Axes'] = None,
                      cmap: str = 'tab20',
                      by_material: bool = False,
                      show_colorbar: bool = False,
                      show_fill: bool = True,
                      show_contours: bool = True,
                      contour_color: str = 'black',
                      contour_width: float = 0.5,
                      show_errors: bool = True,
                      overlap_color: str = 'red',
                      undefined_color: str = 'orange',
                      error_contour_style: str = 'dashed',
                      title: Optional[str] = None,
                      xlabel: Optional[str] = None,
                      ylabel: Optional[str] = None,
                      overlay: Optional[Any] = None,
                      overlay_extent: Optional[Tuple[float, float, float, float]] = None,
                      overlay_cmap: str = 'jet',
                      overlay_norm: Optional[str] = None,
                      overlay_alpha: float = 1.0,
                      overlay_vmin: Optional[float] = None,
                      overlay_vmax: Optional[float] = None,
                      overlay_label: Optional[str] = None,
                      contour_by: Optional[str] = None) -> 'Axes':
    """Plot slice with grid-based contours (recommended approach).

    Contours are detected by finding boundaries between adjacent pixels with
    different cell IDs, similar to plotter.c. This is faster and more accurate
    than using analytical curves for contour display.

    Args:
        grid_result: Result from sys.find_cells_grid_z/y/x (with detect_errors=True for error display)
        ax: Matplotlib axes
        cmap: Colormap name for geometry fill
        by_material: Color by material ID instead of cell ID
        show_colorbar: Show colorbar
        show_fill: If True, show filled regions; if False, show only contours
        show_contours: If True, draw boundary contours between cells
        contour_color: Color for normal contour lines
        contour_width: Width of contour lines
        show_errors: If True and grid_result has 'errors', highlight geometry errors
        overlap_color: Color for overlap boundary lines (error code 1)
        undefined_color: Color for undefined boundary lines (error code 2)
        error_contour_style: Line style for error boundaries ('dashed' or 'solid')
        title, xlabel, ylabel: Labels
        overlay: 2D numpy array of scalar data to plot instead of geometry fill.
            When provided, this replaces the cell/material coloring. Geometry
            contours are still drawn on top. Useful for plotting MCNP mesh tally
            results (flux, dose, heating) on the geometry.
        overlay_extent: (h_min, h_max, v_min, v_max) extent of the overlay data.
            If None, uses the same extent as the geometry grid.
        overlay_cmap: Colormap for overlay data (default 'jet')
        overlay_norm: Normalization for overlay: 'log' for LogNorm, None for linear
        overlay_alpha: Opacity of overlay layer (0-1, default 1.0)
        overlay_vmin: Minimum value for overlay colormap scaling
        overlay_vmax: Maximum value for overlay colormap scaling
        overlay_label: Colorbar label for overlay data
        contour_by: Which ID to use for contour boundary detection:
            'cell' for cell boundaries, 'material' for material boundaries,
            or None to auto-select ('material' when by_material=True or
            overlay is set, otherwise 'cell').

    Returns:
        Matplotlib Axes
    """
    if not HAS_MATPLOTLIB:
        raise ImportError("matplotlib is required for plotting")
    if not HAS_NUMPY:
        raise ImportError("numpy is required for plotting")

    if ax is None:
        fig, ax = plt.subplots()

    # Determine grid dimensions and extent based on which axes are present
    if 'nu' in grid_result and 'nv' in grid_result:
        n_horiz = grid_result['nu']
        n_vert = grid_result['nv']
        h_min = grid_result['u_min']
        h_max = grid_result['u_max']
        v_min = grid_result['v_min']
        v_max = grid_result['v_max']
    elif 'nx' in grid_result and 'ny' in grid_result:
        n_horiz = grid_result['nx']
        n_vert = grid_result['ny']
        h_min = grid_result['x_min']
        h_max = grid_result['x_max']
        v_min = grid_result['y_min']
        v_max = grid_result['y_max']
    elif 'nx' in grid_result and 'nz' in grid_result:
        n_horiz = grid_result['nx']
        n_vert = grid_result['nz']
        h_min = grid_result['x_min']
        h_max = grid_result['x_max']
        v_min = grid_result['z_min']
        v_max = grid_result['z_max']
    elif 'ny' in grid_result and 'nz' in grid_result:
        n_horiz = grid_result['ny']
        n_vert = grid_result['nz']
        h_min = grid_result['y_min']
        h_max = grid_result['y_max']
        v_min = grid_result['z_min']
        v_max = grid_result['z_max']
    else:
        raise ValueError("Cannot determine grid dimensions from result")

    # Reshape ID arrays
    cell_ids = np.array(grid_result['cell_ids']).reshape(n_vert, n_horiz)
    material_ids = np.array(grid_result['material_ids']).reshape(n_vert, n_horiz)

    # Resolve contour_by default
    if contour_by is None:
        contour_by = 'material' if (by_material or overlay is not None) else 'cell'
    if contour_by not in ('cell', 'material'):
        raise ValueError(f"contour_by must be 'cell', 'material', or None, got {contour_by!r}")
    contour_ids = material_ids if contour_by == 'material' else cell_ids

    # Get color data
    if by_material:
        color_data = material_ids
    else:
        color_data = cell_ids.copy()

    extent = [h_min, h_max, v_min, v_max]

    # Set white background for void regions
    ax.set_facecolor('white')

    # Plot overlay data or geometry fill
    if overlay is not None:
        # Overlay mode: show scalar data (e.g. mesh tally) instead of geometry fill
        overlay_data = np.asarray(overlay, dtype=float)
        ov_extent = overlay_extent if overlay_extent is not None else extent

        # Build normalization
        norm = None
        if overlay_norm == 'log':
            from matplotlib.colors import LogNorm
            norm = LogNorm(vmin=overlay_vmin, vmax=overlay_vmax)
        elif overlay_vmin is not None or overlay_vmax is not None:
            from matplotlib.colors import Normalize
            norm = Normalize(vmin=overlay_vmin, vmax=overlay_vmax)

        ov_cmap = plt.colormaps.get_cmap(overlay_cmap).copy()
        im = ax.imshow(overlay_data, origin='lower', extent=ov_extent,
                       aspect='equal', cmap=ov_cmap, interpolation='nearest',
                       alpha=overlay_alpha, norm=norm)

        if show_colorbar:
            plt.colorbar(im, ax=ax, label=overlay_label or '')

    elif show_fill:
        # Normal geometry fill
        # Mask void regions: cell_id < 0 OR material_id == 0 (void material)
        if by_material:
            # When coloring by material, mask both undefined cells and void material
            masked_data = np.ma.masked_where(
                (cell_ids < 0) | (color_data == 0), color_data)
        else:
            # When coloring by cell, just mask undefined cells
            masked_data = np.ma.masked_where(color_data < 0, color_data)

        # Use colormap with white for masked (void) values
        cmap_obj = plt.colormaps.get_cmap(cmap).copy()
        cmap_obj.set_bad(color='white')

        im = ax.imshow(masked_data, origin='lower', extent=extent,
                       aspect='equal', cmap=cmap_obj, interpolation='nearest')

        if show_colorbar:
            plt.colorbar(im, ax=ax, label='Material ID' if by_material else 'Cell ID')

    # Get errors if present
    errors = None
    if 'errors' in grid_result:
        errors = np.array(grid_result['errors']).reshape(n_vert, n_horiz)

    # Draw contours by detecting boundaries between adjacent regions
    if show_contours:
        _draw_grid_contours(ax, contour_ids, errors, n_horiz, n_vert,
                           h_min, h_max, v_min, v_max,
                           contour_color, contour_width,
                           show_errors, overlap_color, undefined_color,
                           error_contour_style)

    ax.set_xlim(h_min, h_max)
    ax.set_ylim(v_min, v_max)

    if title:
        ax.set_title(title)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)

    return ax


def _draw_grid_contours(ax, ids: np.ndarray, errors: Optional[np.ndarray],
                        n_horiz: int, n_vert: int,
                        h_min: float, h_max: float, v_min: float, v_max: float,
                        contour_color: str, contour_width: float,
                        show_errors: bool, overlap_color: str, undefined_color: str,
                        error_contour_style: str) -> None:
    """Draw contours by detecting boundaries in grid (like plotter.c draw_contours_ex).

    This detects horizontal and vertical boundaries between adjacent pixels
    with different IDs (cell or material) and draws line segments.
    """
    # Pixel size in world coordinates
    dx = (h_max - h_min) / n_horiz
    dy = (v_max - v_min) / n_vert

    # Collect boundary segments by type for efficient drawing
    normal_h_segs = []  # Horizontal segments (between rows)
    normal_v_segs = []  # Vertical segments (between columns)
    overlap_h_segs = []
    overlap_v_segs = []
    undefined_h_segs = []
    undefined_v_segs = []

    # Detect horizontal boundaries (between row y and y+1)
    for y in range(n_vert - 1):
        for x in range(n_horiz):
            if ids[y, x] != ids[y + 1, x]:
                # Boundary between (x,y) and (x,y+1)
                x0 = h_min + x * dx
                x1 = h_min + (x + 1) * dx
                y_pos = v_min + (y + 1) * dy

                # Determine error type
                error_type = 0
                if show_errors and errors is not None:
                    e1, e2 = errors[y, x], errors[y + 1, x]
                    if e1 == GRID_ERROR_OVERLAP or e2 == GRID_ERROR_OVERLAP:
                        error_type = GRID_ERROR_OVERLAP
                    elif e1 == GRID_ERROR_UNDEFINED or e2 == GRID_ERROR_UNDEFINED:
                        error_type = GRID_ERROR_UNDEFINED

                if error_type == GRID_ERROR_OVERLAP:
                    overlap_h_segs.append(((x0, y_pos), (x1, y_pos)))
                elif error_type == GRID_ERROR_UNDEFINED:
                    undefined_h_segs.append(((x0, y_pos), (x1, y_pos)))
                else:
                    normal_h_segs.append(((x0, y_pos), (x1, y_pos)))

    # Detect vertical boundaries (between column x and x+1)
    for y in range(n_vert):
        for x in range(n_horiz - 1):
            if ids[y, x] != ids[y, x + 1]:
                # Boundary between (x,y) and (x+1,y)
                x_pos = h_min + (x + 1) * dx
                y0 = v_min + y * dy
                y1 = v_min + (y + 1) * dy

                # Determine error type
                error_type = 0
                if show_errors and errors is not None:
                    e1, e2 = errors[y, x], errors[y, x + 1]
                    if e1 == GRID_ERROR_OVERLAP or e2 == GRID_ERROR_OVERLAP:
                        error_type = GRID_ERROR_OVERLAP
                    elif e1 == GRID_ERROR_UNDEFINED or e2 == GRID_ERROR_UNDEFINED:
                        error_type = GRID_ERROR_UNDEFINED

                if error_type == GRID_ERROR_OVERLAP:
                    overlap_v_segs.append(((x_pos, y0), (x_pos, y1)))
                elif error_type == GRID_ERROR_UNDEFINED:
                    undefined_v_segs.append(((x_pos, y0), (x_pos, y1)))
                else:
                    normal_v_segs.append(((x_pos, y0), (x_pos, y1)))

    # Draw segments using LineCollection for efficiency
    from matplotlib.collections import LineCollection

    # Normal boundaries
    all_normal = normal_h_segs + normal_v_segs
    if all_normal:
        lc = LineCollection(all_normal, colors=contour_color, linewidths=contour_width)
        ax.add_collection(lc)

    # Error boundaries
    linestyle = '--' if error_contour_style == 'dashed' else '-'

    all_overlap = overlap_h_segs + overlap_v_segs
    if all_overlap:
        lc = LineCollection(all_overlap, colors=overlap_color, linewidths=contour_width,
                           linestyles=linestyle)
        ax.add_collection(lc)

    all_undefined = undefined_h_segs + undefined_v_segs
    if all_undefined:
        lc = LineCollection(all_undefined, colors=undefined_color, linewidths=contour_width,
                           linestyles=linestyle)
        ax.add_collection(lc)


# =============================================================================
# Error constants
# =============================================================================

# Special cell IDs for geometry errors (used in curve filtering)
CELL_VOID = -1
CELL_UNDEFINED = -2
CELL_OVERLAP = -3

# Grid error codes (from find_cells_grid_*_ex with detect_errors=True)
GRID_ERROR_OK = 0          # No error - valid cell
GRID_ERROR_OVERLAP = 1     # Multiple cells claim this point
GRID_ERROR_UNDEFINED = 2   # No cell claims this point

# =============================================================================
# Curve filtering - only show true cell boundaries
# =============================================================================


def _get_curve_sample_point(curve: Dict[str, Any]) -> Tuple[float, float]:
    """Get a representative point on the curve for sampling."""
    curve_type = curve.get('type', 'none')

    if curve_type == 'circle':
        # Point on right side of circle
        center = curve['center']
        radius = curve['radius']
        return (center[0] + radius, center[1])

    elif curve_type == 'arc':
        center = curve['center']
        radius = curve['radius']
        theta_start = curve.get('theta_start', 0)
        theta_end = curve.get('theta_end', 2 * np.pi)
        theta_mid = (theta_start + theta_end) / 2
        return (center[0] + radius * np.cos(theta_mid),
                center[1] + radius * np.sin(theta_mid))

    elif curve_type == 'ellipse':
        center = curve['center']
        semi_a = curve['semi_a']
        return (center[0] + semi_a, center[1])

    elif curve_type in ('line', 'line_segment'):
        point = curve['point']
        direction = curve['direction']
        t_mid = (curve.get('t_min', -100) + curve.get('t_max', 100)) / 2
        return (point[0] + t_mid * direction[0],
                point[1] + t_mid * direction[1])

    elif curve_type == 'polygon':
        vertices = curve['vertices']
        if len(vertices) >= 2:
            # Midpoint of first edge
            v0, v1 = vertices[0], vertices[1]
            return ((v0[0] + v1[0]) / 2, (v0[1] + v1[1]) / 2)
        elif vertices:
            return tuple(vertices[0])
        return (0.0, 0.0)

    elif curve_type == 'parallel_lines':
        # Midpoint between the two lines
        p1 = curve['point1']
        p2 = curve['point2']
        direction = curve['direction']
        t_mid = (curve.get('t_min', -100) + curve.get('t_max', 100)) / 2
        mid = ((p1[0] + p2[0]) / 2 + t_mid * direction[0],
               (p1[1] + p2[1]) / 2 + t_mid * direction[1])
        return mid

    return (0.0, 0.0)


def _get_curve_normal(curve: Dict[str, Any],
                      sample_pt: Tuple[float, float]) -> Tuple[float, float]:
    """Get outward normal to the curve at the sample point."""
    curve_type = curve.get('type', 'none')

    if curve_type in ('circle', 'arc'):
        center = curve['center']
        # Normal points outward from center
        dx = sample_pt[0] - center[0]
        dy = sample_pt[1] - center[1]
        length = np.sqrt(dx*dx + dy*dy)
        if length > 1e-10:
            return (dx / length, dy / length)
        return (1.0, 0.0)

    elif curve_type == 'ellipse':
        # For ellipse, normal is more complex, approximate
        center = curve['center']
        dx = sample_pt[0] - center[0]
        dy = sample_pt[1] - center[1]
        length = np.sqrt(dx*dx + dy*dy)
        if length > 1e-10:
            return (dx / length, dy / length)
        return (1.0, 0.0)

    elif curve_type in ('line', 'line_segment'):
        direction = curve['direction']
        # Normal is perpendicular to direction
        return (-direction[1], direction[0])

    elif curve_type == 'polygon':
        vertices = curve['vertices']
        if len(vertices) >= 2:
            v0, v1 = vertices[0], vertices[1]
            dx = v1[0] - v0[0]
            dy = v1[1] - v0[1]
            length = np.sqrt(dx*dx + dy*dy)
            if length > 1e-10:
                return (-dy / length, dx / length)
        return (1.0, 0.0)

    elif curve_type == 'parallel_lines':
        direction = curve['direction']
        return (-direction[1], direction[0])

    return (1.0, 0.0)


def filter_boundary_curves(curves_result: Dict[str, Any],
                           find_cell_func,
                           slice_coord: float,
                           slice_axis: str = 'z',
                           epsilon: float = 0.001) -> Dict[str, Any]:
    """Filter curves to only include true cell boundaries.

    A curve is a boundary if cells on either side are different.

    Args:
        curves_result: Result from get_slice_curves_*
        find_cell_func: Callable(x, y, z) -> cell_id (int or None)
        slice_coord: The coordinate value of the slice (z, y, or x value)
        slice_axis: Which axis is constant ('z', 'y', or 'x')
        epsilon: Small offset for sampling on each side of curve

    Returns:
        Filtered curves_result with only boundary curves.
        Each curve will have 'cell_left' and 'cell_right' attributes.
    """
    if not HAS_NUMPY:
        raise ImportError("numpy is required for curve filtering")

    filtered_curves = []

    for curve in curves_result.get('curves', []):
        # Get sample point on curve (in 2D slice coordinates)
        pt_2d = _get_curve_sample_point(curve)
        normal_2d = _get_curve_normal(curve, pt_2d)

        # Offset points on each side
        left_2d = (pt_2d[0] - normal_2d[0] * epsilon,
                   pt_2d[1] - normal_2d[1] * epsilon)
        right_2d = (pt_2d[0] + normal_2d[0] * epsilon,
                    pt_2d[1] + normal_2d[1] * epsilon)

        # Convert to 3D based on slice axis
        if slice_axis == 'z':
            # XY plane at z=slice_coord
            left_3d = (left_2d[0], left_2d[1], slice_coord)
            right_3d = (right_2d[0], right_2d[1], slice_coord)
        elif slice_axis == 'y':
            # XZ plane at y=slice_coord
            left_3d = (left_2d[0], slice_coord, left_2d[1])
            right_3d = (right_2d[0], slice_coord, right_2d[1])
        elif slice_axis == 'x':
            # YZ plane at x=slice_coord
            left_3d = (slice_coord, left_2d[0], left_2d[1])
            right_3d = (slice_coord, right_2d[0], right_2d[1])
        else:
            # For arbitrary planes, caller must handle 3D conversion
            # Default to treating as XY
            left_3d = (left_2d[0], left_2d[1], slice_coord)
            right_3d = (right_2d[0], right_2d[1], slice_coord)

        # Query cells on each side
        cell_left = find_cell_func(*left_3d)
        cell_right = find_cell_func(*right_3d)

        # Convert None to CELL_UNDEFINED
        if cell_left is None:
            cell_left = CELL_UNDEFINED
        if cell_right is None:
            cell_right = CELL_UNDEFINED

        # Only include if cells differ (it's a true boundary)
        if cell_left != cell_right:
            curve_copy = dict(curve)
            curve_copy['cell_left'] = cell_left
            curve_copy['cell_right'] = cell_right
            filtered_curves.append(curve_copy)

    # Build result with filtered curves
    result = dict(curves_result)
    result['curves'] = filtered_curves
    return result


def filter_boundary_curves_arbitrary(curves_result: Dict[str, Any],
                                     find_cell_func,
                                     origin: Tuple[float, float, float],
                                     u_axis: Tuple[float, float, float],
                                     v_axis: Tuple[float, float, float],
                                     epsilon: float = 0.001) -> Dict[str, Any]:
    """Filter curves for arbitrary plane slices.

    Args:
        curves_result: Result from get_slice_curves
        find_cell_func: Callable(x, y, z) -> cell_id
        origin: Origin point of the slice plane
        u_axis: Unit vector for U (horizontal) direction in 3D
        v_axis: Unit vector for V (vertical) direction in 3D
        epsilon: Sampling offset

    Returns:
        Filtered curves_result with only boundary curves.
    """
    if not HAS_NUMPY:
        raise ImportError("numpy is required for curve filtering")

    origin = np.array(origin)
    u_axis = np.array(u_axis)
    v_axis = np.array(v_axis)
    # Normal to plane
    normal = np.cross(u_axis, v_axis)
    normal = normal / np.linalg.norm(normal)

    filtered_curves = []

    for curve in curves_result.get('curves', []):
        pt_2d = _get_curve_sample_point(curve)
        normal_2d = _get_curve_normal(curve, pt_2d)

        # Offset in 2D
        left_2d = (pt_2d[0] - normal_2d[0] * epsilon,
                   pt_2d[1] - normal_2d[1] * epsilon)
        right_2d = (pt_2d[0] + normal_2d[0] * epsilon,
                    pt_2d[1] + normal_2d[1] * epsilon)

        # Convert to 3D: point = origin + u*u_axis + v*v_axis
        left_3d = origin + left_2d[0] * u_axis + left_2d[1] * v_axis
        right_3d = origin + right_2d[0] * u_axis + right_2d[1] * v_axis

        cell_left = find_cell_func(*left_3d)
        cell_right = find_cell_func(*right_3d)

        if cell_left is None:
            cell_left = CELL_UNDEFINED
        if cell_right is None:
            cell_right = CELL_UNDEFINED

        if cell_left != cell_right:
            curve_copy = dict(curve)
            curve_copy['cell_left'] = cell_left
            curve_copy['cell_right'] = cell_right
            filtered_curves.append(curve_copy)

    result = dict(curves_result)
    result['curves'] = filtered_curves
    return result


# =============================================================================
# Counting functions - extract cell/surface/material counts from slice results
# =============================================================================

def count_surfaces(curves_result: Dict[str, Any]) -> int:
    """Count unique surfaces in a slice from analytical curves.

    Args:
        curves_result: Result from get_slice_curves_*

    Returns:
        Number of unique surfaces in the slice
    """
    surface_ids: Set[int] = set()
    for curve in curves_result.get('curves', []):
        surface_id = curve.get('surface_id')
        if surface_id is not None:
            surface_ids.add(surface_id)
    return len(surface_ids)


def count_cells(grid_result: Dict[str, Any]) -> int:
    """Count unique cells in a slice from grid sampling.

    Args:
        grid_result: Result from find_cells_grid_*

    Returns:
        Number of unique cells in the slice (excluding void)
    """
    cell_ids = grid_result.get('cell_ids', [])
    unique = set(cell_ids)
    unique.discard(-1)  # Remove void
    return len(unique)


def count_materials(grid_result: Dict[str, Any]) -> int:
    """Count unique materials in a slice from grid sampling.

    Args:
        grid_result: Result from find_cells_grid_*

    Returns:
        Number of unique materials in the slice (excluding void)
    """
    material_ids = grid_result.get('material_ids', [])
    unique = set(material_ids)
    unique.discard(-1)  # Remove void (cell not found)
    unique.discard(0)   # Remove void (material 0)
    return len(unique)


def get_slice_stats(curves_result: Dict[str, Any],
                    grid_result: Dict[str, Any]) -> Dict[str, int]:
    """Get cell, surface, and material counts from slice results.

    Args:
        curves_result: Result from get_slice_curves_*
        grid_result: Result from find_cells_grid_*

    Returns:
        Dict with 'cells', 'surfaces', 'materials' counts
    """
    return {
        'cells': count_cells(grid_result),
        'surfaces': count_surfaces(curves_result),
        'materials': count_materials(grid_result),
    }


# =============================================================================
# High-level plot orchestration
# =============================================================================

def plot(model: 'Model',
         z: Optional[float] = None, y: Optional[float] = None,
         x: Optional[float] = None,
         origin: Optional[Tuple[float, float, float]] = None,
         normal: Optional[Tuple[float, float, float]] = None,
         up: Optional[Tuple[float, float, float]] = None,
         bounds: Optional[Tuple[float, float, float, float]] = None,
         filled: bool = True,
         by_material: bool = False,
         show_colorbar: bool = False,
         show_contours: bool = True,
         detect_errors: bool = False,
         universe_depth: int = -1,
         save: Optional[str] = None,
         ax=None,
         overlay=None,
         overlay_extent: Optional[Tuple[float, float, float, float]] = None,
         overlay_cmap: str = 'jet',
         overlay_norm: Optional[str] = None,
         overlay_alpha: float = 1.0,
         overlay_vmin: Optional[float] = None,
         overlay_vmax: Optional[float] = None,
         overlay_label: Optional[str] = None,
         contour_by: Optional[str] = None,
         **kwargs):
    """Plot geometry slice using grid-based rendering.

    Specify the slice plane using ONE of:
    - z: XY slice at z=constant
    - y: XZ slice at y=constant
    - x: YZ slice at x=constant
    - origin+normal: arbitrary plane

    If no plane is specified, defaults to XY slice at z=0.
    If bounds are not specified, auto-detects from geometry.

    Args:
        model: The geometry model.
        z, y, x: Coordinate for axis-aligned slices (use only one)
        origin: Point on arbitrary slice plane (x, y, z)
        normal: Normal vector for arbitrary plane
        up: Up vector for arbitrary plane orientation
        bounds: (min1, max1, min2, max2) slice bounds. Auto-detected if None.
        filled: If True, show filled regions; if False, contours only
        by_material: Color by material ID instead of cell ID
        show_colorbar: Show colorbar
        show_contours: If True, draw boundary contours between cells
        detect_errors: If True, detect geometry errors (overlaps, undefined)
        universe_depth: Universe hierarchy depth (-1 for innermost)
        save: Save to file instead of displaying
        ax: Matplotlib axes (creates new if None)
        overlay: 2D numpy array of scalar data to plot on the geometry.
        overlay_extent: (h_min, h_max, v_min, v_max) of the overlay data.
        overlay_cmap: Colormap for overlay (default 'jet')
        overlay_norm: 'log' for log scale, None for linear
        overlay_alpha: Overlay opacity (0-1, default 1.0)
        overlay_vmin: Min value for overlay color scaling
        overlay_vmax: Max value for overlay color scaling
        overlay_label: Colorbar label for overlay
        contour_by: Which ID to use for contour boundaries: 'cell',
            'material', or None (auto)
        **kwargs: Additional plotting options (resolution, contour_color, etc.)

    Returns:
        Matplotlib axes
    """
    if not HAS_MATPLOTLIB:
        raise ImportError("matplotlib is required for plotting. Install with: pip install matplotlib")

    model._rebuild_if_needed()

    # Default to z=0 if no plane specified
    if z is None and y is None and x is None and origin is None:
        z = 0.0

    # Auto-detect bounds if not specified
    if bounds is None:
        bounds = (-10.0, 10.0, -10.0, 10.0)

    # Create figure if needed
    own_fig = ax is None
    if own_fig:
        fig, ax = plt.subplots(figsize=(8, 8))
    else:
        fig = ax.get_figure()

    # Get grid resolution
    resolution = kwargs.pop('resolution', (100, 100))

    # Get grid data for the slice — use Model's axis-specific methods
    if z is not None:
        grid = model.find_cells_grid_z(z, bounds, resolution,
                                       universe_depth=universe_depth,
                                       detect_errors=detect_errors)
        xlabel, ylabel = 'X', 'Y'
        title = f'Z = {z}'
    elif y is not None:
        grid = model.find_cells_grid_y(y, bounds, resolution,
                                       universe_depth=universe_depth,
                                       detect_errors=detect_errors)
        xlabel, ylabel = 'X', 'Z'
        title = f'Y = {y}'
    elif x is not None:
        grid = model.find_cells_grid_x(x, bounds, resolution,
                                       universe_depth=universe_depth,
                                       detect_errors=detect_errors)
        xlabel, ylabel = 'Y', 'Z'
        title = f'X = {x}'
    elif origin is not None and normal is not None:
        if up is None:
            up = (0, 0, 1)
        grid = model.find_cells_grid(origin, normal, up, bounds, resolution,
                                     universe_depth=universe_depth,
                                     detect_errors=detect_errors)
        xlabel, ylabel = 'U', 'V'
        title = f'Plane at {origin}'
    else:
        raise ValueError("Must specify z, y, x, or origin+normal")

    # Plot using grid-based rendering
    ax = plot_slice_filled(
        grid,
        ax=ax,
        by_material=by_material,
        show_colorbar=show_colorbar,
        show_fill=filled,
        show_contours=show_contours,
        show_errors=detect_errors,
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
        overlay=overlay,
        overlay_extent=overlay_extent,
        overlay_cmap=overlay_cmap,
        overlay_norm=overlay_norm,
        overlay_alpha=overlay_alpha,
        overlay_vmin=overlay_vmin,
        overlay_vmax=overlay_vmax,
        overlay_label=overlay_label,
        contour_by=contour_by,
        **kwargs
    )

    if save:
        fig.savefig(save, dpi=150, bbox_inches='tight')
    elif own_fig:
        plt.show()

    return ax


def plot_views(model: 'Model',
               bounds: Optional[Tuple[float, float, float, float, float, float]] = None,
               save: Optional[str] = None,
               **kwargs):
    """Plot three orthogonal views (XY, XZ, YZ).

    Args:
        model: The geometry model.
        bounds: (xmin, xmax, ymin, ymax, zmin, zmax). Auto-detected if None.
        save: Save to file instead of displaying.
        **kwargs: Additional plotting options.

    Returns:
        Figure with three subplots.
    """
    if not HAS_MATPLOTLIB:
        raise ImportError("matplotlib is required for plotting")

    model._rebuild_if_needed()

    if bounds is None:
        bounds = (-10, 10, -10, 10, -10, 10)

    xmin, xmax, ymin, ymax, zmin, zmax = bounds

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # XY at z=0
    plot(model, z=0, bounds=(xmin, xmax, ymin, ymax), ax=axes[0], **kwargs)
    axes[0].set_title('XY (z=0)')

    # XZ at y=0
    plot(model, y=0, bounds=(xmin, xmax, zmin, zmax), ax=axes[1], **kwargs)
    axes[1].set_title('XZ (y=0)')

    # YZ at x=0
    plot(model, x=0, bounds=(ymin, ymax, zmin, zmax), ax=axes[2], **kwargs)
    axes[2].set_title('YZ (x=0)')

    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches='tight')
    else:
        plt.show()

    return fig
