#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""
Command-line tool to plot MCNP/OpenMC geometry slices.

Equivalent to the C plotter example, with support for:
- Axis-aligned slices (X, Y, Z planes)
- Arbitrary planes defined by origin + normal + up
- Planes defined by 3 points
- Cell/material labels
- Color by cell ID or material ID
- Axis tick marks
- Error detection (overlaps, undefined regions)
- Batch mode for multiple plots
- PNG/PDF output

Usage (single plot):
    python plot_geometry.py input.inp Z 0 -50 50 -50 50 800 output.png
    python plot_geometry.py input.inp Y 0 -50 50 -50 50 800x600 output.png --labels=cells
    python plot_geometry.py input.inp X 0 -50 50 -50 50 800 output.png --color=materials

Usage (batch mode):
    python plot_geometry.py input.inp --batch=plots.txt

Batch file format (one plot per line):
    Z value u_min u_max v_min v_max WxH output.png [options]
    Y value u_min u_max v_min v_max WxH output.png [options]
    X value u_min u_max v_min v_max WxH output.png [options]
    PLANE ox oy oz nx ny nz ux uy uz u_min u_max v_min v_max WxH output [options]
    PLANE3 x1 y1 z1 x2 y2 z2 x3 y3 z3 u_min u_max v_min v_max WxH output [options]

Options:
    --labels=cells          Show cell ID labels
    --labels=materials      Show material ID labels
    --labels=all            Show all labels
    --color=cells           Color by cell ID
    --color=materials       Color by material ID (default)
    --cmap=<name>           Matplotlib colormap (default: tab20)
    --no-ticks              Hide axis tick marks (shown by default)
    --errors                Detect and highlight geometry errors
    --no-contours           Hide boundary contours
"""

import argparse
import sys
import time
from pathlib import Path
from dataclasses import dataclass
from typing import Optional, Tuple, List
import numpy as np


@dataclass
class PlotParams:
    """Parameters for a single plot."""
    slice_type: str  # 'Z', 'Y', 'X', or 'PLANE'
    value: float = 0.0  # Axis value for Z/Y/X slices
    u_min: float = -10.0
    u_max: float = 10.0
    v_min: float = -10.0
    v_max: float = 10.0
    width: int = 800
    height: int = 800
    output: str = "plot.png"

    # For arbitrary planes
    origin: Optional[Tuple[float, float, float]] = None
    normal: Optional[Tuple[float, float, float]] = None
    up: Optional[Tuple[float, float, float]] = None

    # Options
    show_cell_labels: bool = False
    show_material_labels: bool = False
    show_surface_labels: bool = False
    color_by_material: bool = True  # Default: color by material
    show_ticks: bool = True  # Default: show ticks
    detect_errors: bool = False
    show_contours: bool = True
    cmap: str = 'tab20'  # Colormap name


def parse_resolution(res_str: str) -> Tuple[int, int]:
    """Parse resolution string like '800' or '800x600'."""
    if 'x' in res_str.lower():
        parts = res_str.lower().split('x')
        return int(parts[0]), int(parts[1])
    else:
        size = int(res_str)
        return size, size


def parse_labels(label_str: str) -> Tuple[bool, bool, bool]:
    """Parse labels option string."""
    cells, materials, surfaces = False, False, False
    for part in label_str.lower().split(','):
        part = part.strip()
        if part == 'cells':
            cells = True
        elif part == 'materials':
            materials = True
        elif part == 'surfaces':
            surfaces = True
        elif part == 'all':
            cells = materials = surfaces = True
    return cells, materials, surfaces


def plane_from_3_points(p1: Tuple[float, float, float],
                        p2: Tuple[float, float, float],
                        p3: Tuple[float, float, float]) -> Tuple[Tuple, Tuple, Tuple]:
    """Compute plane parameters (origin, normal, up) from 3 points."""
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)

    # Origin = P1
    origin = tuple(p1)

    # U = P2 - P1 (horizontal direction)
    u = p2 - p1
    u_len = np.linalg.norm(u)
    if u_len > 1e-10:
        u = u / u_len

    # V' = P3 - P1
    v_prime = p3 - p1

    # Normal = U x V'
    normal = np.cross(u, v_prime)
    normal_len = np.linalg.norm(normal)
    if normal_len > 1e-10:
        normal = normal / normal_len

    # Up = use V' direction (will be orthogonalized)
    up = tuple(v_prime)

    return origin, tuple(normal), up


def parse_batch_line(line: str) -> Optional[PlotParams]:
    """Parse a single line from batch file."""
    line = line.strip()

    # Skip comments and empty lines
    if not line or line.startswith('#'):
        return None

    parts = line.split()
    if len(parts) < 2:
        return None

    params = PlotParams(slice_type='Z')

    try:
        slice_type = parts[0].upper()

        if slice_type in ('Z', 'Y', 'X'):
            # Axis-aligned slice: TYPE value u_min u_max v_min v_max WxH output [options]
            if len(parts) < 8:
                return None
            params.slice_type = slice_type
            params.value = float(parts[1])
            params.u_min = float(parts[2])
            params.u_max = float(parts[3])
            params.v_min = float(parts[4])
            params.v_max = float(parts[5])
            params.width, params.height = parse_resolution(parts[6])
            params.output = parts[7]
            option_start = 8

        elif slice_type == 'PLANE':
            # PLANE ox oy oz nx ny nz ux uy uz u_min u_max v_min v_max WxH output [options]
            if len(parts) < 16:
                return None
            params.slice_type = 'PLANE'
            params.origin = (float(parts[1]), float(parts[2]), float(parts[3]))
            params.normal = (float(parts[4]), float(parts[5]), float(parts[6]))
            params.up = (float(parts[7]), float(parts[8]), float(parts[9]))
            params.u_min = float(parts[10])
            params.u_max = float(parts[11])
            params.v_min = float(parts[12])
            params.v_max = float(parts[13])
            params.width, params.height = parse_resolution(parts[14])
            params.output = parts[15]
            option_start = 16

        elif slice_type == 'PLANE3':
            # PLANE3 x1 y1 z1 x2 y2 z2 x3 y3 z3 u_min u_max v_min v_max WxH output [options]
            if len(parts) < 16:
                return None
            params.slice_type = 'PLANE'
            p1 = (float(parts[1]), float(parts[2]), float(parts[3]))
            p2 = (float(parts[4]), float(parts[5]), float(parts[6]))
            p3 = (float(parts[7]), float(parts[8]), float(parts[9]))
            params.origin, params.normal, params.up = plane_from_3_points(p1, p2, p3)
            params.u_min = float(parts[10])
            params.u_max = float(parts[11])
            params.v_min = float(parts[12])
            params.v_max = float(parts[13])
            params.width, params.height = parse_resolution(parts[14])
            params.output = parts[15]
            option_start = 16
        else:
            return None

        # Parse options
        for opt in parts[option_start:]:
            if opt.startswith('labels='):
                cells, mats, surfs = parse_labels(opt[7:])
                params.show_cell_labels = cells
                params.show_material_labels = mats
                params.show_surface_labels = surfs
            elif opt.startswith('color='):
                params.color_by_material = opt[6:].lower() not in ('cells', 'cell')
            elif opt.startswith('cmap='):
                params.cmap = opt[5:]
            elif opt == 'no-ticks':
                params.show_ticks = False
            elif opt == 'errors':
                params.detect_errors = True
            elif opt == 'no-contours':
                params.show_contours = False

    except (ValueError, IndexError):
        return None

    return params


def render_plot(model, params: PlotParams, verbose: bool = False) -> bool:
    """Render a single plot."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgba

    try:
        import aleathor as ath
    except ImportError:
        print("Error: aleathor package not found", file=sys.stderr)
        return False

    t0 = time.time()

    if verbose:
        if params.slice_type == 'PLANE':
            print(f"  Plane origin={params.origin} normal={params.normal}, "
                  f"{params.width}x{params.height} -> {params.output}")
        else:
            print(f"  {params.slice_type}={params.value:.4g}, "
                  f"bounds=[{params.u_min:.4g},{params.u_max:.4g}]x"
                  f"[{params.v_min:.4g},{params.v_max:.4g}], "
                  f"{params.width}x{params.height} -> {params.output}")

    bounds = (params.u_min, params.u_max, params.v_min, params.v_max)
    resolution = (params.width, params.height)

    # Get grid data (this is the only query needed - contours are derived from grid)
    if params.slice_type == 'Z':
        grid = model.find_cells_grid_z(params.value, bounds, resolution,
                                       detect_errors=params.detect_errors)
        xlabel, ylabel = 'X', 'Y'
        title = f'Z = {params.value:.4g}'
    elif params.slice_type == 'Y':
        grid = model.find_cells_grid_y(params.value, bounds, resolution,
                                       detect_errors=params.detect_errors)
        xlabel, ylabel = 'X', 'Z'
        title = f'Y = {params.value:.4g}'
    elif params.slice_type == 'X':
        grid = model.find_cells_grid_x(params.value, bounds, resolution,
                                       detect_errors=params.detect_errors)
        xlabel, ylabel = 'Y', 'Z'
        title = f'X = {params.value:.4g}'
    elif params.slice_type == 'PLANE':
        grid = model.find_cells_grid(params.origin, params.normal, params.up,
                                     bounds, resolution, detect_errors=params.detect_errors)
        xlabel, ylabel = 'U', 'V'
        title = f'Plane at {params.origin}'
    else:
        print(f"Error: Unknown slice type {params.slice_type}", file=sys.stderr)
        return False

    t1 = time.time()
    if verbose:
        n_pixels = params.width * params.height
        print(f"    Grid query: {(t1-t0)*1000:.1f} ms ({n_pixels/(t1-t0)/1e6:.2f} Mpx/s)")

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 10 * params.height / params.width))

    # Plot filled regions with grid-based contours (like plotter.c)
    ath.plot_slice_filled(
        grid,
        ax=ax,
        cmap=params.cmap,
        by_material=params.color_by_material,
        show_contours=params.show_contours,
        show_errors=params.detect_errors,
        contour_color='black',
        contour_width=0.5
    )

    # Add labels if requested
    if params.show_cell_labels or params.show_material_labels:
        labels = model.find_label_positions(grid, min_pixels=100,
                                            by_material=params.show_material_labels)

        # Determine grid dimensions
        if 'nx' in grid and 'ny' in grid:
            n_horiz, n_vert = grid['nx'], grid['ny']
            h_min, h_max = grid['x_min'], grid['x_max']
            v_min, v_max = grid['y_min'], grid['y_max']
        elif 'nx' in grid and 'nz' in grid:
            n_horiz, n_vert = grid['nx'], grid['nz']
            h_min, h_max = grid['x_min'], grid['x_max']
            v_min, v_max = grid['z_min'], grid['z_max']
        elif 'ny' in grid and 'nz' in grid:
            n_horiz, n_vert = grid['ny'], grid['nz']
            h_min, h_max = grid['y_min'], grid['y_max']
            v_min, v_max = grid['z_min'], grid['z_max']
        else:
            n_horiz, n_vert = grid['nu'], grid['nv']
            h_min, h_max = grid['u_min'], grid['u_max']
            v_min, v_max = grid['v_min'], grid['v_max']

        for label in labels:
            # Convert pixel to world coordinates
            x = h_min + (label['px'] + 0.5) / n_horiz * (h_max - h_min)
            y = v_min + (label['py'] + 0.5) / n_vert * (v_max - v_min)

            prefix = 'M' if params.show_material_labels else 'C'
            text = f"{prefix}{label['id']}"

            ax.text(x, y, text, ha='center', va='center',
                   fontsize=8, color='white', fontweight='bold',
                   bbox=dict(boxstyle='round,pad=0.2', facecolor='black', alpha=0.6))

    # Add tick marks if requested
    if params.show_ticks:
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.tick_params(axis='both', which='major', labelsize=8)
    else:
        ax.set_xticks([])
        ax.set_yticks([])

    ax.set_title(title)
    ax.set_aspect('equal')

    # Save
    output_path = Path(params.output)
    fig.savefig(output_path, dpi=150, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close(fig)

    t2 = time.time()
    if verbose:
        print(f"    Total: {(t2-t0)*1000:.1f} ms")

    return True


def run_batch(batch_file: str, model, verbose: bool = False) -> int:
    """Run batch processing from file."""
    try:
        with open(batch_file, 'r') as f:
            lines = f.readlines()
    except IOError as e:
        print(f"Error reading batch file: {e}", file=sys.stderr)
        return 1

    print(f"Processing batch file: {batch_file}")
    print()

    plot_count = 0
    error_count = 0

    for line_num, line in enumerate(lines, 1):
        params = parse_batch_line(line)
        if params is None:
            continue

        plot_count += 1
        print(f"Plot {plot_count}:")

        if not render_plot(model, params, verbose=True):
            error_count += 1

    print()
    print(f"Batch complete: {plot_count} plots, {error_count} errors")

    return 1 if error_count > 0 else 0


def main():
    parser = argparse.ArgumentParser(
        description="Plot MCNP/OpenMC geometry slices (equivalent to plotter.c)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Single plot examples:
  %(prog)s input.inp Z 0 -50 50 -50 50 800 output.png
  %(prog)s input.inp Y 0 -50 50 -50 50 800x600 slice_y.png --labels=cells
  %(prog)s input.inp X 0 -50 50 -50 50 1000 slice_x.png --color=materials --ticks

Arbitrary plane (origin + normal + up):
  %(prog)s input.inp PLANE 0 0 0  1 1 0  0 0 1  -50 50 -50 50 800 diagonal.png

Plane from 3 points:
  %(prog)s input.inp PLANE3 0 0 0  10 0 0  0 10 0  -50 50 -50 50 800 plane.png

Batch mode:
  %(prog)s input.inp --batch=plots.txt

Batch file format (one plot per line):
  Z value u_min u_max v_min v_max WxH output.png [options]
  Y value u_min u_max v_min v_max WxH output.png [options]
  X value u_min u_max v_min v_max WxH output.png [options]
  PLANE ox oy oz nx ny nz ux uy uz u_min u_max v_min v_max WxH output [options]
  PLANE3 x1 y1 z1 x2 y2 z2 x3 y3 z3 u_min u_max v_min v_max WxH output [options]

Options (for command line use --option, in batch file use option without --):
  --labels=cells|materials|all    Show labels
  --color=cells|materials         Color mode (materials by default)
  --cmap=<name>                   Colormap (default: tab20)
  --no-ticks                      Hide axis ticks (shown by default)
  --errors                        Detect geometry errors
  --no-contours                   Hide boundary lines
"""
    )

    parser.add_argument("input", help="Input geometry file (MCNP or OpenMC)")
    parser.add_argument("args", nargs='*', help="Slice parameters (see examples)")

    parser.add_argument("--batch", metavar="FILE",
                        help="Run batch mode from file")
    parser.add_argument("--labels", metavar="TYPE",
                        help="Show labels: cells, materials, surfaces, or all")
    parser.add_argument("--color", metavar="MODE", default="materials",
                        help="Color mode: cells or materials (default)")
    parser.add_argument("--cmap", metavar="NAME", default="tab20",
                        help="Matplotlib colormap name (default: tab20)")
    parser.add_argument("--no-ticks", action="store_true",
                        help="Hide axis tick marks (shown by default)")
    parser.add_argument("--errors", action="store_true",
                        help="Detect and highlight geometry errors")
    parser.add_argument("--no-contours", action="store_true",
                        help="Hide boundary contours")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose output")

    args = parser.parse_args()

    # Import aleathor
    try:
        import aleathor as ath
    except ImportError:
        print("Error: aleathor package not found.", file=sys.stderr)
        print("Install with: pip install aleathor", file=sys.stderr)
        return 1

    # Check matplotlib
    try:
        import matplotlib
    except ImportError:
        print("Error: matplotlib is required for plotting.", file=sys.stderr)
        print("Install with: pip install matplotlib", file=sys.stderr)
        return 1

    # Check input file
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: File not found: {args.input}", file=sys.stderr)
        return 1

    # Load geometry
    print(f"Loading geometry: {args.input}")
    t0 = time.time()

    try:
        model = ath.load(args.input)
    except Exception as e:
        print(f"Error loading geometry: {e}", file=sys.stderr)
        return 1

    t1 = time.time()
    print(f"  Loaded in {(t1-t0)*1000:.1f} ms")
    print(f"  Cells: {len(model.cells)}")

    # Build spatial index for faster queries
    print("Building spatial index...")
    t0 = time.time()
    model.build_spatial_index()
    t1 = time.time()
    print(f"  Built in {(t1-t0)*1000:.1f} ms")
    print()

    # Batch mode
    if args.batch:
        return run_batch(args.batch, model, args.verbose)

    # Single plot mode - parse positional arguments
    pos_args = args.args
    if len(pos_args) < 8:
        print("Error: Not enough arguments for single plot mode.", file=sys.stderr)
        print("Usage: plot_geometry.py input.inp <axis> <value> <u_min> <u_max> <v_min> <v_max> <WxH> [output.png] [options]", file=sys.stderr)
        print("       plot_geometry.py input.inp --batch=<file.txt>", file=sys.stderr)
        return 1

    params = PlotParams(slice_type='Z')

    try:
        slice_type = pos_args[0].upper()

        if slice_type in ('Z', 'Y', 'X'):
            params.slice_type = slice_type
            params.value = float(pos_args[1])
            params.u_min = float(pos_args[2])
            params.u_max = float(pos_args[3])
            params.v_min = float(pos_args[4])
            params.v_max = float(pos_args[5])
            params.width, params.height = parse_resolution(pos_args[6])
            params.output = pos_args[7] if len(pos_args) > 7 else "plot.png"

        elif slice_type == 'PLANE':
            if len(pos_args) < 15:
                print("Error: PLANE requires: ox oy oz nx ny nz ux uy uz u_min u_max v_min v_max WxH [output]", file=sys.stderr)
                return 1
            params.slice_type = 'PLANE'
            params.origin = (float(pos_args[1]), float(pos_args[2]), float(pos_args[3]))
            params.normal = (float(pos_args[4]), float(pos_args[5]), float(pos_args[6]))
            params.up = (float(pos_args[7]), float(pos_args[8]), float(pos_args[9]))
            params.u_min = float(pos_args[10])
            params.u_max = float(pos_args[11])
            params.v_min = float(pos_args[12])
            params.v_max = float(pos_args[13])
            params.width, params.height = parse_resolution(pos_args[14])
            params.output = pos_args[15] if len(pos_args) > 15 else "plot.png"

        elif slice_type == 'PLANE3':
            if len(pos_args) < 15:
                print("Error: PLANE3 requires: x1 y1 z1 x2 y2 z2 x3 y3 z3 u_min u_max v_min v_max WxH [output]", file=sys.stderr)
                return 1
            params.slice_type = 'PLANE'
            p1 = (float(pos_args[1]), float(pos_args[2]), float(pos_args[3]))
            p2 = (float(pos_args[4]), float(pos_args[5]), float(pos_args[6]))
            p3 = (float(pos_args[7]), float(pos_args[8]), float(pos_args[9]))
            params.origin, params.normal, params.up = plane_from_3_points(p1, p2, p3)
            params.u_min = float(pos_args[10])
            params.u_max = float(pos_args[11])
            params.v_min = float(pos_args[12])
            params.v_max = float(pos_args[13])
            params.width, params.height = parse_resolution(pos_args[14])
            params.output = pos_args[15] if len(pos_args) > 15 else "plot.png"
        else:
            print(f"Error: Unknown slice type '{slice_type}'. Use Z, Y, X, PLANE, or PLANE3.", file=sys.stderr)
            return 1

    except (ValueError, IndexError) as e:
        print(f"Error parsing arguments: {e}", file=sys.stderr)
        return 1

    # Apply command-line options
    if args.labels:
        cells, mats, surfs = parse_labels(args.labels)
        params.show_cell_labels = cells
        params.show_material_labels = mats
        params.show_surface_labels = surfs

    params.color_by_material = args.color.lower() not in ('cells', 'cell')
    params.cmap = args.cmap
    params.show_ticks = not args.no_ticks
    params.detect_errors = args.errors
    params.show_contours = not args.no_contours

    # Render
    print("Rendering...")
    if render_plot(model, params, verbose=True):
        print(f"\nDone: {params.output}")
        return 0
    else:
        return 1


if __name__ == "__main__":
    sys.exit(main())
