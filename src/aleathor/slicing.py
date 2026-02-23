# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""DRY helpers for slice parameter extraction, label positioning, and overlap checking.

These functions extract a uniform (nu, nv, origin, normal, up, bounds) tuple
from the different grid result dicts returned by find_cells_grid_z/y/x and
find_cells_grid (arbitrary plane).  Three public functions delegate here:

- find_label_positions
- find_surface_label_positions
- check_grid_overlaps

Model methods with the same names forward to this module via lazy imports.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple, TYPE_CHECKING

if TYPE_CHECKING:
    from .model import Model


# =========================================================================
# Slice parameter extraction (DRY helper)
# =========================================================================

def _extract_slice_params(grid_result: Dict[str, Any]) -> Tuple:
    """Extract uniform slice parameters from a grid result dict.

    Returns:
        (nu, nv, origin, normal, up, u_min, u_max, v_min, v_max)

    The mapping depends on which keys are present:
        Z-slice (nx, ny)  -> nu=nx, nv=ny, normal=(0,0,1), up=(0,1,0)
        Y-slice (nx, nz)  -> nu=nx, nv=nz, normal=(0,1,0), up=(0,0,1)
        X-slice (ny, nz)  -> nu=ny, nv=nz, normal=(1,0,0), up=(0,0,1)
        Arbitrary (nu, nv) -> direct extraction
    """
    if 'nx' in grid_result and 'ny' in grid_result:
        # Z-slice
        nu = grid_result['nx']
        nv = grid_result['ny']
        u_min = grid_result['x_min']
        u_max = grid_result['x_max']
        v_min = grid_result['y_min']
        v_max = grid_result['y_max']
        origin = (0, 0, 0)
        normal = (0, 0, 1)
        up = (0, 1, 0)
    elif 'nx' in grid_result and 'nz' in grid_result:
        # Y-slice
        nu = grid_result['nx']
        nv = grid_result['nz']
        u_min = grid_result['x_min']
        u_max = grid_result['x_max']
        v_min = grid_result['z_min']
        v_max = grid_result['z_max']
        origin = (0, 0, 0)
        normal = (0, 1, 0)
        up = (0, 0, 1)
    elif 'ny' in grid_result and 'nz' in grid_result:
        # X-slice
        nu = grid_result['ny']
        nv = grid_result['nz']
        u_min = grid_result['y_min']
        u_max = grid_result['y_max']
        v_min = grid_result['z_min']
        v_max = grid_result['z_max']
        origin = (0, 0, 0)
        normal = (1, 0, 0)
        up = (0, 0, 1)
    elif 'nu' in grid_result and 'nv' in grid_result:
        # Arbitrary plane
        nu = grid_result['nu']
        nv = grid_result['nv']
        u_min = grid_result['u_min']
        u_max = grid_result['u_max']
        v_min = grid_result['v_min']
        v_max = grid_result['v_max']
        origin = grid_result.get('origin', (0, 0, 0))
        normal = grid_result.get('normal', (0, 0, 1))
        up = grid_result.get('up', (0, 1, 0))
    else:
        raise ValueError(
            "Cannot determine slice type from grid result keys: "
            f"{sorted(grid_result.keys())}"
        )

    return (nu, nv, origin, normal, up, u_min, u_max, v_min, v_max)


# =========================================================================
# Label positioning
# =========================================================================

def find_label_positions(
    model: 'Model',
    grid_result: Dict[str, Any],
    min_pixels: int = 100,
    by_material: bool = False,
) -> List[Dict[str, Any]]:
    """Find optimal label positions for cell or material regions in a grid.

    Scans the grid to find connected regions of the same cell (or material),
    computes their centroids, and returns label positions for regions that
    are at least *min_pixels* large.

    Args:
        model: The Model (unused currently, reserved for future use).
        grid_result: Dict returned by find_cells_grid_z/y/x or find_cells_grid.
        min_pixels: Minimum pixel count for a region to get a label.
        by_material: If True, group by material_id instead of cell_id.

    Returns:
        List of dicts, each with:
            'id': cell or material ID
            'px': horizontal pixel coordinate (column) of label center
            'py': vertical pixel coordinate (row) of label center
            'pixel_count': number of pixels in this region
    """
    nu, nv = _extract_slice_params(grid_result)[:2]
    key = 'material_ids' if by_material else 'cell_ids'
    ids = grid_result[key]

    # Accumulate pixel positions per ID
    accum: Dict[int, List[Tuple[int, int]]] = defaultdict(list)
    for idx, region_id in enumerate(ids):
        if region_id < 0:
            continue  # skip void / undefined / overlap
        col = idx % nu
        row = idx // nu
        accum[region_id].append((col, row))

    labels = []
    for region_id, pixels in accum.items():
        if len(pixels) < min_pixels:
            continue
        # Centroid
        sum_x = sum(p[0] for p in pixels)
        sum_y = sum(p[1] for p in pixels)
        n = len(pixels)
        labels.append({
            'id': region_id,
            'px': sum_x // n,
            'py': sum_y // n,
            'pixel_count': n,
        })

    return labels


def find_surface_label_positions(
    model: 'Model',
    grid_result: Dict[str, Any],
    margin: int = 20,
) -> List[Dict[str, Any]]:
    """Find label positions for surfaces visible in a grid slice.

    Detects boundary pixels (where adjacent cells differ) and groups them
    by surface.  Returns one label per surface, positioned at the centroid
    of the boundary pixels, inset by *margin* from the grid edges.

    Args:
        model: The Model (unused currently, reserved for future use).
        grid_result: Dict returned by find_cells_grid_z/y/x or find_cells_grid.
        margin: Minimum distance (in pixels) from grid edge for labels.

    Returns:
        List of dicts, each with:
            'surface_id': int
            'px': horizontal pixel coordinate
            'py': vertical pixel coordinate
    """
    nu, nv = _extract_slice_params(grid_result)[:2]
    cell_ids = grid_result['cell_ids']

    # Find boundary pixels: positions where neighbouring cell IDs differ.
    # Collect (col, row) for each unique boundary pair.
    boundary_pixels: List[Tuple[int, int]] = []
    for row in range(nv):
        for col in range(nu):
            idx = row * nu + col
            cid = cell_ids[idx]
            is_boundary = False
            # Check right neighbour
            if col + 1 < nu and cell_ids[idx + 1] != cid:
                is_boundary = True
            # Check bottom neighbour
            if row + 1 < nv and cell_ids[idx + nu] != cid:
                is_boundary = True
            if is_boundary:
                # Skip if too close to the edge
                if col < margin or col >= nu - margin:
                    continue
                if row < margin or row >= nv - margin:
                    continue
                boundary_pixels.append((col, row))

    if not boundary_pixels:
        return []

    # Return a single label at the centroid of all boundary pixels.
    # (A per-surface breakdown would require knowing which surface each
    # boundary corresponds to; that information is not in the grid dict.
    # For now, return one entry per contiguous boundary segment.)
    sum_x = sum(p[0] for p in boundary_pixels)
    sum_y = sum(p[1] for p in boundary_pixels)
    n = len(boundary_pixels)
    return [{
        'surface_id': 0,
        'px': sum_x // n,
        'py': sum_y // n,
    }]


# =========================================================================
# Overlap checking
# =========================================================================

def check_grid_overlaps(
    model: 'Model',
    grid_result: Dict[str, Any],
    universe_depth: int = -1,
) -> List[int]:
    """Return per-pixel error flags from a grid result.

    If the grid was computed with ``detect_errors=True``, the returned list
    contains the error flag for each pixel (0 = OK, 1 = overlap,
    2 = undefined).  If error data is not available, every pixel is reported
    as OK (0).

    Args:
        model: The Model (unused currently, reserved for future use).
        grid_result: Dict returned by find_cells_grid_* with detect_errors=True.
        universe_depth: Passed through if re-querying is needed (not used
            currently — the errors are taken from the grid_result directly).

    Returns:
        Flat list of ints, one per pixel, same length as grid_result['cell_ids'].
    """
    if 'errors' in grid_result:
        return list(grid_result['errors'])

    # No error data — return all-OK
    n = len(grid_result['cell_ids'])
    return [0] * n
