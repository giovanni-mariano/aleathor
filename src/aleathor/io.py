# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""
Input/Output functions for reading and writing geometry files.
"""

from __future__ import annotations
from pathlib import Path
from typing import Union, Optional

from .model import Model, _CellData as Cell, Material
from .geometry import Region, Halfspace, Intersection, Union as UnionRegion, Complement
from .surfaces import (
    Surface, Plane, XPlane, YPlane, ZPlane, Sphere,
    CylinderX, CylinderY, CylinderZ, Box, RCC
)

try:
    import _alea
except ImportError:
    _alea = None


def read_mcnp(filename: Union[str, Path]) -> Model:
    """Read MCNP input file and create a Model.

    Args:
        filename: Path to MCNP input file.

    Returns:
        Model with loaded geometry.

    Raises:
        IOError: If file cannot be read.
        ValueError: If file contains parse errors.
    """
    if _alea is None:
        raise RuntimeError("C extension not available")

    path = Path(filename)
    if not path.exists():
        raise IOError(f"File not found: {filename}")

    # Load via C library
    sys = _alea.load_mcnp(str(path))

    # Build model from loaded system
    model = Model(title=path.stem)
    model._sys = sys
    model._dirty = False

    # Extract cells from C system
    sys.build_universe_index()
    cells = sys.get_cells()

    for cell_info in cells:
        # Create imported region with reference to C system for point queries
        cell = Cell(
            id=cell_info['cell_id'],
            region=_ImportedRegion(sys, cell_info['root_node']),
            material=cell_info['material_id'],
            density=cell_info['density'],
            density_unit="g/cm3" if cell_info.get('is_mass_density', True) else "atoms/b-cm",
            universe=cell_info['universe_id'],
            fill=cell_info['fill_universe'] if cell_info['fill_universe'] >= 0 else None,
        )
        model._cells[cell.id] = cell

    return model


def read_openmc(filename: Union[str, Path]) -> Model:
    """Read OpenMC XML geometry file and create a Model.

    Args:
        filename: Path to OpenMC XML geometry file.

    Returns:
        Model with loaded geometry.

    Raises:
        IOError: If file cannot be read.
        ValueError: If file contains parse errors.
    """
    if _alea is None:
        raise RuntimeError("C extension not available")

    path = Path(filename)
    if not path.exists():
        raise IOError(f"File not found: {filename}")

    # Load via C library
    sys = _alea.load_openmc(str(path))

    # Build model from loaded system
    model = Model(title=path.stem)
    model._sys = sys
    model._dirty = False

    # Extract cells from C system
    sys.build_universe_index()
    cells = sys.get_cells()

    for cell_info in cells:
        # Create imported region with reference to C system for point queries
        cell = Cell(
            id=cell_info['cell_id'],
            region=_ImportedRegion(sys, cell_info['root_node']),
            material=cell_info['material_id'],
            density=cell_info['density'],
            density_unit="g/cm3" if cell_info.get('is_mass_density', True) else "atoms/b-cm",
            universe=cell_info['universe_id'],
            fill=cell_info['fill_universe'] if cell_info['fill_universe'] >= 0 else None,
        )
        model._cells[cell.id] = cell

    return model


def read_mcnp_string(content: str) -> Model:
    """Read MCNP input from string and create a Model.

    Args:
        content: MCNP input file contents.

    Returns:
        Model with loaded geometry.
    """
    if _alea is None:
        raise RuntimeError("C extension not available")

    sys = _alea.load_mcnp_string(content)

    model = Model(title="MCNP Model")
    model._sys = sys
    model._dirty = False

    sys.build_universe_index()
    cells = sys.get_cells()

    for cell_info in cells:
        cell = Cell(
            id=cell_info['cell_id'],
            region=_ImportedRegion(sys, cell_info['root_node']),
            material=cell_info['material_id'],
            density=cell_info['density'],
            density_unit="g/cm3" if cell_info.get('is_mass_density', True) else "atoms/b-cm",
            universe=cell_info['universe_id'],
            fill=cell_info['fill_universe'] if cell_info['fill_universe'] >= 0 else None,
        )
        model._cells[cell.id] = cell

    return model


def write_mcnp(model: Model, filename: Union[str, Path],
               deduplicate: bool = True) -> None:
    """Write model to MCNP format.

    Args:
        model: Model to export.
        filename: Output file path.
        deduplicate: Deduplicate surfaces.
    """
    model.export_mcnp(filename, deduplicate=deduplicate)


def write_openmc(model: Model, filename: Union[str, Path]) -> None:
    """Write model to OpenMC XML format.

    Args:
        model: Model to export.
        filename: Output file path.
    """
    model.export_openmc(filename)


class _SurfaceRef:
    """Lightweight reference to a surface by MCNP ID (for imported models)."""
    __slots__ = ('id',)

    def __init__(self, surface_id: int):
        self.id = surface_id

    def __repr__(self):
        return str(self.id)

    def __eq__(self, other):
        return isinstance(other, _SurfaceRef) and self.id == other.id

    def __hash__(self):
        return hash(self.id)


_OP_PRIMITIVE = 0
_OP_UNION = 1
_OP_INTERSECTION = 2
_OP_DIFFERENCE = 3
_OP_COMPLEMENT = 4


def _unpack_csg_tree(tree):
    """Convert nested tuple from C node_tree() into a Python Region tree."""
    op = tree[0]

    if op == _OP_PRIMITIVE:
        _, surface_id, sense = tree
        return Halfspace(_SurfaceRef(surface_id), positive=(sense > 0))

    if op == _OP_COMPLEMENT:
        return Complement(_unpack_csg_tree(tree[1]))

    left = _unpack_csg_tree(tree[1])
    right = _unpack_csg_tree(tree[2])

    if op == _OP_INTERSECTION:
        return Intersection(left, right)
    if op == _OP_UNION:
        return UnionRegion(left, right)
    if op == _OP_DIFFERENCE:
        return Intersection(left, Complement(right))

    raise ValueError(f"Unknown CSG operation: {op}")


class _ImportedRegion(Region):
    """Region imported from MCNP/C library.

    Holds a reference to the C system and node for fast point containment
    queries. The CSG tree can be reconstructed lazily via the ``tree``
    property, which walks the C tree once and caches the result.
    """

    def __init__(self, sys, node_id: int):
        self._sys = sys
        self._node_id = node_id
        self._cached_tree = None

    @property
    def tree(self):
        """Lazily reconstruct the Python Region tree from the C CSG tree."""
        if self._cached_tree is None:
            raw = self._sys.node_tree(self._node_id)
            self._cached_tree = _unpack_csg_tree(raw)
        return self._cached_tree

    def __contains__(self, point):
        """Point containment via C library."""
        x, y, z = point
        return self._sys.point_inside(self._node_id, x, y, z)

    def get_surfaces(self):
        """Get all surfaces referenced in this region's CSG tree."""
        return self.tree.get_surfaces()

    def _to_csg(self, model):
        return self._node_id

    def __repr__(self):
        return repr(self.tree)


# Backward compatibility alias
_PlaceholderRegion = _ImportedRegion
