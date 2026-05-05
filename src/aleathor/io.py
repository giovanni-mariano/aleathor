# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""
Input/Output functions for reading and writing geometry files.
"""

from __future__ import annotations
from pathlib import Path
from typing import Union, Optional

from .model import Model
from .geometry import Region, Halfspace, Intersection, Union as UnionRegion, Complement
from .surfaces import (
    Surface, Plane, XPlane, YPlane, ZPlane, Sphere,
    XCylinder, YCylinder, ZCylinder, XCone, YCone, ZCone,
    XTorus, YTorus, ZTorus, RPP, Quadric,
    RCC, Box, TRC, ELL, REC, WED, RHP,
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
    _populate_regions(model, sys)

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
    _populate_regions(model, sys)

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
    _populate_regions(model, sys)

    return model


def read_openmc_string(content: str) -> Model:
    """Read OpenMC XML input from string and create a Model.

    Args:
        content: OpenMC XML file contents.

    Returns:
        Model with loaded geometry.
    """
    if _alea is None:
        raise RuntimeError("C extension not available")

    sys = _alea.load_openmc_string(content)

    model = Model(title="OpenMC Model")
    model._sys = sys
    _populate_regions(model, sys)

    return model


# Map ALEA_PRIMITIVE_* enum values (alea_types.h) to constructors that
# rebuild a Python Surface from the dict returned by node_primitive_data.
# Note: libalea's Plane stores ax+by+cz+d=0; aleathor's Plane uses
# ax+by+cz=d, so d is negated on the round-trip.
_PRIMITIVE_BUILDERS = {
    1:  lambda d, sid: Plane(d['a'], d['b'], d['c'], -d['d'], surface_id=sid),
    2:  lambda d, sid: Sphere(d['center_x'], d['center_y'], d['center_z'], d['radius'], surface_id=sid),
    3:  lambda d, sid: XCylinder(d['center_y'], d['center_z'], d['radius'], surface_id=sid),
    4:  lambda d, sid: YCylinder(d['center_x'], d['center_z'], d['radius'], surface_id=sid),
    5:  lambda d, sid: ZCylinder(d['center_x'], d['center_y'], d['radius'], surface_id=sid),
    6:  lambda d, sid: XCone(d['apex_x'], d['apex_y'], d['apex_z'], d['tan_angle_sq'], surface_id=sid),
    7:  lambda d, sid: YCone(d['apex_x'], d['apex_y'], d['apex_z'], d['tan_angle_sq'], surface_id=sid),
    8:  lambda d, sid: ZCone(d['apex_x'], d['apex_y'], d['apex_z'], d['tan_angle_sq'], surface_id=sid),
    9:  lambda d, sid: RPP(d['min_x'], d['max_x'], d['min_y'], d['max_y'],
                            d['min_z'], d['max_z'], surface_id=sid),
    10: lambda d, sid: Quadric(*d['coeffs'], surface_id=sid),
    11: lambda d, sid: XTorus(d['center_x'], d['center_y'], d['center_z'],
                               d['major_radius'], d['minor_radius'], surface_id=sid),
    12: lambda d, sid: YTorus(d['center_x'], d['center_y'], d['center_z'],
                               d['major_radius'], d['minor_radius'], surface_id=sid),
    13: lambda d, sid: ZTorus(d['center_x'], d['center_y'], d['center_z'],
                               d['major_radius'], d['minor_radius'], surface_id=sid),
    14: lambda d, sid: RCC(d['base_x'], d['base_y'], d['base_z'],
                            d['height_x'], d['height_y'], d['height_z'],
                            d['radius'], surface_id=sid),
    15: lambda d, sid: Box(d['corner_x'], d['corner_y'], d['corner_z'],
                            d['v1_x'], d['v1_y'], d['v1_z'],
                            d['v2_x'], d['v2_y'], d['v2_z'],
                            d['v3_x'], d['v3_y'], d['v3_z'], surface_id=sid),
    17: lambda d, sid: TRC(d['base_x'], d['base_y'], d['base_z'],
                            d['height_x'], d['height_y'], d['height_z'],
                            d['base_radius'], d['top_radius'], surface_id=sid),
    18: lambda d, sid: ELL(d['v1_x'], d['v1_y'], d['v1_z'],
                            d['v2_x'], d['v2_y'], d['v2_z'],
                            d['major_axis_len'], surface_id=sid),
    19: lambda d, sid: REC(d['base_x'], d['base_y'], d['base_z'],
                            d['height_x'], d['height_y'], d['height_z'],
                            d['axis1_x'], d['axis1_y'], d['axis1_z'],
                            d['axis2_x'], d['axis2_y'], d['axis2_z'], surface_id=sid),
    20: lambda d, sid: WED(d['vertex_x'], d['vertex_y'], d['vertex_z'],
                            d['v1_x'], d['v1_y'], d['v1_z'],
                            d['v2_x'], d['v2_y'], d['v2_z'],
                            d['v3_x'], d['v3_y'], d['v3_z'], surface_id=sid),
}


def _build_surface_from_node(sys, surface_id: int) -> Optional[Surface]:
    """Reconstruct a Python Surface object from the C system by surface ID."""
    node_id = sys.surface_node(surface_id, 1)
    if node_id is None:
        node_id = sys.surface_node(surface_id, -1)
    if node_id is None:
        return None
    data = sys.node_primitive_data(node_id)
    builder = _PRIMITIVE_BUILDERS.get(data.get('type'))
    if builder is None:
        return None
    return builder(data, surface_id)


def _populate_regions(model: Model, sys) -> None:
    """Attach region wrappers from the C system into the model.

    Surfaces are NOT materialized here: doing so would force every cell's
    CSG tree through `node_tree()` and `_unpack_csg_tree`, which dominates
    load time for large MCNP files. `Model.surfaces` lazily reconstructs
    the dict on first access.
    """
    sys.build_universe_index()
    for cell_info in sys.get_cells():
        cell_id = cell_info['cell_id']
        region = _ImportedRegion(sys, cell_info['root_node'])
        model._regions[cell_id] = region


def write_mcnp(model: Model, filename: Union[str, Path]) -> None:
    """Write model to MCNP format.

    Args:
        model: Model to export.
        filename: Output file path.
    """
    model.export_mcnp(filename)


def write_openmc(model: Model, filename: Union[str, Path]) -> None:
    """Write model to OpenMC XML format.

    Args:
        model: Model to export.
        filename: Output file path.
    """
    model.export_openmc(filename)


def write_serpent(model: Model, filename: Union[str, Path]) -> None:
    """Write model to Serpent input format.

    Args:
        model: Model to export.
        filename: Output file path.
    """
    model.export_serpent(filename)


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


_PlaceholderRegion = _ImportedRegion
