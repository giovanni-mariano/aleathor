# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""
Model, Cell, Material, and Universe classes for geometry definition.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple, Set, Union, TYPE_CHECKING, overload
from pathlib import Path
import math

if TYPE_CHECKING:
    from .geometry import Region
    from .surfaces import Surface

# Import the C extension
try:
    import _alea
except ImportError:
    _alea = None

# Import collection classes
from .collections import CellCollection, Cell, TraceResult, TraceSegment, VoidResult


_ELEMENT_SYMBOLS = (
    None,
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
)


def _element_symbol(z: int) -> str:
    if 0 < z < len(_ELEMENT_SYMBOLS):
        return _ELEMENT_SYMBOLS[z] or f"Z{z}"
    return f"Z{z}"


class Material:
    """Material definition backed by the C system.

    Provides methods to define nuclide/element composition and query
    material properties. All mutations are pushed to C immediately.

    Example:
        mat = model.add_material(1, name="Steel")
        mat.add_element(26, 0.70)   # Fe 70%
        mat.add_element(24, 0.18)   # Cr 18%
        mat.add_element(28, 0.12)   # Ni 12%
        mat.density = 7.8
        mat.weight_fractions = True
    """

    def __init__(self, model: 'Model', mat_index: int, material_id: int,
                 name: Optional[str] = None):
        self._model = model
        self._index = mat_index
        self._id = material_id
        self._name = name

    @property
    def id(self) -> int:
        """Material MCNP ID."""
        return self._id

    @property
    def name(self) -> Optional[str]:
        """Material name."""
        return self._name

    @name.setter
    def name(self, value: Optional[str]) -> None:
        self._name = value

    @property
    def density(self) -> Optional[float]:
        """Material standard density, or None if not set."""
        return self._model._sys.material_get_density(self._index)

    @density.setter
    def density(self, value: float) -> None:
        self._model._sys.material_set_density(self._index, value)

    @property
    def weight_fractions(self) -> bool:
        """True if fractions are weight fractions, False for atom fractions."""
        return self._model._sys.material_is_weight_fraction(self._index)

    @weight_fractions.setter
    def weight_fractions(self, value: bool) -> None:
        self._model._sys.material_set_weight_fraction(self._index, value)

    def add_nuclide(self, zaid: int, fraction: float,
                    library: Optional[str] = None) -> 'Material':
        """Add a nuclide to this material.

        Args:
            zaid: ZAID (ZZAAA format, e.g. 92235 for U-235)
            fraction: Atom or weight fraction (positive)
            library: Cross-section library suffix (e.g. ".80c"), or None

        Returns:
            self (for chaining)
        """
        self._model._sys.material_add_nuclide(self._index, zaid, library, fraction)
        return self

    def add_element(self, Z: int, fraction: float,
                    library: Optional[str] = None) -> 'Material':
        """Add an element with natural isotopic composition.

        Args:
            Z: Atomic number (1-118)
            fraction: Atom or weight fraction (positive)
            library: Library suffix for generated nuclides, or None

        Returns:
            self (for chaining)
        """
        self._model._sys.material_add_element(self._index, Z, library, fraction)
        return self

    def expand_elements(self) -> 'Material':
        """Expand element entries to explicit nuclides (natural abundances).

        Returns:
            self (for chaining)
        """
        self._model._sys.material_expand_elements(self._index)
        return self

    @property
    def nuclides(self) -> List[Dict[str, Any]]:
        """List of nuclides: [{zaid, library, fraction}, ...]."""
        return self._model._sys.material_get_nuclides(self._index)

    @property
    def elements(self) -> List[Dict[str, Any]]:
        """List of elements: [{Z, library, fraction}, ...]."""
        return self._model._sys.material_get_elements(self._index)

    @staticmethod
    def _format_composition_entry(kind: str, item: Dict[str, Any]) -> str:
        if kind == "nuclide":
            zaid = int(item["zaid"])
            z = zaid // 1000
            a = zaid % 1000
            label = f"{_element_symbol(z)}{a}"
        else:
            label = f"{_element_symbol(int(item['Z']))}0"

        return f"{label}={item['fraction']:g}"

    def _format_composition(self, limit: int = 6) -> str:
        entries = [
            self._format_composition_entry("nuclide", item)
            for item in self.nuclides
        ]
        entries.extend(
            self._format_composition_entry("element", item)
            for item in self.elements
        )

        if len(entries) > limit:
            shown = entries[:limit]
            shown.append(f"... +{len(entries) - limit} more")
            entries = shown

        return "[" + ", ".join(entries) + "]"

    def __repr__(self) -> str:
        fields = [f"id={self._id}"]
        if self._name:
            fields.append(f"name={self._name!r}")
        d = self.density
        if d is not None:
            fields.append(f"density={d:g}")
        fraction_kind = "weight" if self.weight_fractions else "atom"
        fields.append(f"fractions={fraction_kind!r}")
        nn = self._model._sys.material_nuclide_count(self._index)
        ne = self._model._sys.material_element_count(self._index)
        if nn > 0:
            fields.append(f"nuclides={nn}")
        if ne > 0:
            fields.append(f"elements={ne}")
        if nn > 0 or ne > 0:
            fields.append(f"composition={self._format_composition()}")
        return f"Material({', '.join(fields)})"

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Material):
            return self._model is other._model and self._id == other._id
        return False

    def __hash__(self) -> int:
        return hash((id(self._model), self._id))


@dataclass
class _CellData:
    """Internal cell data storage.

    A cell is a region of space with a specific material assignment.
    This is the internal dataclass; users interact via :class:`Cell`
    (in ``collections.py``).

    Attributes:
        id: Cell number (must be positive).
        region: CSG region defining the cell geometry.
        material: Material number (0 for void).
        density: Material density (always positive).
        density_unit: Density unit — ``"g/cm3"`` or ``"atoms/b-cm"``.
        name: Optional cell name.
        universe: Universe this cell belongs to.
        fill: Universe to fill this cell with.
        fill_transform: Transform to apply to fill.
        importance: Importance for particle transport.
    """
    id: int
    region: 'Region'
    material: int = 0
    density: float = 0.0
    density_unit: str = "g/cm3"
    name: Optional[str] = None
    universe: int = 0
    fill: Optional[int] = None
    fill_transform: Optional[int] = None
    importance: float = 1.0

    def __post_init__(self):
        if self.id <= 0:
            raise ValueError("Cell ID must be positive")

    @property
    def is_void(self) -> bool:
        """Check if cell is void (no material)."""
        return self.material == 0

    @property
    def is_filled(self) -> bool:
        """Check if cell is filled with another universe."""
        return self.fill is not None


class SliceApi:
    """Slice operations for a model."""

    def __init__(self, model: 'Model'):
        self._model = model

    @staticmethod
    def _axis_value(axis: Optional[str], value: Optional[float]) -> Tuple[Optional[str], Optional[float]]:
        if axis is None:
            return None, value
        axis = axis.lower()
        if axis not in ("x", "y", "z"):
            raise ValueError("axis must be 'x', 'y', or 'z'")
        if value is None:
            raise ValueError("value is required when axis is specified")
        return axis, value

    def curves(self, *, axis: Optional[str] = None, value: Optional[float] = None,
               origin: Optional[Tuple[float, float, float]] = None,
               normal: Optional[Tuple[float, float, float]] = None,
               up: Optional[Tuple[float, float, float]] = None,
               bounds: Optional[Tuple[float, float, float, float]] = None) -> Dict[str, Any]:
        """Get analytical curves for an axis-aligned or arbitrary slice."""
        model = self._model
        model._ensure_query_caches()
        if bounds is None:
            raise ValueError("bounds is required")

        axis, value = self._axis_value(axis, value)
        has_axis = axis is not None
        has_plane = origin is not None or normal is not None
        if int(has_axis) + int(has_plane) != 1:
            raise ValueError("Specify either axis+value or origin+normal")

        if axis == "z":
            x_min, x_max, y_min, y_max = bounds
            return model._sys.get_slice_curves_z(value, x_min, x_max, y_min, y_max)
        if axis == "y":
            x_min, x_max, z_min, z_max = bounds
            return model._sys.get_slice_curves_y(value, x_min, x_max, z_min, z_max)
        if axis == "x":
            y_min, y_max, z_min, z_max = bounds
            return model._sys.get_slice_curves_x(value, y_min, y_max, z_min, z_max)

        if origin is None or normal is None:
            raise ValueError("origin and normal are required for an arbitrary slice")
        if up is None:
            up = (0, 0, 1)
        u_min, u_max, v_min, v_max = bounds
        return model._sys.get_slice_curves(
            origin, normal, up, u_min, u_max, v_min, v_max
        )

    def grid(self, *, axis: Optional[str] = None, value: Optional[float] = None,
             origin: Optional[Tuple[float, float, float]] = None,
             normal: Optional[Tuple[float, float, float]] = None,
             up: Optional[Tuple[float, float, float]] = None,
             bounds: Optional[Tuple[float, float, float, float]] = None,
             resolution: Tuple[int, int] = (100, 100),
             universe_depth: int = -1,
             detect_errors: bool = False) -> Dict[str, Any]:
        """Sample cells/materials on an axis-aligned or arbitrary slice grid."""
        model = self._model
        model._ensure_query_caches()
        if bounds is None:
            raise ValueError("bounds is required")

        axis, value = self._axis_value(axis, value)
        has_axis = axis is not None
        has_plane = origin is not None or normal is not None
        if int(has_axis) + int(has_plane) != 1:
            raise ValueError("Specify either axis+value or origin+normal")

        if axis == "z":
            x_min, x_max, y_min, y_max = bounds
            nx, ny = resolution
            return model._sys.find_cells_grid_z(
                value, x_min, x_max, y_min, y_max, nx, ny,
                universe_depth=universe_depth,
                detect_errors=detect_errors,
            )
        if axis == "y":
            x_min, x_max, z_min, z_max = bounds
            nx, nz = resolution
            return model._sys.find_cells_grid_y(
                value, x_min, x_max, z_min, z_max, nx, nz,
                universe_depth=universe_depth,
                detect_errors=detect_errors,
            )
        if axis == "x":
            y_min, y_max, z_min, z_max = bounds
            ny, nz = resolution
            return model._sys.find_cells_grid_x(
                value, y_min, y_max, z_min, z_max, ny, nz,
                universe_depth=universe_depth,
                detect_errors=detect_errors,
            )

        if origin is None or normal is None:
            raise ValueError("origin and normal are required for an arbitrary slice")
        if up is None:
            up = (0, 0, 1)
        u_min, u_max, v_min, v_max = bounds
        nu, nv = resolution
        return model._sys.find_cells_grid(
            origin, normal, up, u_min, u_max, v_min, v_max, nu, nv,
            universe_depth=universe_depth,
            detect_errors=detect_errors,
        )

    def labels(self, grid_result, min_pixels=100, by_material=False):
        """Find optimal cell/material label positions for a slice grid."""
        from .slicing import find_label_positions
        return find_label_positions(self._model, grid_result, min_pixels, by_material)

    def surface_labels(self, grid_result, margin=20):
        """Find optimal surface label positions for a slice grid."""
        from .slicing import find_surface_label_positions
        return find_surface_label_positions(self._model, grid_result, margin)

    def check_overlaps(self, grid_result, universe_depth=-1):
        """Run comprehensive overlap detection on a slice grid."""
        from .slicing import check_grid_overlaps
        return check_grid_overlaps(self._model, grid_result, universe_depth)


class MeshApi:
    """Structured mesh operations for a model."""

    def __init__(self, model: 'Model'):
        self._model = model

    def export(self, filename: Union[str, Path], nx: int = 10, ny: int = 10,
               nz: int = 10,
               bounds: Optional[Tuple[float, float, float, float, float, float]] = None,
               format: str = "gmsh", void_material_id: int = 0,
               auto_pad: float = 0.01) -> None:
        """Export geometry as a structured mesh."""
        model = self._model
        model._ensure_sys()
        kwargs = dict(filename=str(filename), nx=nx, ny=ny, nz=nz,
                      format=format, void_material_id=void_material_id,
                      auto_pad=auto_pad)
        if bounds is not None:
            kwargs.update(x_min=bounds[0], x_max=bounds[1],
                          y_min=bounds[2], y_max=bounds[3],
                          z_min=bounds[4], z_max=bounds[5])
        model._sys.mesh_export(**kwargs)

    def sample(self, nx: int = 10, ny: int = 10, nz: int = 10,
               bounds: Optional[Tuple[float, float, float, float, float, float]] = None,
               void_material_id: int = 0,
               auto_pad: float = 0.01) -> Dict[str, Any]:
        """Sample geometry on a structured mesh."""
        model = self._model
        model._ensure_sys()
        kwargs = dict(nx=nx, ny=ny, nz=nz,
                      void_material_id=void_material_id,
                      auto_pad=auto_pad)
        if bounds is not None:
            kwargs.update(x_min=bounds[0], x_max=bounds[1],
                          y_min=bounds[2], y_max=bounds[3],
                          z_min=bounds[4], z_max=bounds[5])
        return model._sys.mesh_sample(**kwargs)


class IdApi:
    """ID renumbering and offset operations for a model."""

    def __init__(self, model: 'Model'):
        self._model = model

    def renumber_cells(self, start_id=1):
        self._model._ensure_sys()
        return self._model._sys.renumber_cells(start_id)

    def renumber_surfaces(self, start_id=1):
        self._model._ensure_sys()
        return self._model._sys.renumber_surfaces(start_id)

    def offset_cells(self, offset):
        self._model._ensure_sys()
        self._model._sys.offset_cell_ids(offset)

    def offset_surfaces(self, offset):
        self._model._ensure_sys()
        self._model._sys.offset_surface_ids(offset)

    def offset_materials(self, offset):
        self._model._ensure_sys()
        self._model._sys.offset_material_ids(offset)


class RepairApi:
    """Geometry cleanup and simplification operations for a model."""

    def __init__(self, model: 'Model'):
        self._model = model

    def flatten_universe(self, universe_id: int) -> None:
        self._model._ensure_sys()
        self._model._sys.flatten_universe(universe_id)

    def simplify(self) -> Dict[str, Any]:
        self._model._ensure_sys()
        return self._model._sys.flatten_all_cells()

    def split_union_cells(self):
        self._model._ensure_sys()
        return self._model._sys.split_union_cells()

    def expand_macrobodies(self):
        self._model._ensure_sys()
        return self._model._sys.expand_macrobodies()

    def tighten_bboxes(self, tolerance=1.0):
        self._model._ensure_sys()
        return self._model._sys.tighten_all_bboxes(tolerance)

    def tighten_cell_bbox(self, cell_id, tolerance=1.0):
        self._model._ensure_sys()
        idx = self._model._sys.cell_find(cell_id)
        if idx is None:
            raise KeyError(f"Cell {cell_id} not found")
        return self._model._sys.tighten_cell_bbox(idx, tolerance)

    def tighten_bbox_numerical(self, cell_id: int) -> None:
        self._model._ensure_sys()
        idx = self._model._sys.cell_find(cell_id)
        if idx is None:
            raise KeyError(f"Cell {cell_id} not found")
        self._model._sys.tighten_cell_bbox_numerical(idx)


class VoidApi:
    """Void generation operations for a model."""

    def __init__(self, model: 'Model'):
        self._model = model

    def generate(self,
                 bounds: Optional[Tuple[float, float, float, float, float, float]] = None,
                 max_depth: int = 8,
                 min_size: float = 0.1,
                 probes_per_axis: int = 3,
                 region: Optional['Region'] = None) -> 'VoidResult':
        """Find empty space in the geometry using octree decomposition."""
        if _alea is None:
            raise RuntimeError("C extension not available")
        model = self._model
        model._ensure_sys()

        if bounds is not None and region is not None:
            raise ValueError("bounds and region are mutually exclusive")

        bounds_region = None
        if region is not None:
            from .geometry import Region
            if not isinstance(region, Region):
                raise TypeError("region must be a geometry Region")
            bounds_region = region._to_csg(model)

        c_result = _alea.generate_void(
            model._sys,
            bounds=bounds,
            max_depth=max_depth,
            min_size=min_size,
            probes_per_axis=probes_per_axis,
            bounds_region=bounds_region,
        )
        return VoidResult(model, c_result)

    def add(self, voids: 'VoidResult') -> int:
        """Commit generated void regions as material-0 cells."""
        model = self._model
        model._ensure_sys()
        if not isinstance(voids, VoidResult):
            raise TypeError("add() expects a VoidResult")
        if voids._model is not model:
            raise ValueError("VoidResult belongs to a different Model")
        return voids._c_result.add_cells()

    def add_graveyard(self, voids: 'VoidResult') -> int:
        """Add a graveyard cell enclosing generated void bounds."""
        model = self._model
        model._ensure_sys()
        if not isinstance(voids, VoidResult):
            raise TypeError("add_graveyard() expects a VoidResult")
        if voids._model is not model:
            raise ValueError("VoidResult belongs to a different Model")
        return voids._c_result.add_graveyard()


class AnalysisApi:
    """Diagnostics and analysis operations for a model."""

    def __init__(self, model: 'Model'):
        self._model = model

    def find_overlaps(self, max_pairs: int = 100) -> List[Tuple[Cell, Cell]]:
        model = self._model
        model._ensure_sys()
        pairs = model._sys.find_overlaps(max_pairs)
        return [(Cell(model, idx1), Cell(model, idx2)) for idx1, idx2 in pairs]

    def bounding_sphere(self, tolerance=1.0):
        self._model._ensure_sys()
        return self._model._sys.compute_bounding_sphere(tolerance)

    def estimate_cell_volumes(self, n_rays=100000, center=None, radius=None):
        model = self._model
        model._ensure_query_caches()
        if center is None or radius is None:
            cx, cy, cz, r = self.bounding_sphere()
            if center is None:
                center = (cx, cy, cz)
            if radius is None:
                radius = r
        ox, oy, oz = center
        return model._sys.estimate_cell_volumes(ox, oy, oz, radius, n_rays)

    def estimate_instance_volumes(self, n_rays=100000):
        self._model._ensure_query_caches()
        return self._model._sys.estimate_instance_volumes(n_rays)

    def remove_cells_by_volume(self, volumes, threshold):
        self._model._ensure_sys()
        return self._model._sys.remove_cells_by_volume(volumes, threshold)


class BackendApi:
    """Low-level backend controls for a model."""

    def __init__(self, model: 'Model'):
        self._model = model

    @property
    def config(self) -> Dict[str, Any]:
        self._model._ensure_sys()
        return self._model._sys.get_config()

    @config.setter
    def config(self, settings: Dict[str, Any]) -> None:
        self._model._ensure_sys()
        self._model._sys.set_config(settings)

    def set_verbose(self, enabled: bool) -> None:
        self._model._ensure_sys()
        self._model._sys.set_verbose(enabled)

    def build_spatial_index(self) -> 'Model':
        self._model._ensure_sys()
        self._model._sys.build_spatial_index()
        return self._model

    @property
    def spatial_index_instance_count(self) -> int:
        if self._model._sys is None:
            return 0
        try:
            return self._model._sys.spatial_index_instance_count
        except AttributeError:
            return 0


@dataclass
class Universe:
    """Universe definition.

    A universe is a collection of cells that can be placed inside
    other cells via the FILL mechanism.

    Attributes:
        id: Universe number.
        name: Optional universe name.
        cells: Cells in this universe.
    """
    id: int
    name: Optional[str] = None
    cells: List[_CellData] = field(default_factory=list)

    def add_cell(self, cell: _CellData) -> None:
        """Add a cell to this universe."""
        cell.universe = self.id
        self.cells.append(cell)


class Model:
    """Complete geometry model.

    The Model class is the main container for a CSG geometry. It holds
    cells, materials, universes, and provides methods for export to
    various formats.

    Example:
        model = Model()

        # Add surfaces (optional - can use inline)
        fuel_surf = Sphere(0, 0, 0, radius=0.5)

        # Add cells
        model.add_cell(
            region=-fuel_surf,  # Inside sphere
            material=1,
            density=10.5,
            name="fuel"
        )

        # Export
        model.export_mcnp("geometry.inp")
    """

    def __init__(self, title: str = "aleathor Model"):
        """
        Args:
            title: Model title (used in exports).
        """
        self.title = title
        self._materials: Dict[int, Material] = {}
        self._universes: Dict[int, Universe] = {}
        self._surfaces: Dict[int, 'Surface'] = {}

        # Python-side metadata (not stored in C)
        self._names: Dict[int, Optional[str]] = {}       # cell_id -> name
        self._importances: Dict[int, float] = {}          # cell_id -> importance
        self._regions: Dict[int, 'Region'] = {}           # cell_id -> Region (optional)

        # Material MCNP ID -> C material index mapping
        self._material_index_map: Dict[int, int] = {}
        self._mixture_names: Dict[int, str] = {}  # mixture_id -> name

        # Internal C system (created lazily)
        self._sys = None
        # Cache keyed by (surface_id, positive) since sense is stored on node
        self._halfspace_node_cache: Dict[Tuple[int, bool], int] = {}

    def _ensure_sys(self):
        """Ensure internal C system exists."""
        if self._sys is None:
            if _alea is None:
                raise RuntimeError("C extension not available")
            self._sys = _alea.System()

    def _ensure_query_caches(self):
        """Ensure C system exists and query acceleration caches are built.

        Called by query methods (point lookups, raycasts, slices, volume
        estimation). The C entry point is idempotent: once a cache bit is
        set, subsequent calls reduce to atomic checks.
        """
        self._ensure_sys()
        self._sys.prepare_query_acceleration()

    def _is_mixture(self, material_id: int) -> bool:
        """Check if a material ID refers to a mixture."""
        if self._sys is None:
            return False
        return self._sys.find_mixture_by_id(material_id) is not None

    def _ensure_material(self, material_id: int) -> int:
        """Ensure a material is registered in the C system.

        Returns the C material index for the given MCNP material ID.
        Returns -1 for void (material_id == 0) or mixtures.
        """
        if material_id == 0:
            return -1
        # Mixtures are separate — don't register as material
        if self._is_mixture(material_id):
            return -1
        if material_id not in self._material_index_map:
            self._ensure_sys()
            self._material_index_map[material_id] = self._sys.add_material(material_id)
        return self._material_index_map[material_id]

    def _get_or_create_halfspace_node(self, surface: 'Surface', positive: bool) -> int:
        """Get or create a CSG node for a halfspace (surface + sense).

        Uses the csg_thor_*_surface functions which create surfaces with both
        positive and negative halfspace nodes, properly registered for raycast.
        """
        self._ensure_sys()

        cache_key = (surface.id, positive)
        if cache_key in self._halfspace_node_cache:
            return self._halfspace_node_cache[cache_key]

        # Check if surface already created (with opposite sense)
        opposite_key = (surface.id, not positive)
        if opposite_key in self._halfspace_node_cache:
            # Surface already exists, get nodes from cache
            # (both were cached when surface was created)
            return self._halfspace_node_cache[cache_key]

        # Create surface using *_surface functions (returns index, pos_node, neg_node)
        from .surfaces import (
            Plane, Sphere, XCylinder, YCylinder, ZCylinder,
            XCone, YCone, ZCone, Box, RCC,
            XTorus, YTorus, ZTorus, Quadric, TRC, ELL, REC, WED, RHP, GeneralBox
        )

        if isinstance(surface, Sphere):
            _, pos_node, neg_node = self._sys.sphere_surface(
                surface.id, surface.x0, surface.y0, surface.z0, surface.radius
            )
        elif isinstance(surface, ZCylinder):
            _, pos_node, neg_node = self._sys.cylinder_z_surface(
                surface.id, surface.x0, surface.y0, surface.radius
            )
        elif isinstance(surface, XCylinder):
            _, pos_node, neg_node = self._sys.cylinder_x_surface(
                surface.id, surface.y0, surface.z0, surface.radius
            )
        elif isinstance(surface, YCylinder):
            _, pos_node, neg_node = self._sys.cylinder_y_surface(
                surface.id, surface.x0, surface.z0, surface.radius
            )
        elif isinstance(surface, ZCone):
            _, pos_node, neg_node = self._sys.cone_z_surface(
                surface.id, surface.x0, surface.y0, surface.z0, surface.t_sq
            )
        elif isinstance(surface, XCone):
            _, pos_node, neg_node = self._sys.cone_x_surface(
                surface.id, surface.x0, surface.y0, surface.z0, surface.t_sq
            )
        elif isinstance(surface, YCone):
            _, pos_node, neg_node = self._sys.cone_y_surface(
                surface.id, surface.x0, surface.y0, surface.z0, surface.t_sq
            )
        elif isinstance(surface, ZTorus):
            _, pos_node, neg_node = self._sys.torus_z_surface(
                surface.id, surface.x0, surface.y0, surface.z0,
                surface.major_radius, surface.minor_radius
            )
        elif isinstance(surface, XTorus):
            _, pos_node, neg_node = self._sys.torus_x_surface(
                surface.id, surface.x0, surface.y0, surface.z0,
                surface.major_radius, surface.minor_radius
            )
        elif isinstance(surface, YTorus):
            _, pos_node, neg_node = self._sys.torus_y_surface(
                surface.id, surface.x0, surface.y0, surface.z0,
                surface.major_radius, surface.minor_radius
            )
        elif isinstance(surface, Box):
            _, pos_node, neg_node = self._sys.box_surface(
                surface.id,
                surface.xmin, surface.xmax,
                surface.ymin, surface.ymax,
                surface.zmin, surface.zmax
            )
        elif isinstance(surface, RCC):
            _, pos_node, neg_node = self._sys.rcc_surface(
                surface.id,
                surface.base_x, surface.base_y, surface.base_z,
                surface.height_x, surface.height_y, surface.height_z,
                surface.radius
            )
        elif isinstance(surface, TRC):
            _, pos_node, neg_node = self._sys.trc_surface(
                surface.id,
                surface.base_x, surface.base_y, surface.base_z,
                surface.height_x, surface.height_y, surface.height_z,
                surface.base_radius, surface.top_radius
            )
        elif isinstance(surface, ELL):
            _, pos_node, neg_node = self._sys.ell_surface(
                surface.id,
                surface.v1_x, surface.v1_y, surface.v1_z,
                surface.v2_x, surface.v2_y, surface.v2_z,
                surface.major_axis_len
            )
        elif isinstance(surface, REC):
            _, pos_node, neg_node = self._sys.rec_surface(
                surface.id,
                surface.base_x, surface.base_y, surface.base_z,
                surface.height_x, surface.height_y, surface.height_z,
                surface.axis1_x, surface.axis1_y, surface.axis1_z,
                surface.axis2_x, surface.axis2_y, surface.axis2_z
            )
        elif isinstance(surface, WED):
            _, pos_node, neg_node = self._sys.wed_surface(
                surface.id,
                surface.vertex_x, surface.vertex_y, surface.vertex_z,
                surface.v1_x, surface.v1_y, surface.v1_z,
                surface.v2_x, surface.v2_y, surface.v2_z,
                surface.v3_x, surface.v3_y, surface.v3_z
            )
        elif isinstance(surface, RHP):
            _, pos_node, neg_node = self._sys.rhp_surface(
                surface.id,
                surface.base_x, surface.base_y, surface.base_z,
                surface.height_x, surface.height_y, surface.height_z,
                surface.r1_x, surface.r1_y, surface.r1_z,
                surface.r2_x, surface.r2_y, surface.r2_z,
                surface.r3_x, surface.r3_y, surface.r3_z
            )
        elif isinstance(surface, GeneralBox):
            _, pos_node, neg_node = self._sys.box_general_surface(
                surface.id,
                surface.corner_x, surface.corner_y, surface.corner_z,
                surface.v1_x, surface.v1_y, surface.v1_z,
                surface.v2_x, surface.v2_y, surface.v2_z,
                surface.v3_x, surface.v3_y, surface.v3_z
            )
        elif isinstance(surface, Quadric):
            _, pos_node, neg_node = self._sys.quadric_surface(
                surface.id,
                surface.A, surface.B, surface.C,
                surface.D, surface.E, surface.F,
                surface.G, surface.H, surface.J, surface.K
            )
        elif isinstance(surface, Plane):
            _, pos_node, neg_node = self._sys.plane_surface(
                surface.id, surface.a, surface.b, surface.c, -surface.d
            )
        else:
            raise NotImplementedError(f"Surface type {type(surface)} not yet supported")

        # Cache both nodes for this surface
        self._halfspace_node_cache[(surface.id, True)] = pos_node
        self._halfspace_node_cache[(surface.id, False)] = neg_node
        self._surfaces[surface.id] = surface

        return pos_node if positive else neg_node

    # =========================================================================
    # Cell Management
    # =========================================================================

    def add_cell(self, region: 'Region', cell_id: Optional[int] = None,
                 material: Union[int, Material] = 0, density: float = 0.0,
                 density_unit: str = "g/cm3",
                 name: Optional[str] = None, universe: int = 0,
                 fill: Optional[int] = None, importance: float = 1.0) -> Cell:
        """Add a cell to the model.

        Pushes the cell to the C backend immediately.

        Args:
            region: CSG region defining the cell.
            cell_id: Cell number (auto-assigned if None).
            material: Material object or material ID (0 for void).
            density: Material density (positive value).
            density_unit: Density unit — ``"g/cm3"`` (default) or
                ``"atoms/b-cm"``.
            name: Optional cell name.
            universe: Universe this cell belongs to.
            fill: Universe to fill this cell with.
            importance: Particle importance.

        Returns:
            Cell wrapper for the created cell.
        """
        if density_unit not in ("g/cm3", "atoms/b-cm"):
            raise ValueError(
                f"density_unit must be 'g/cm3' or 'atoms/b-cm', got {density_unit!r}"
            )
        if fill is not None and fill <= 0:
            raise ValueError("fill universe must be a positive integer, or None")

        # Accept Material object or int
        if isinstance(material, Material):
            material = material.id

        self._ensure_sys()

        # Auto-assign cell_id
        if cell_id is None:
            # Use C cell count + existing IDs to find next available
            existing_ids = set(self._names.keys())
            for i in range(self._sys.cell_count):
                existing_ids.add(self._sys.get_cell_id(i))
            cell_id = max(existing_ids, default=0) + 1

        # Register material in C
        mat_index = self._ensure_material(material)

        # Build CSG node from Python Region
        node_id = region._to_csg(self)

        # Signed density for C
        signed_density = -abs(density) if density_unit == "g/cm3" else abs(density)

        # Push to C
        index = self._sys.add_cell(
            cell_id, node_id,
            material_index=mat_index,
            density=signed_density,
            universe_id=universe
        )

        # Apply mixture if material is a mixture
        if isinstance(material, int) and material != 0 and self._is_mixture(material):
            self._sys.cell_set_mixture(index, material)

        # Apply fill
        if fill is not None:
            self._sys.set_fill(index, fill, 0)

        # Store Python-side metadata
        self._names[cell_id] = name
        self._importances[cell_id] = importance
        self._regions[cell_id] = region

        # Track surfaces
        for surface in region.get_surfaces():
            self._surfaces[surface.id] = surface

        return Cell(self, index)

    def get_cell(self, cell_id: int) -> Cell:
        """Get cell by ID.

        Args:
            cell_id: Cell ID (MCNP cell number)

        Returns:
            Cell for the cell

        Raises:
            KeyError: If cell not found
        """
        self._ensure_sys()
        idx = self._sys.cell_find(cell_id)
        if idx is None:
            raise KeyError(f"Cell {cell_id} not found")
        return Cell(self, idx)

    def __getitem__(self, cell_id: int) -> Cell:
        """Get cell by ID using subscript notation.

        Args:
            cell_id: Cell ID (MCNP cell number)

        Returns:
            Cell for the cell

        Raises:
            KeyError: If cell not found

        Example:
            cell = model[10]
            print(cell.material)
        """
        return self.get_cell(cell_id)

    def remove_cell(self, cell_id: int) -> None:
        """Remove cell by ID."""
        self._ensure_sys()
        idx = self._sys.cell_find(cell_id)
        if idx is not None:
            self._sys.cell_remove(idx)
            self._names.pop(cell_id, None)
            self._importances.pop(cell_id, None)
            self._regions.pop(cell_id, None)

    def update_cell(self, cell_id: int,
                    material: Optional[Union[int, Material]] = None,
                    density: Optional[float] = None,
                    density_unit: Optional[str] = None,
                    importance: Optional[float] = None) -> Cell:
        """Update cell properties directly in the C backend.

        Args:
            cell_id: Cell ID to update
            material: New Material object or ID (None to keep current)
            density: New density (None to keep current)
            density_unit: New density unit (None to keep current)
            importance: New importance (None to keep current)

        Returns:
            Fresh Cell reflecting the updated cell

        Raises:
            KeyError: If cell not found

        Example:
            model.update_cell(10, material=2, density=5.0)
        """
        self._ensure_sys()
        idx = self._sys.cell_find(cell_id)
        if idx is None:
            raise KeyError(f"Cell {cell_id} not found")

        if material is not None:
            if isinstance(material, Material):
                material = material.id
            mat_index = self._ensure_material(material)
            self._sys.cell_set_material(idx, mat_index)
        if density is not None or density_unit is not None:
            # Need current values to compute signed density
            info = self._sys.get_cell_by_index(idx)
            d = abs(density) if density is not None else abs(info['density'])
            if density_unit is not None:
                if density_unit not in ("g/cm3", "atoms/b-cm"):
                    raise ValueError(
                        f"density_unit must be 'g/cm3' or 'atoms/b-cm', got {density_unit!r}"
                    )
                is_mass = (density_unit == "g/cm3")
            else:
                is_mass = info.get('is_mass_density', True)
            signed = -d if is_mass else d
            self._sys.cell_set_density(idx, signed)
        if importance is not None:
            self._importances[cell_id] = importance

        return Cell(self, idx)

    @property
    def cells(self) -> CellCollection:
        """Get cell collection with filtering support.

        Returns a CellCollection that supports iteration, indexing, and filtering.

        Example:
            # Iterate all cells
            for cell in model.cells:
                print(cell.id, cell.material)

            # Filter by material
            fuel_cells = model.cells.by_material(92)

            # Filter by universe
            u3_cells = model.cells.by_universe(3)

            # Chain filters
            model.cells.by_material(92).by_universe(0)

            # Count
            len(model.cells)
        """
        self._ensure_sys()
        return CellCollection(self)

    @property
    def surfaces(self) -> Dict[int, 'Surface']:
        """Get all surfaces by ID.

        Returns a dictionary mapping surface IDs to Surface objects.
        """
        return dict(self._surfaces)

    @property
    def slice(self) -> SliceApi:
        """Slice operations."""
        return SliceApi(self)

    @property
    def mesh(self) -> MeshApi:
        """Structured mesh operations."""
        return MeshApi(self)

    @property
    def ids(self) -> IdApi:
        """ID renumbering and offset operations."""
        return IdApi(self)

    @property
    def repair(self) -> RepairApi:
        """Geometry cleanup and simplification operations."""
        return RepairApi(self)

    @property
    def void(self) -> VoidApi:
        """Void generation operations."""
        return VoidApi(self)

    @property
    def analysis(self) -> AnalysisApi:
        """Diagnostics and analysis operations."""
        return AnalysisApi(self)

    @property
    def backend(self) -> BackendApi:
        """Low-level backend controls."""
        return BackendApi(self)

    # =========================================================================
    # Material Management
    # =========================================================================

    def add_material(self, material_id: int, name: Optional[str] = None,
                     density: Optional[float] = None) -> Material:
        """Add a material to the model.

        Args:
            material_id: Material number (positive integer).
            name: Optional material name.
            density: Optional standard density (g/cm3 if positive).

        Returns:
            Material object for adding nuclides/elements.

        Example:
            mat = model.add_material(1, name="Steel")
            mat.add_element(26, 0.70)   # Fe
            mat.add_element(24, 0.18)   # Cr
            mat.add_element(28, 0.12)   # Ni
            mat.density = 7.8
        """
        if material_id <= 0:
            raise ValueError("Material ID must be positive")
        if material_id in self._materials:
            raise ValueError(f"Material {material_id} already exists")

        mat_index = self._ensure_material(material_id)
        mat = Material(self, mat_index, material_id, name=name)

        if density is not None:
            mat.density = density

        self._materials[material_id] = mat
        return mat

    def get_material(self, material_id: int) -> Material:
        """Get material by ID.

        Works for both Python-created and C-loaded materials.
        """
        if material_id in self._materials:
            return self._materials[material_id]

        # Try to find in C system
        self._ensure_sys()
        idx = self._sys.find_material_by_id(material_id)
        if idx is None:
            raise KeyError(f"Material {material_id} not found")
        mat = Material(self, idx, material_id)
        self._materials[material_id] = mat
        return mat

    @property
    def materials(self) -> List[Material]:
        """Get all materials.

        Returns materials from both Python-created and C-loaded sources.
        """
        self._ensure_sys()
        count = self._sys.material_count()
        result = []
        for i in range(count):
            mid = self._sys.material_get_id(i)
            if mid in self._materials:
                result.append(self._materials[mid])
            else:
                mat = Material(self, i, mid)
                self._materials[mid] = mat
                result.append(mat)
        return result

    def create_mixture(self, material_ids: List[int], fractions: List[float],
                       new_id: int = 0, name: Optional[str] = None) -> int:
        """Create a mixture of existing materials.

        Creates a new material as a weighted combination of existing materials.
        Fractions are normalized automatically. The mixture is stored in the
        C backend and will be exported as a combined material card.

        Args:
            material_ids: List of material IDs to mix
            fractions: List of weight fractions (will be normalized)
            new_id: ID for new mixture (0 for auto-assign)
            name: Optional name for the mixture

        Returns:
            Assigned mixture ID

        Example:
            mix_id = model.create_mixture([1, 2], [0.7, 0.3], name="fuel-clad-mix")
        """
        self._ensure_sys()
        mix_id = self._sys.create_mixture(material_ids, fractions, new_id)
        if name:
            self._mixture_names[mix_id] = name
        return mix_id

    # =========================================================================
    # Universe Management
    # =========================================================================

    def add_universe(self, universe_id: int, name: Optional[str] = None) -> Universe:
        """Add a universe to the model.

        Args:
            universe_id: Universe number.
            name: Optional universe name.

        Returns:
            The created Universe object.
        """
        if universe_id in self._universes:
            raise ValueError(f"Universe {universe_id} already exists")

        univ = Universe(id=universe_id, name=name)
        self._universes[universe_id] = univ
        return univ

    def get_universe(self, universe_id: int) -> Universe:
        """Get universe by ID."""
        if universe_id not in self._universes:
            raise KeyError(f"Universe {universe_id} not found")
        return self._universes[universe_id]

    @property
    def universes(self) -> List[Universe]:
        """Get all universes."""
        return list(self._universes.values())

    # =========================================================================
    # Queries
    # =========================================================================

    def cell_path_at(self, x: float, y: float, z: float) -> List[Cell]:
        """Return the containment path for a point across hierarchy depths.

        Returns cells at all levels of the universe hierarchy, ordered by
        depth (depth=0 outermost). This is useful for debugging hierarchical
        FILL geometries where a point is contained by parent cells before the
        terminal cell returned by ``cell_at``.

        Args:
            x, y, z: Point coordinates

        Returns:
            List of Cell objects, ordered by depth. Each Cell
            has its depth property set. Empty list if point is in void.

        Example:
            cells = model.cell_path_at(0, 0, 0)
            for cell in cells:
                print(f"depth={cell.depth}: cell {cell.id}, material {cell.material}")
        """
        self._ensure_sys()
        hits = self._sys.find_all_cells(x, y, z)
        result = []
        for hit in hits:
            cv = Cell(self, hit['cell_index'], _depth=hit['depth'])
            result.append(cv)
        return result

    def contains_point(self, x: float, y: float, z: float) -> bool:
        """Check if point is inside any cell (not void)."""
        return self.cell_at(x, y, z) is not None

    # =========================================================================
    # Export
    # =========================================================================

    def _expand_all_elements(self) -> None:
        """Expand element entries to nuclides for all materials before export."""
        count = self._sys.material_count()
        for i in range(count):
            if self._sys.material_element_count(i) > 0:
                self._sys.material_expand_elements(i)

    def export_mcnp(self, filename: Union[str, Path]) -> None:
        """Export model to MCNP format.

        Args:
            filename: Output file path.
        """
        self._ensure_sys()
        self._expand_all_elements()
        self._sys.export_mcnp(str(filename))

    def export_openmc(self, filename: Union[str, Path]) -> None:
        """Export model to OpenMC XML format.

        Args:
            filename: Output file path (geometry.xml).
        """
        self._ensure_sys()
        self._expand_all_elements()
        self._sys.export_openmc(str(filename))

    def export_serpent(self, filename: Union[str, Path]) -> None:
        """Export model to Serpent input format.

        Args:
            filename: Output file path.
        """
        self._ensure_sys()
        self._expand_all_elements()
        self._sys.export_serpent(str(filename))

    def to_mcnp_string(self) -> str:
        """Generate MCNP input as string.

        Returns:
            MCNP input file contents.
        """
        self._ensure_sys()
        # Use temporary file approach for now
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.inp', delete=False) as f:
            self.export_mcnp(f.name)
            f.flush()
            with open(f.name, 'r') as rf:
                content = rf.read()
            import os
            os.unlink(f.name)
            return content

    # =========================================================================
    # Mesh Export
    # =========================================================================

    # =========================================================================
    # Validation
    # =========================================================================

    def validate(self) -> List[str]:
        """Validate model for common issues.

        Returns:
            List of warning/error messages.
        """
        issues = []

        # Check for undefined materials in cells
        if self._sys:
            for i in range(self._sys.cell_count):
                info = self._sys.get_cell_by_index(i)
                mat_id = info.get('material_id', 0)
                if mat_id > 0 and mat_id not in self._materials:
                    issues.append(f"Cell {info['cell_id']} uses undefined material {mat_id}")

        # Check for overlaps
        overlaps = self.analysis.find_overlaps()
        for c1, c2 in overlaps:
            issues.append(f"Cells {c1.id} and {c2.id} overlap")

        # Check C library validation
        self._ensure_sys()
        n_issues = self._sys.validate()
        if n_issues > 0:
            issues.append(f"Internal validation found {n_issues} issues")

        return issues

    # =========================================================================
    # Utilities
    # =========================================================================

    def print_summary(self) -> None:
        """Print model summary."""
        print(str(self))

    def __str__(self) -> str:
        """Human-readable string representation."""
        n_cells = self._sys.cell_count if self._sys else 0
        n_surfaces = self._sys.surface_count if self._sys else len(self._surfaces)
        n_universes = len(self._universes) if self._universes else 0
        return f"Model: {n_cells} cells, {n_surfaces} surfaces, {n_universes} universes"

    def __repr__(self) -> str:
        n_cells = self._sys.cell_count if self._sys else 0
        n_surfaces = self._sys.surface_count if self._sys else len(self._surfaces)
        return f"Model(title='{self.title}', cells={n_cells}, surfaces={n_surfaces})"

    # =========================================================================
    # Cell Filtering (C-backed for performance)
    # =========================================================================

    def get_cells_by_material(self, material_id: int) -> List[int]:
        """Get indices of all cells with the given material ID.

        Args:
            material_id: Material ID to filter by

        Returns:
            List of cell indices
        """
        self._ensure_sys()
        return self._sys.get_cells_by_material(material_id)

    def get_cells_by_universe(self, universe_id: int) -> List[int]:
        """Get indices of all cells in the given universe.

        Args:
            universe_id: Universe ID to filter by

        Returns:
            List of cell indices
        """
        self._ensure_sys()
        return self._sys.get_cells_by_universe(universe_id)

    def get_cells_filling_universe(self, universe_id: int) -> List[int]:
        """Get indices of all cells that FILL the given universe.

        Returns cells where fill_universe == universe_id, i.e., cells that
        contain this universe as their fill content.

        Args:
            universe_id: Universe ID to search for

        Returns:
            List of cell indices
        """
        self._ensure_sys()
        if universe_id <= 0:
            raise ValueError("fill universe must be a positive integer")
        return self._sys.get_cells_filling_universe(universe_id)

    def get_cells_in_bbox(self, bounds: Tuple[float, float, float, float, float, float]) -> List[int]:
        """Get indices of cells whose bounding boxes intersect the given region.

        Args:
            bounds: (x_min, x_max, y_min, y_max, z_min, z_max)

        Returns:
            List of cell indices
        """
        self._ensure_query_caches()
        x_min, x_max, y_min, y_max, z_min, z_max = bounds
        return self._sys.get_cells_in_bbox(x_min, x_max, y_min, y_max, z_min, z_max)

    # =========================================================================
    # Extraction Operations
    # =========================================================================

    def extract_universe(self, universe_id: int) -> 'Model':
        """Extract a universe and all universes it references into a new Model.

        Args:
            universe_id: Universe ID to extract

        Returns:
            New Model with the extracted universe
        """
        self._ensure_sys()
        new_sys = self._sys.extract_universe(universe_id)
        new_model = Model(title=f"Universe {universe_id}")
        new_model._sys = new_sys
        return new_model

    def extract_region(self, bounds: Tuple[float, float, float, float, float, float]) -> 'Model':
        """Extract cells in a bounding box region into a new Model.

        Args:
            bounds: (x_min, x_max, y_min, y_max, z_min, z_max)

        Returns:
            New Model with the extracted cells
        """
        self._ensure_sys()
        x_min, x_max, y_min, y_max, z_min, z_max = bounds
        new_sys = self._sys.extract_region(x_min, x_max, y_min, y_max, z_min, z_max)
        new_model = Model(title=f"Region [{x_min:.1f},{x_max:.1f}]x[{y_min:.1f},{y_max:.1f}]x[{z_min:.1f},{z_max:.1f}]")
        new_model._sys = new_sys
        return new_model

    # =========================================================================
    # Design Document API - Queries
    # =========================================================================

    def cell_at(self, x: float, y: float, z: float) -> Optional[Cell]:
        """Find the cell containing a point.

        Args:
            x, y, z: Point coordinates

        Returns:
            Cell for the cell containing the point, or None if in void

        Example:
            cell = model.cell_at(0, 0, 0)
            if cell:
                print(f"Point is in cell {cell.id}, material {cell.material}")
        """
        self._ensure_query_caches()
        result = self._sys.find_cell(x, y, z)
        if result is None:
            return None
        cell_idx, _ = result
        return Cell(self, cell_idx)

    def trace(self, origin: Optional[Tuple[float, float, float]] = None,
              direction: Optional[Tuple[float, float, float]] = None,
              start: Optional[Tuple[float, float, float]] = None,
              end: Optional[Tuple[float, float, float]] = None,
              max_distance: float = 0,
              cell_aware: bool = False) -> TraceResult:
        """Trace a ray through the geometry.

        Can be called in two ways:
        1. With origin and direction: trace(origin=(0,0,0), direction=(1,0,0))
        2. With start and end points: trace(start=(0,0,0), end=(100,0,0))

        Args:
            origin: Ray origin point (x, y, z)
            direction: Ray direction vector (dx, dy, dz)
            start: Start point for point-to-point trace
            end: End point for point-to-point trace
            max_distance: Maximum trace distance (0 = infinite)
            cell_aware: If True, use cell-aware raycast (more efficient,
                tests only surfaces belonging to each cell).

        Returns:
            TraceResult containing segments along the ray

        Example:
            # Trace with direction
            result = model.trace(origin=(0,0,0), direction=(1,0,0))

            # Trace between points
            result = model.trace(start=(0,0,0), end=(100,0,0))

            # Cell-aware trace (faster)
            result = model.trace(origin=(0,0,0), direction=(1,0,0), cell_aware=True)

            # Iterate segments
            for seg in result:
                print(f"{seg.cell}: {seg.length} cm")
        """
        self._ensure_sys()

        # Determine origin and direction
        if start is not None and end is not None:
            # Point-to-point mode
            ox, oy, oz = start
            ex, ey, ez = end
            dx, dy, dz = ex - ox, ey - oy, ez - oz
            dist = math.sqrt(dx*dx + dy*dy + dz*dz)
            if dist < 1e-10:
                return TraceResult(self, [])
            # Normalize direction
            dx, dy, dz = dx/dist, dy/dist, dz/dist
            max_distance = dist
        elif origin is not None and direction is not None:
            ox, oy, oz = origin
            dx, dy, dz = direction
            # Normalize direction
            norm = math.sqrt(dx*dx + dy*dy + dz*dz)
            if norm < 1e-10:
                return TraceResult(self, [])
            dx, dy, dz = dx/norm, dy/norm, dz/norm
        else:
            raise ValueError("Must provide either (origin, direction) or (start, end)")

        self._ensure_query_caches()

        # Call C raycast
        if cell_aware:
            segments = self._sys.raycast_cell_aware(ox, oy, oz, dx, dy, dz, t_max=max_distance)
        else:
            segments = self._sys.raycast(ox, oy, oz, dx, dy, dz, max_distance)
        return TraceResult(self, segments)

    # =========================================================================
    # Design Document API - Export
    # =========================================================================

    def save(self, filename: Union[str, Path], format: Optional[str] = None) -> None:
        """Save model to file.

        Auto-detects format from extension or uses explicit format.

        Args:
            filename: Output file path
            format: Export format ('mcnp', 'openmc', or 'serpent').
                If None, detected from extension.

        Example:
            model.save("output.inp")              # MCNP (default)
            model.save("geometry.xml")            # OpenMC (from extension)
            model.save("geometry.serp")           # Serpent (from extension)
            model.save("output", format="mcnp")   # Explicit format
        """
        path = Path(filename)

        # Determine format
        if format is None:
            ext = path.suffix.lower()
            if ext in ('.xml',):
                format = 'openmc'
            elif ext in ('.serp', '.sss', '.serpent'):
                format = 'serpent'
            else:
                format = 'mcnp'

        if format == 'openmc':
            self.export_openmc(filename)
        elif format == 'serpent':
            self.export_serpent(filename)
        else:
            self.export_mcnp(filename)

    # =========================================================================
    # Design Document API - Plotting (delegates to plotting module)
    # =========================================================================

    def plot(self, z=None, y=None, x=None, **kwargs):
        """Plot geometry slice using grid-based rendering."""
        from .plotting import plot as _plot
        return _plot(self, z=z, y=y, x=x, **kwargs)

    def plot_views(self, bounds=None, save=None, **kwargs):
        """Plot three orthogonal views (XY, XZ, YZ)."""
        from .plotting import plot_views as _plot_views
        return _plot_views(self, bounds=bounds, save=save, **kwargs)
