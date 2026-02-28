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
from .collections import CellCollection, Cell, TraceResult, TraceSegment


@dataclass
class Material:
    """Material definition.

    Attributes:
        id: Material number (must be positive).
        name: Optional material name.
        density: Material density.
        composition: Optional composition dict {nuclide: fraction}.
    """
    id: int
    name: Optional[str] = None
    density: float = 0.0
    composition: Dict[str, float] = field(default_factory=dict)

    def __post_init__(self):
        if self.id <= 0:
            raise ValueError("Material ID must be positive")


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

    def __init__(self, title: str = "AleaTHOR Model"):
        """
        Args:
            title: Model title (used in exports).
        """
        self.title = title
        self._cells: Dict[int, _CellData] = {}
        self._materials: Dict[int, Material] = {}
        self._universes: Dict[int, Universe] = {}
        self._surfaces: Dict[int, 'Surface'] = {}

        # Internal C system (created lazily)
        self._sys = None
        # Cache keyed by (surface_id, positive) since sense is stored on node
        self._halfspace_node_cache: Dict[Tuple[int, bool], int] = {}
        self._dirty = True  # Track if model needs rebuild

    def _ensure_sys(self):
        """Ensure internal C system exists."""
        if self._sys is None:
            if _alea is None:
                raise RuntimeError("C extension not available")
            self._sys = _alea.System()

    def _rebuild_if_needed(self):
        """Rebuild internal system if dirty."""
        if not self._dirty:
            return

        self._ensure_sys()
        self._sys.reset()
        self._halfspace_node_cache.clear()

        # Register materials in C system and build id→index map
        material_index_map = {}  # MCNP material_id → C material_index
        for cell in self._cells.values():
            mat_id = cell.material
            if mat_id != 0 and mat_id not in material_index_map:
                material_index_map[mat_id] = self._sys.add_material(mat_id)

        # Build cells
        for cell in self._cells.values():
            node_id = cell.region._to_csg(self)
            # C expects signed density: negative = g/cm³, positive = atoms/b-cm
            signed_density = -cell.density if cell.density_unit == "g/cm3" else cell.density
            mat_index = material_index_map.get(cell.material, -1)  # -1 = void
            self._sys.add_cell(
                cell.id, node_id,
                material_index=mat_index,
                density=signed_density,
                universe_id=cell.universe
            )

        self._sys.build_universe_index()

        # Apply fills from Python Cells (add_cell doesn't pass fill to C)
        for cell in self._cells.values():
            if cell.fill is not None:
                idx = self._sys.cell_find(cell.id)
                if idx is not None:
                    transform = cell.fill_transform if cell.fill_transform is not None else 0
                    self._sys.set_fill(idx, cell.fill, transform)

        self._dirty = False

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
            Plane, Sphere, CylinderX, CylinderY, CylinderZ,
            ConeX, ConeY, ConeZ, Box, RCC,
            TorusX, TorusY, TorusZ, Quadric, TRC, ELL, REC, WED, RHP, GeneralBox
        )

        if isinstance(surface, Sphere):
            _, pos_node, neg_node = self._sys.sphere_surface(
                surface.id, surface.x0, surface.y0, surface.z0, surface.radius
            )
        elif isinstance(surface, CylinderZ):
            _, pos_node, neg_node = self._sys.cylinder_z_surface(
                surface.id, surface.x0, surface.y0, surface.radius
            )
        elif isinstance(surface, CylinderX):
            _, pos_node, neg_node = self._sys.cylinder_x_surface(
                surface.id, surface.y0, surface.z0, surface.radius
            )
        elif isinstance(surface, CylinderY):
            _, pos_node, neg_node = self._sys.cylinder_y_surface(
                surface.id, surface.x0, surface.z0, surface.radius
            )
        elif isinstance(surface, ConeZ):
            _, pos_node, neg_node = self._sys.cone_z_surface(
                surface.id, surface.x0, surface.y0, surface.z0, surface.t_sq
            )
        elif isinstance(surface, ConeX):
            _, pos_node, neg_node = self._sys.cone_x_surface(
                surface.id, surface.x0, surface.y0, surface.z0, surface.t_sq
            )
        elif isinstance(surface, ConeY):
            _, pos_node, neg_node = self._sys.cone_y_surface(
                surface.id, surface.x0, surface.y0, surface.z0, surface.t_sq
            )
        elif isinstance(surface, TorusZ):
            _, pos_node, neg_node = self._sys.torus_z_surface(
                surface.id, surface.x0, surface.y0, surface.z0,
                surface.major_radius, surface.minor_radius
            )
        elif isinstance(surface, TorusX):
            _, pos_node, neg_node = self._sys.torus_x_surface(
                surface.id, surface.x0, surface.y0, surface.z0,
                surface.major_radius, surface.minor_radius
            )
        elif isinstance(surface, TorusY):
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
                 material: int = 0, density: float = 0.0,
                 density_unit: str = "g/cm3",
                 name: Optional[str] = None, universe: int = 0,
                 fill: Optional[int] = None, importance: float = 1.0) -> _CellData:
        """Add a cell to the model.

        Args:
            region: CSG region defining the cell.
            cell_id: Cell number (auto-assigned if None).
            material: Material number (0 for void).
            density: Material density (positive value).
            density_unit: Density unit — ``"g/cm3"`` (default) or
                ``"atoms/b-cm"``.
            name: Optional cell name.
            universe: Universe this cell belongs to.
            fill: Universe to fill this cell with.
            importance: Particle importance.

        Returns:
            The created _CellData object.
        """
        if density_unit not in ("g/cm3", "atoms/b-cm"):
            raise ValueError(
                f"density_unit must be 'g/cm3' or 'atoms/b-cm', got {density_unit!r}"
            )

        if cell_id is None:
            cell_id = max(self._cells.keys(), default=0) + 1

        if cell_id in self._cells:
            raise ValueError(f"Cell {cell_id} already exists")

        cell = _CellData(
            id=cell_id,
            region=region,
            material=material,
            density=abs(density),
            density_unit=density_unit,
            name=name,
            universe=universe,
            fill=fill,
            importance=importance
        )

        self._cells[cell_id] = cell
        self._dirty = True

        # Track surfaces
        for surface in region.get_surfaces():
            self._surfaces[surface.id] = surface

        return cell

    def get_cell(self, cell_id: int) -> Cell:
        """Get cell by ID.

        Args:
            cell_id: Cell ID (MCNP cell number)

        Returns:
            Cell for the cell

        Raises:
            KeyError: If cell not found
        """
        self._rebuild_if_needed()
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
        if cell_id in self._cells:
            del self._cells[cell_id]
            self._dirty = True

    def update_cell(self, cell_id: int, material: Optional[int] = None,
                    density: Optional[float] = None,
                    density_unit: Optional[str] = None,
                    importance: Optional[float] = None) -> Cell:
        """Update cell properties.

        Modifies the Python Cell object and marks the model dirty so that
        the C system is rebuilt on the next query.

        Args:
            cell_id: Cell ID to update
            material: New material ID (None to keep current)
            density: New density (None to keep current)
            density_unit: New density unit (None to keep current)
            importance: New importance (None to keep current)

        Returns:
            Fresh Cell reflecting the updated cell

        Raises:
            KeyError: If cell not found in Python model

        Example:
            model.update_cell(10, material=2, density=5.0)
        """
        if cell_id not in self._cells:
            raise KeyError(f"Cell {cell_id} not found")

        cell = self._cells[cell_id]
        if material is not None:
            cell.material = material
        if density is not None:
            cell.density = abs(density)
        if density_unit is not None:
            if density_unit not in ("g/cm3", "atoms/b-cm"):
                raise ValueError(
                    f"density_unit must be 'g/cm3' or 'atoms/b-cm', got {density_unit!r}"
                )
            cell.density_unit = density_unit
        if importance is not None:
            cell.importance = importance

        self._dirty = True
        return self.get_cell(cell_id)

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
        self._rebuild_if_needed()
        return CellCollection(self)

    @property
    def surfaces(self) -> Dict[int, 'Surface']:
        """Get all surfaces by ID.

        Returns a dictionary mapping surface IDs to Surface objects.
        """
        return dict(self._surfaces)

    # =========================================================================
    # Material Management
    # =========================================================================

    def add_material(self, material_id: int, name: Optional[str] = None,
                     density: float = 0.0,
                     composition: Optional[Dict[str, float]] = None) -> Material:
        """Add a material to the model.

        Args:
            material_id: Material number.
            name: Optional material name.
            density: Default density.
            composition: Nuclide composition {nuclide: fraction}.

        Returns:
            The created Material object.
        """
        if material_id in self._materials:
            raise ValueError(f"Material {material_id} already exists")

        mat = Material(
            id=material_id,
            name=name,
            density=density,
            composition=composition or {}
        )

        self._materials[material_id] = mat
        return mat

    def get_material(self, material_id: int) -> Material:
        """Get material by ID."""
        if material_id not in self._materials:
            raise KeyError(f"Material {material_id} not found")
        return self._materials[material_id]

    @property
    def materials(self) -> List[Material]:
        """Get all materials."""
        return list(self._materials.values())

    def create_mixture(self, material_ids: List[int], fractions: List[float],
                       new_id: int = 0, name: Optional[str] = None) -> int:
        """Create a mixture of existing materials.

        Creates a new material as a weighted combination of existing materials.
        Fractions are normalized automatically. The mixture stores references
        to the original materials, which is preserved in MCNP export comments.

        Args:
            material_ids: List of material IDs to mix
            fractions: List of weight fractions (will be normalized)
            new_id: ID for new mixture (0 for auto-assign)
            name: Optional name for the mixture

        Returns:
            Assigned material ID

        Example:
            # Create 70% M1 + 30% M2 mixture
            mix_id = model.create_mixture([1, 2], [0.7, 0.3], name="fuel-clad-mix")
        """
        self._ensure_sys()
        self._rebuild_if_needed()
        mix_id = self._sys.create_mixture(material_ids, fractions, new_id)

        # Track in Python materials dict too
        if name:
            mat = Material(id=mix_id, name=name)
            self._materials[mix_id] = mat

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

    def cells_at(self, x: float, y: float, z: float) -> List[Cell]:
        """Find all cells containing a point across hierarchy depths.

        Returns cells at all levels of the universe hierarchy, ordered
        by depth (depth=0 outermost). This is essential for hierarchical
        FILL geometries where multiple cells contain the same point at
        different nesting levels.

        Args:
            x, y, z: Point coordinates

        Returns:
            List of Cell objects, ordered by depth. Each Cell
            has its depth property set. Empty list if point is in void.

        Example:
            cells = model.cells_at(0, 0, 0)
            for cell in cells:
                print(f"depth={cell.depth}: cell {cell.id}, material {cell.material}")
        """
        self._rebuild_if_needed()
        hits = self._sys.find_all_cells(x, y, z)
        result = []
        for hit in hits:
            cv = Cell(self, hit['cell_index'], _depth=hit['depth'])
            result.append(cv)
        return result

    def contains_point(self, x: float, y: float, z: float) -> bool:
        """Check if point is inside any cell (not void)."""
        return self.cell_at(x, y, z) is not None

    def find_overlaps(self, max_pairs: int = 100) -> List[Tuple[Cell, Cell]]:
        """Find overlapping cell pairs.

        Returns:
            List of (Cell, Cell) tuples for overlapping cells.
        """
        self._rebuild_if_needed()
        pairs = self._sys.find_overlaps(max_pairs)
        result = []
        for idx1, idx2 in pairs:
            result.append((Cell(self, idx1), Cell(self, idx2)))
        return result

    # =========================================================================
    # Export
    # =========================================================================

    def export_mcnp(self, filename: Union[str, Path]) -> None:
        """Export model to MCNP format.

        Args:
            filename: Output file path.
        """
        self._rebuild_if_needed()
        self._sys.export_mcnp(str(filename))

    def export_openmc(self, filename: Union[str, Path]) -> None:
        """Export model to OpenMC XML format.

        Args:
            filename: Output file path (geometry.xml).
        """
        self._rebuild_if_needed()
        self._sys.export_openmc(str(filename))

    def to_mcnp_string(self) -> str:
        """Generate MCNP input as string.

        Returns:
            MCNP input file contents.
        """
        self._rebuild_if_needed()
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

    def export_mesh(self, filename: Union[str, Path], nx: int = 10, ny: int = 10,
                    nz: int = 10, bounds: Optional[Tuple[float, float, float, float, float, float]] = None,
                    format: str = "gmsh", void_material_id: int = 0,
                    auto_pad: float = 0.01) -> None:
        """Export geometry as structured mesh (Gmsh or VTK).

        Args:
            filename: Output file path (.msh for Gmsh, .vtk for VTK)
            nx, ny, nz: Number of elements per axis
            bounds: (x_min, x_max, y_min, y_max, z_min, z_max) or None for auto
            format: 'gmsh' or 'vtk'
            void_material_id: Material ID for void regions (default 0)
            auto_pad: Fractional padding for auto-bounds (default 0.01)

        Example:
            model.export_mesh("output.msh", nx=50, ny=50, nz=50)
        """
        self._rebuild_if_needed()
        kwargs = dict(filename=str(filename), nx=nx, ny=ny, nz=nz,
                      format=format, void_material_id=void_material_id,
                      auto_pad=auto_pad)
        if bounds is not None:
            kwargs.update(x_min=bounds[0], x_max=bounds[1],
                          y_min=bounds[2], y_max=bounds[3],
                          z_min=bounds[4], z_max=bounds[5])
        self._sys.mesh_export(**kwargs)

    def sample_mesh(self, nx: int = 10, ny: int = 10, nz: int = 10,
                    bounds: Optional[Tuple[float, float, float, float, float, float]] = None,
                    void_material_id: int = 0,
                    auto_pad: float = 0.01) -> Dict[str, Any]:
        """Sample geometry on a structured mesh.

        Args:
            nx, ny, nz: Number of elements per axis
            bounds: (x_min, x_max, y_min, y_max, z_min, z_max) or None for auto
            void_material_id: Material ID for void regions (default 0)
            auto_pad: Fractional padding for auto-bounds (default 0.01)

        Returns:
            Dict with 'material_ids', 'cell_ids' (flat lists, Z-major order),
            'x_nodes', 'y_nodes', 'z_nodes' (node positions), 'nx', 'ny', 'nz'.

        Example:
            result = model.sample_mesh(nx=50, ny=50, nz=50)
            materials = result['material_ids']
        """
        self._rebuild_if_needed()
        kwargs = dict(nx=nx, ny=ny, nz=nz,
                      void_material_id=void_material_id,
                      auto_pad=auto_pad)
        if bounds is not None:
            kwargs.update(x_min=bounds[0], x_max=bounds[1],
                          y_min=bounds[2], y_max=bounds[3],
                          z_min=bounds[4], z_max=bounds[5])
        return self._sys.mesh_sample(**kwargs)

    # =========================================================================
    # Analytical Slice API
    # =========================================================================

    def get_slice_curves_z(self, z: float, bounds: Tuple[float, float, float, float]) -> Dict[str, Any]:
        """Get analytical curves for XY slice at z=constant."""
        self._rebuild_if_needed()
        x_min, x_max, y_min, y_max = bounds
        return self._sys.get_slice_curves_z(z, x_min, x_max, y_min, y_max)

    def get_slice_curves_y(self, y: float, bounds: Tuple[float, float, float, float]) -> Dict[str, Any]:
        """Get analytical curves for XZ slice at y=constant."""
        self._rebuild_if_needed()
        x_min, x_max, z_min, z_max = bounds
        return self._sys.get_slice_curves_y(y, x_min, x_max, z_min, z_max)

    def get_slice_curves_x(self, x: float, bounds: Tuple[float, float, float, float]) -> Dict[str, Any]:
        """Get analytical curves for YZ slice at x=constant."""
        self._rebuild_if_needed()
        y_min, y_max, z_min, z_max = bounds
        return self._sys.get_slice_curves_x(x, y_min, y_max, z_min, z_max)

    def get_slice_curves(self, origin: Tuple[float, float, float],
                         normal: Tuple[float, float, float],
                         up: Tuple[float, float, float],
                         bounds: Tuple[float, float, float, float]) -> Dict[str, Any]:
        """Get analytical curves for arbitrary plane slice."""
        self._rebuild_if_needed()
        u_min, u_max, v_min, v_max = bounds
        return self._sys.get_slice_curves(origin, normal, up, u_min, u_max, v_min, v_max)

    def find_cells_grid_z(self, z, bounds, resolution=(100, 100),
                          universe_depth=-1, detect_errors=False):
        """Sample cells on a grid for XY slice at z=constant."""
        self._rebuild_if_needed()
        x_min, x_max, y_min, y_max = bounds
        nx, ny = resolution
        return self._sys.find_cells_grid_z(z, x_min, x_max, y_min, y_max, nx, ny,
                                           universe_depth=universe_depth,
                                           detect_errors=detect_errors)

    def find_cells_grid_y(self, y, bounds, resolution=(100, 100),
                          universe_depth=-1, detect_errors=False):
        """Sample cells on a grid for XZ slice at y=constant."""
        self._rebuild_if_needed()
        x_min, x_max, z_min, z_max = bounds
        nx, nz = resolution
        return self._sys.find_cells_grid_y(y, x_min, x_max, z_min, z_max, nx, nz,
                                           universe_depth=universe_depth,
                                           detect_errors=detect_errors)

    def find_cells_grid_x(self, x, bounds, resolution=(100, 100),
                          universe_depth=-1, detect_errors=False):
        """Sample cells on a grid for YZ slice at x=constant."""
        self._rebuild_if_needed()
        y_min, y_max, z_min, z_max = bounds
        ny, nz = resolution
        return self._sys.find_cells_grid_x(x, y_min, y_max, z_min, z_max, ny, nz,
                                           universe_depth=universe_depth,
                                           detect_errors=detect_errors)

    def find_cells_grid(self, origin, normal, up, bounds,
                        resolution=(100, 100), universe_depth=-1, detect_errors=False):
        """Sample cells on a grid for arbitrary plane slice."""
        self._rebuild_if_needed()
        u_min, u_max, v_min, v_max = bounds
        nu, nv = resolution
        return self._sys.find_cells_grid(origin, normal, up, u_min, u_max, v_min, v_max, nu, nv,
                                         universe_depth=universe_depth,
                                         detect_errors=detect_errors)

    def find_label_positions(self, grid_result, min_pixels=100, by_material=False):
        """Find optimal label positions for regions in a grid."""
        from .slicing import find_label_positions
        return find_label_positions(self, grid_result, min_pixels, by_material)

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
        for cell in self._cells.values():
            if cell.material > 0 and cell.material not in self._materials:
                issues.append(f"Cell {cell.id} uses undefined material {cell.material}")

        # Check for overlaps
        overlaps = self.find_overlaps()
        for c1, c2 in overlaps:
            issues.append(f"Cells {c1.id} and {c2.id} overlap")

        # Check C library validation
        self._rebuild_if_needed()
        n_issues = self._sys.validate()
        if n_issues > 0:
            issues.append(f"Internal validation found {n_issues} issues")

        return issues

    # =========================================================================
    # Utilities
    # =========================================================================

    def set_verbose(self, enabled: bool) -> None:
        """Enable/disable verbose output."""
        self._ensure_sys()
        self._sys.set_verbose(enabled)

    def print_summary(self) -> None:
        """Print model summary."""
        print(str(self))

    def __str__(self) -> str:
        """Human-readable string representation."""
        n_cells = len(self._cells) if self._cells else (self._sys.cell_count if self._sys else 0)
        n_surfaces = len(self._surfaces) if self._surfaces else (self._sys.surface_count if self._sys else 0)
        n_universes = len(self._universes) if self._universes else 0
        return f"Model: {n_cells} cells, {n_surfaces} surfaces, {n_universes} universes"

    def __repr__(self) -> str:
        n_cells = len(self._cells) if self._cells else (self._sys.cell_count if self._sys else 0)
        n_surfaces = len(self._surfaces) if self._surfaces else (self._sys.surface_count if self._sys else 0)
        return f"Model(title='{self.title}', cells={n_cells}, surfaces={n_surfaces})"

    # =========================================================================
    # Spatial Indexing
    # =========================================================================

    def build_spatial_index(self) -> 'Model':
        """Build spatial index for fast queries.

        Builds a BVH over all cell instances (flattened bboxes, not geometry).
        This enables fast spatial queries and slice rendering without memory
        explosion from flattening large models with many FILL references.

        Call this once after loading a geometry with FILLs. The spatial index
        is used automatically by slice rendering functions.

        Returns:
            self for chaining

        Example:
            model = ath.load("big_model.inp")
            model.build_spatial_index()  # Fast slicing now works
            model.plot(z=0, bounds=(-10, 10, -10, 10))
        """
        self._rebuild_if_needed()
        self._sys.build_spatial_index()
        return self

    @property
    def spatial_index_instance_count(self) -> int:
        """Number of cell instances in spatial index (for diagnostics)."""
        if self._sys is None:
            return 0
        try:
            return self._sys.spatial_index_instance_count
        except AttributeError:
            return 0

    def flatten_universe(self, universe_id: int) -> None:
        """Flatten a universe by expanding all FILLs.

        Args:
            universe_id: Universe ID to flatten
        """
        self._rebuild_if_needed()
        self._sys.flatten_universe(universe_id)

    def simplify(self) -> Dict[str, Any]:
        """Run full CSG simplification on all cells.

        Applies complement elimination, double-negation removal,
        idempotent/absorption reductions, subtree deduplication, and more.

        Returns:
            Dict with simplification statistics (nodes_before, nodes_after, etc.)

        Example:
            stats = model.simplify()
            print(f"Reduced {stats['nodes_before']} -> {stats['nodes_after']} nodes")
        """
        self._rebuild_if_needed()
        return self._sys.flatten_all_cells()

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
        self._rebuild_if_needed()
        return self._sys.get_cells_by_material(material_id)

    def get_cells_by_universe(self, universe_id: int) -> List[int]:
        """Get indices of all cells in the given universe.

        Args:
            universe_id: Universe ID to filter by

        Returns:
            List of cell indices
        """
        self._rebuild_if_needed()
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
        self._rebuild_if_needed()
        return self._sys.get_cells_filling_universe(universe_id)

    def get_cells_in_bbox(self, bounds: Tuple[float, float, float, float, float, float]) -> List[int]:
        """Get indices of cells whose bounding boxes intersect the given region.

        Args:
            bounds: (x_min, x_max, y_min, y_max, z_min, z_max)

        Returns:
            List of cell indices
        """
        self._rebuild_if_needed()
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
        self._rebuild_if_needed()
        new_sys = self._sys.extract_universe(universe_id)
        new_model = Model.__new__(Model)
        new_model.title = f"Universe {universe_id}"
        new_model._sys = new_sys
        new_model._cells = {}
        new_model._surfaces = {}
        new_model._materials = {}
        new_model._universes = {}
        new_model._next_cell_id = 1
        new_model._next_surface_id = 1
        new_model._dirty = False
        return new_model

    def extract_region(self, bounds: Tuple[float, float, float, float, float, float]) -> 'Model':
        """Extract cells in a bounding box region into a new Model.

        Args:
            bounds: (x_min, x_max, y_min, y_max, z_min, z_max)

        Returns:
            New Model with the extracted cells
        """
        self._rebuild_if_needed()
        x_min, x_max, y_min, y_max, z_min, z_max = bounds
        new_sys = self._sys.extract_region(x_min, x_max, y_min, y_max, z_min, z_max)
        new_model = Model.__new__(Model)
        new_model.title = f"Region [{x_min:.1f},{x_max:.1f}]x[{y_min:.1f},{y_max:.1f}]x[{z_min:.1f},{z_max:.1f}]"
        new_model._sys = new_sys
        new_model._cells = {}
        new_model._surfaces = {}
        new_model._materials = {}
        new_model._universes = {}
        new_model._next_cell_id = 1
        new_model._next_surface_id = 1
        new_model._dirty = False
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
        self._rebuild_if_needed()
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
        self._rebuild_if_needed()

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
            format: Export format ('mcnp' or 'openmc'). If None, detected from extension.

        Example:
            model.save("output.inp")              # MCNP (default)
            model.save("geometry.xml")            # OpenMC (from extension)
            model.save("output", format="mcnp")   # Explicit format
        """
        path = Path(filename)

        # Determine format
        if format is None:
            ext = path.suffix.lower()
            if ext in ('.xml',):
                format = 'openmc'
            else:
                format = 'mcnp'

        if format == 'openmc':
            self.export_openmc(filename)
        else:
            self.export_mcnp(filename)

    # =========================================================================
    # Design Document API - Plotting (delegates to plotting module)
    # =========================================================================

    def plot(self, z=None, y=None, x=None, **kwargs):
        """Plot geometry slice using grid-based rendering."""
        from .plotting import plot as _plot
        return _plot(self, z=z, y=y, x=x, **kwargs)

    # =========================================================================
    # Volume Estimation
    # =========================================================================

    def compute_bounding_sphere(self, tolerance=1.0):
        """Compute a tight bounding sphere for the entire model."""
        self._rebuild_if_needed()
        return self._sys.compute_bounding_sphere(tolerance)

    def estimate_cell_volumes(self, n_rays=100000, center=None, radius=None):
        """Estimate cell volumes using random ray tracing."""
        self._rebuild_if_needed()
        if center is None or radius is None:
            cx, cy, cz, r = self.compute_bounding_sphere()
            if center is None:
                center = (cx, cy, cz)
            if radius is None:
                radius = r
        ox, oy, oz = center
        return self._sys.estimate_cell_volumes(ox, oy, oz, radius, n_rays)

    def estimate_instance_volumes(self, n_rays=100000):
        """Estimate volumes per cell instance (spatial-index aware)."""
        self._rebuild_if_needed()
        return self._sys.estimate_instance_volumes(n_rays)

    def remove_cells_by_volume(self, volumes, threshold):
        """Remove cells whose estimated volume is below a threshold."""
        self._rebuild_if_needed()
        return self._sys.remove_cells_by_volume(volumes, threshold)

    # =========================================================================
    # Geometry Transforms
    # =========================================================================

    def renumber_cells(self, start_id=1):
        """Renumber all cells with consecutive IDs."""
        self._rebuild_if_needed()
        return self._sys.renumber_cells(start_id)

    def renumber_surfaces(self, start_id=1):
        """Renumber all surfaces with consecutive IDs."""
        self._rebuild_if_needed()
        return self._sys.renumber_surfaces(start_id)

    def offset_cell_ids(self, offset):
        """Add an offset to all cell IDs."""
        self._rebuild_if_needed()
        self._sys.offset_cell_ids(offset)

    def offset_surface_ids(self, offset):
        """Add an offset to all surface IDs."""
        self._rebuild_if_needed()
        self._sys.offset_surface_ids(offset)

    def offset_material_ids(self, offset):
        """Add an offset to all material IDs."""
        self._rebuild_if_needed()
        self._sys.offset_material_ids(offset)

    def split_union_cells(self):
        """Split cells with top-level unions into multiple simpler cells."""
        self._rebuild_if_needed()
        return self._sys.split_union_cells()

    def expand_macrobodies(self):
        """Expand all macrobodies to primitive surfaces."""
        self._rebuild_if_needed()
        return self._sys.expand_macrobodies()

    def tighten_bboxes(self, tolerance=1.0):
        """Tighten all cell bounding boxes via interval arithmetic."""
        self._rebuild_if_needed()
        return self._sys.tighten_all_bboxes(tolerance)

    def tighten_cell_bbox(self, cell_id, tolerance=1.0):
        """Tighten a single cell's bounding box."""
        self._rebuild_if_needed()
        idx = self._sys.cell_find(cell_id)
        if idx is None:
            raise KeyError(f"Cell {cell_id} not found")
        return self._sys.tighten_cell_bbox(idx, tolerance)

    def tighten_bbox_numerical(self, cell_id: int) -> None:
        """Tighten a cell's bounding box using numerical sampling.

        Fallback for cells where interval arithmetic gives loose bounds.

        Args:
            cell_id: Cell ID (MCNP cell number)

        Raises:
            KeyError: If cell not found
            RuntimeError: On failure
        """
        self._rebuild_if_needed()
        idx = self._sys.cell_find(cell_id)
        if idx is None:
            raise KeyError(f"Cell {cell_id} not found")
        self._sys.tighten_cell_bbox_numerical(idx)

    # =========================================================================
    # Configuration
    # =========================================================================

    @property
    def config(self) -> Dict[str, Any]:
        """Get current system configuration as a dictionary.

        Returns:
            Dict with configuration keys like 'abs_tol', 'rel_tol',
            'log_level', 'dedup', 'export_materials', etc.
        """
        self._ensure_sys()
        return self._sys.get_config()

    @config.setter
    def config(self, settings: Dict[str, Any]) -> None:
        """Update system configuration from a dictionary.

        Only provided keys are changed; others keep their current values.

        Args:
            settings: Dict of configuration keys to update.

        Example:
            model.config = {'log_level': 3, 'abs_tol': 1e-8}
        """
        self._ensure_sys()
        self._sys.set_config(settings)

    # =========================================================================
    # Cell Fill Mutation
    # =========================================================================

    def set_fill(self, cell_id, fill_universe, transform=0):
        """Set the fill universe for a cell."""
        self._rebuild_if_needed()
        idx = self._sys.cell_find(cell_id)
        if idx is None:
            raise KeyError(f"Cell {cell_id} not found")
        self._sys.set_fill(idx, fill_universe, transform)
        # Keep Python Cell in sync so fills survive rebuild
        py_cell = self._cells.get(cell_id)
        if py_cell is not None:
            py_cell.fill = fill_universe if fill_universe >= 0 else None
            py_cell.fill_transform = transform

    # =========================================================================
    # Comprehensive Overlap Check (delegates to slicing module)
    # =========================================================================

    def check_grid_overlaps(self, grid_result, universe_depth=-1):
        """Run comprehensive overlap detection on a previously computed grid."""
        from .slicing import check_grid_overlaps
        return check_grid_overlaps(self, grid_result, universe_depth)

    def find_surface_label_positions(self, grid_result, margin=20):
        """Find optimal label positions for surfaces on a slice plane."""
        from .slicing import find_surface_label_positions
        return find_surface_label_positions(self, grid_result, margin)

    def plot_views(self, bounds=None, save=None, **kwargs):
        """Plot three orthogonal views (XY, XZ, YZ)."""
        from .plotting import plot_views as _plot_views
        return _plot_views(self, bounds=bounds, save=save, **kwargs)
