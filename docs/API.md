<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO

SPDX-License-Identifier: MPL-2.0
-->

# API Reference

Every public class, method, and function in AleaTHOR, grouped by what you're trying to do.

**Conventions**:
- All geometry operations require the C extension (`_alea`). If not installed, `RuntimeError` is raised.
- Cell IDs are MCNP-numbered (positive integers). Cell indices are 0-based internal positions.
- Methods that query the geometry automatically rebuild the internal C system if changes are pending.

---

## Loading Models

### ath.load

```python
ath.load(filename: str) -> Model
```

Load a geometry file, auto-detecting format from extension. `.xml` files are read as OpenMC. Everything else is read as MCNP (`.inp`, `.i`, `.mcnp`, or no extension).

### ath.read_mcnp

```python
ath.read_mcnp(filename: str) -> Model
```

Parse an MCNP input file. Handles cell cards, surface cards, data cards (materials, transforms), `LIKE BUT`, cell complements (`#cell`), macrobodies, universe fills, and lattices.

### ath.read_mcnp_string

```python
ath.read_mcnp_string(content: str) -> Model
```

Parse MCNP input from a string instead of a file.

### ath.read_openmc

```python
ath.read_openmc(filename: str) -> Model
```

Parse an OpenMC XML geometry file.

---

## Model

### Constructor

```python
Model(title: str = "AleaTHOR Model")
```

Create an empty model. Cells, materials, and surfaces are added programmatically.

### Properties

```python
model.title         # str: Model title
model.cells         # CellCollection: all cells, with filtering support
model.surfaces      # Dict[int, Surface]: all surfaces by ID
model.materials     # List[Material]: all materials
model.universes     # List[Universe]: all universes
model.config        # Dict[str, Any]: system configuration (getter/setter)
model.spatial_index_instance_count  # int: cell instances in spatial index
```

---

## Cell Management

### model.add_cell

```python
model.add_cell(
    region: Region,
    cell_id: int = None,      # Auto-assigned if None
    material: int = 0,         # 0 = void
    density: float = 0.0,
    name: str = None,
    universe: int = 0,
    fill: int = None,          # Universe to fill with
    importance: float = 1.0,
) -> Cell
```

Add a cell to the model. Returns the created `Cell` dataclass. Raises `ValueError` if `cell_id` already exists.

### model.get_cell

```python
model.get_cell(cell_id: int) -> Cell
model[cell_id] -> Cell    # Subscript notation
```

Get cell by MCNP cell ID. Returns a `Cell`. Raises `KeyError` if not found.

### model.update_cell

```python
model.update_cell(
    cell_id: int,
    material: int = None,     # None = keep current
    density: float = None,
    importance: float = None,
) -> Cell
```

Update cell properties. Only provided fields are changed. Marks the model dirty. Returns a fresh `Cell`. Raises `KeyError` if not found.

### model.remove_cell

```python
model.remove_cell(cell_id: int) -> None
```

Remove cell by ID. No-op if not found.

---

## Material Management

### model.add_material

```python
model.add_material(
    material_id: int,
    name: str = None,
    density: float = 0.0,
    composition: Dict[str, float] = None,
) -> Material
```

Add a material. Raises `ValueError` if ID already exists.

### model.get_material

```python
model.get_material(material_id: int) -> Material
```

Get material by ID. Raises `KeyError` if not found.

### model.create_mixture

```python
model.create_mixture(
    material_ids: List[int],
    fractions: List[float],
    new_id: int = 0,          # 0 = auto-assign
    name: str = None,
) -> int
```

Create a weighted mixture of existing materials. Fractions are normalized automatically. Returns the assigned material ID.

---

## Universe Management

### model.add_universe

```python
model.add_universe(universe_id: int, name: str = None) -> Universe
```

### model.get_universe

```python
model.get_universe(universe_id: int) -> Universe
```

---

## Point Queries

### model.cell_at

```python
model.cell_at(x: float, y: float, z: float) -> Optional[Cell]
```

Find the innermost cell containing the point. Traverses the full universe hierarchy. Returns `None` if the point is in void or undefined.

### model.cells_at

```python
model.cells_at(x: float, y: float, z: float) -> List[Cell]
```

Find **all** cells containing the point across all hierarchy depths, ordered by depth (0 = outermost). Each `Cell` has its `depth` property set.

### model.contains_point

```python
model.contains_point(x: float, y: float, z: float) -> bool
```

Check if point is inside any cell (not void).

### model.find_overlaps

```python
model.find_overlaps(max_pairs: int = 100) -> List[Tuple[Cell, Cell]]
```

Find overlapping cell pairs by statistical sampling.

---

## Ray Tracing

### model.trace

```python
model.trace(
    origin: Tuple[float, float, float] = None,
    direction: Tuple[float, float, float] = None,
    start: Tuple[float, float, float] = None,
    end: Tuple[float, float, float] = None,
    max_distance: float = 0,       # 0 = infinite
    cell_aware: bool = False,      # Faster for large models
) -> TraceResult
```

Trace a ray through the geometry. Call in one of two ways:

- **Direction mode**: `trace(origin=(0,0,0), direction=(1,0,0))`
- **Point-to-point mode**: `trace(start=(0,0,0), end=(100,0,0))`

The direction is normalized internally. In point-to-point mode, `max_distance` is set automatically.

Cell-aware tracing tests only surfaces belonging to the current cell at each step, which is faster for models with many surfaces.

---

## Cell Filtering

C-backed for performance. Return lists of cell **indices** (not IDs).

```python
model.get_cells_by_material(material_id: int) -> List[int]
model.get_cells_by_universe(universe_id: int) -> List[int]
model.get_cells_filling_universe(universe_id: int) -> List[int]
model.get_cells_in_bbox(bounds: Tuple[float, float, float, float, float, float]) -> List[int]
```

For most use cases, prefer `model.cells.by_material()` etc. which return `CellCollection` objects with `Cell` iteration.

---

## Extraction

### model.extract_universe

```python
model.extract_universe(universe_id: int) -> Model
```

Extract a universe and all universes it references into a new `Model`.

### model.extract_region

```python
model.extract_region(bounds: Tuple[float, ...]) -> Model
```

Extract cells whose bounding boxes intersect the given `(xmin, xmax, ymin, ymax, zmin, zmax)` region into a new `Model`.

---

## Slice API

All slice functions take a `bounds` tuple of `(min1, max1, min2, max2)` defining the viewport in the slice plane.

### Analytical Curves

```python
model.get_slice_curves_z(z: float, bounds) -> dict
model.get_slice_curves_y(y: float, bounds) -> dict
model.get_slice_curves_x(x: float, bounds) -> dict
model.get_slice_curves(origin, normal, up, bounds) -> dict
```

Return a dict with `'curves'` (list of curve dicts, each with `'type'`, `'surface_id'`, and geometry-specific fields) and bounds info.

### Grid Queries

```python
model.find_cells_grid_z(z, bounds, resolution=(100,100),
                         universe_depth=-1, detect_errors=False) -> dict
model.find_cells_grid_y(y, bounds, ...) -> dict
model.find_cells_grid_x(x, bounds, ...) -> dict
model.find_cells_grid(origin, normal, up, bounds, ...) -> dict
```

Sample cell/material IDs on a 2D pixel grid. Returns a dict with:

| Key | Type | Description |
|-----|------|-------------|
| `cell_ids` | `list[int]` | Cell ID at each pixel (-1 = void) |
| `material_ids` | `list[int]` | Material ID at each pixel (0 = void) |
| `errors` | `list[int]` | Error flags (only if `detect_errors=True`) |
| `nx`, `ny` (or `nu`, `nv`) | `int` | Grid dimensions |
| `x_min`, `x_max`, ... | `float` | Viewport bounds |

**`universe_depth`**: `-1` = innermost cell, `0` = root only, `N` = cells at depth N.

**`detect_errors`**: when True, boundary pixels are rechecked for overlaps. Error values: 0=ok, 1=overlap, 2=undefined.

### Label Positioning

```python
model.find_label_positions(
    grid_result: dict,
    min_pixels: int = 100,
    by_material: bool = False,
) -> List[dict]
```

Find optimal label positions for regions in a grid. Returns list of `{'id', 'px', 'py', 'pixel_count'}`. The `min_pixels` parameter filters out tiny regions.

```python
model.find_surface_label_positions(
    grid_result: dict,
    margin: int = 20,
) -> List[dict]
```

Find label positions for surfaces on a slice plane. Returns list of `{'id', 'px', 'py'}`.

### Overlap Checking

```python
model.check_grid_overlaps(
    grid_result: dict,
    universe_depth: int = -1,
) -> List[int]
```

Comprehensive overlap detection: re-queries every non-void pixel. Returns updated error list (0=ok, 1=overlap, 2=undefined). This is O(area) — use only when thorough validation is needed.

---

## Plotting

Requires `matplotlib` and `numpy`. Import check: `ath.HAS_PLOTTING`.

### model.plot

```python
model.plot(
    z=None, y=None, x=None,          # Axis-aligned slice (use one)
    origin=None, normal=None, up=None, # Arbitrary plane
    bounds=None,                       # (min1, max1, min2, max2)
    filled: bool = True,
    by_material: bool = False,
    show_colorbar: bool = False,
    show_contours: bool = True,
    detect_errors: bool = False,
    universe_depth: int = -1,
    save: str = None,                  # Save to file
    ax=None,                           # Existing matplotlib Axes
    overlay=None,                      # 2D array for tally overlay
    overlay_extent=None,
    overlay_cmap: str = 'jet',
    overlay_norm: str = None,          # 'log' for logarithmic
    overlay_alpha: float = 1.0,
    overlay_vmin=None, overlay_vmax=None,
    overlay_label: str = None,
    contour_by: str = None,            # 'cell', 'material', or None (auto)
) -> Axes
```

Plot geometry slice. If no axis is specified, defaults to `z=0`. Returns the matplotlib `Axes`.

When `by_material=True` or `overlay` is set, contours automatically follow material boundaries. Override with `contour_by='cell'`.

Also available as a module-level function:

```python
from aleathor.plotting import plot
ax = plot(model, z=0, bounds=(-10, 10, -10, 10))
```

### model.plot_views

```python
model.plot_views(
    bounds=None,    # (xmin, xmax, ymin, ymax, zmin, zmax) — 6-element
    save: str = None,
    **kwargs,       # Passed to plot()
) -> Figure
```

Plot three orthogonal views (XY at z=0, XZ at y=0, YZ at x=0). Returns the matplotlib `Figure`.

### Low-level plotting functions

```python
from aleathor.plotting import (
    plot_slice_curves,     # Vector plot of analytical curves
    plot_slice_filled,     # Grid-based filled plot with contours
    plot_curve,            # Plot a single curve
    plot_ray_path,         # 1D bar chart of ray trace segments
    filter_boundary_curves,            # Filter to true cell boundaries
    filter_boundary_curves_arbitrary,  # Same for arbitrary planes
    count_surfaces,        # Count unique surfaces in curves result
    count_cells,           # Count unique cells in grid result
    count_materials,       # Count unique materials in grid result
    get_slice_stats,       # All three counts at once
)
```

### plot_slice_curves

```python
plot_slice_curves(
    curves_result: dict,
    ax=None,
    color: str = 'black',
    linewidth: float = 0.5,
    title: str = None,
    xlabel: str = None, ylabel: str = None,
    xlim=None, ylim=None,
) -> Axes
```

Plot analytical surface boundaries as vector graphics.

### plot_slice_filled

```python
plot_slice_filled(
    grid_result: dict,
    ax=None,
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
    title: str = None,
    overlay=None,                    # 2D array for tally data
    overlay_cmap: str = 'jet',
    overlay_norm: str = None,
    overlay_alpha: float = 1.0,
    overlay_vmin=None, overlay_vmax=None,
    overlay_label: str = None,
    contour_by: str = None,
) -> Axes
```

Grid-based filled plot with contour lines. The recommended way to visualize geometry slices.

### plot_ray_path

```python
plot_ray_path(
    segments: list,
    ax=None,
    show_materials: bool = True,
    t_max: float = None,
) -> Axes
```

Plot ray path as a 1D bar chart showing material transitions.

---

## Export

### model.export_mcnp

```python
model.export_mcnp(filename: str) -> None
```

### model.export_openmc

```python
model.export_openmc(filename: str) -> None
```

### model.save

```python
model.save(filename: str, format: str = None) -> None
```

Auto-detects format from extension. `.xml` = OpenMC, everything else = MCNP. Override with `format='mcnp'` or `format='openmc'`.

### model.to_mcnp_string

```python
model.to_mcnp_string() -> str
```

Generate MCNP input file contents as a string.

### File-level I/O functions

```python
ath.write_mcnp(model, filename, deduplicate=True) -> None
ath.write_openmc(model, filename) -> None
```

---

## Mesh Export

### model.export_mesh

```python
model.export_mesh(
    filename: str,
    nx: int = 10, ny: int = 10, nz: int = 10,
    bounds: Tuple[float, ...] = None,   # (xmin, xmax, ymin, ymax, zmin, zmax), None = auto
    format: str = "gmsh",                # "gmsh" or "vtk"
    void_material_id: int = 0,
    auto_pad: float = 0.01,
) -> None
```

Export the geometry as a structured hexahedral mesh. Each element is assigned the material ID at its center. Bounds are auto-detected from cell bounding boxes if not given.

Supported formats:
- `"gmsh"`: Gmsh `.msh` v2.2 ASCII — can be opened in Gmsh for visualization
- `"vtk"`: VTK legacy `.vtk` ASCII — can be opened in ParaView

### model.sample_mesh

```python
model.sample_mesh(
    nx: int = 10, ny: int = 10, nz: int = 10,
    bounds: Tuple[float, ...] = None,
    void_material_id: int = 0,
    auto_pad: float = 0.01,
) -> dict
```

Sample the geometry on a structured mesh without writing to file. Returns a dict:

| Key | Type | Description |
|-----|------|-------------|
| `material_ids` | `list[int]` | Material at each element (nx*ny*nz, Z-major order) |
| `cell_ids` | `list[int]` | Cell ID at each element |
| `x_nodes` | `list[float]` | X node positions (nx+1 values) |
| `y_nodes` | `list[float]` | Y node positions (ny+1 values) |
| `z_nodes` | `list[float]` | Z node positions (nz+1 values) |
| `nx`, `ny`, `nz` | `int` | Grid dimensions |

---

## Geometry Transforms

### Renumbering

```python
model.renumber_cells(start_id: int = 1) -> int
model.renumber_surfaces(start_id: int = 1) -> int
```

Reassign IDs starting from `start_id`, preserving order. Returns the count.

### Offsetting

```python
model.offset_cell_ids(offset: int) -> None
model.offset_surface_ids(offset: int) -> None
model.offset_material_ids(offset: int) -> None
```

Add a fixed offset to all IDs. Useful when merging models.

### Splitting and expanding

```python
model.split_union_cells() -> int
```

Split cells with top-level unions into multiple simpler cells. Returns number of new cells created.

```python
model.expand_macrobodies() -> int
```

Expand all macrobodies to primitive surfaces. Returns number expanded.

### Bounding boxes

```python
model.tighten_bboxes(tolerance: float = 1.0) -> int
```

Tighten all cell bounding boxes via interval arithmetic. Returns number tightened.

```python
model.tighten_cell_bbox(cell_id: int, tolerance: float = 1.0)
    -> Tuple[float, float, float, float, float, float]
```

Tighten a single cell's bounding box. Returns `(xmin, xmax, ymin, ymax, zmin, zmax)`. Raises `KeyError` if cell not found.

```python
model.tighten_bbox_numerical(cell_id: int) -> None
```

Tighten a cell's bounding box using numerical sampling. This is a fallback for cells where interval arithmetic gives loose bounds (e.g., cells with complex boolean expressions). Raises `KeyError` if cell not found, `RuntimeError` on failure.

### Fill mutation

```python
model.set_fill(cell_id: int, fill_universe: int, transform: int = 0) -> None
```

Set the fill universe for a cell. Operates on the C system directly. Raises `KeyError` if cell not found.

---

## Volume Estimation

### model.compute_bounding_sphere

```python
model.compute_bounding_sphere(tolerance: float = 1.0)
    -> Tuple[float, float, float, float]
```

Compute a tight bounding sphere. Returns `(cx, cy, cz, radius)`.

### model.estimate_cell_volumes

```python
model.estimate_cell_volumes(
    n_rays: int = 100000,
    center: Tuple[float, float, float] = None,  # Auto-computed if None
    radius: float = None,                         # Auto-computed if None
) -> dict
```

Estimate cell volumes using random ray tracing (Cauchy-Crofton method). Returns `{'volumes': list, 'rel_errors': list}` indexed by cell index.

If `center` or `radius` are not given, a bounding sphere is computed automatically.

### model.estimate_instance_volumes

```python
model.estimate_instance_volumes(n_rays: int = 100000) -> dict
```

Estimate volumes per cell instance (spatial-index aware). Requires `build_spatial_index()` first.

### model.remove_cells_by_volume

```python
model.remove_cells_by_volume(volumes: list, threshold: float) -> int
```

Remove cells whose estimated volume is below threshold. Returns number removed.

---

## Spatial Indexing

### model.build_spatial_index

```python
model.build_spatial_index() -> Model
```

Build a KD-tree over cell instances for fast queries. Returns `self` for chaining. **Call this once after loading large models with FILLs** before doing slice rendering.

### model.flatten_universe

```python
model.flatten_universe(universe_id: int) -> None
```

Expand all fills in a universe. After flattening, all cells are in the specified universe with no fills.

### model.simplify

```python
model.simplify() -> dict
```

Run full CSG simplification on all cells. Applies complement elimination, double-negation removal, idempotent/absorption reductions, subtree deduplication, and more. Returns a dict with statistics:

| Key | Type | Meaning |
|-----|------|---------|
| `nodes_before` | `int` | Total CSG nodes before simplification |
| `nodes_after` | `int` | Total CSG nodes after simplification |
| `complements_eliminated` | `int` | Complement nodes removed |
| `double_negations` | `int` | Double negations eliminated |
| `idempotent_reductions` | `int` | Idempotent reductions (A & A = A) |
| `absorption_reductions` | `int` | Absorption reductions (A & universe = A) |
| `subtrees_deduplicated` | `int` | Duplicate subtrees merged |
| `cell_complements_expanded` | `int` | Cell complements expanded inline |
| `contradictions_found` | `int` | Contradictions reduced to empty |
| `tautologies_found` | `int` | Tautologies reduced to universe |
| `empty_cells_removed` | `int` | Empty cells removed |
| `union_branches_absorbed` | `int` | Union branches absorbed |
| `union_common_factors` | `int` | Common factors extracted from unions |
| `union_branches_subsumed` | `int` | Union branches subsumed |

---

## Validation

### model.validate

```python
model.validate() -> List[str]
```

Check for common issues: undefined materials in cells, overlapping cells, internal validation. Returns list of warning/error messages (empty = valid).

---

## Configuration

### model.config

```python
# Read
cfg = model.config    # Returns dict

# Write (only provided keys are changed)
model.config = {'log_level': 3, 'abs_tol': 1e-8}
```

### model.set_verbose

```python
model.set_verbose(enabled: bool) -> None
```

Enable/disable verbose output from the C library.

---

## Utilities

```python
model.print_summary() -> None
str(model)        # "Model: 45 cells, 89 surfaces, 3 universes"
repr(model)       # "Model(title='...', cells=45, surfaces=89)"
```

---

## Cell

Mutable view of a cell in the C system. Returned by `model.get_cell()`, `model.cell_at()`, and iteration over `model.cells`.

### Properties

```python
cell.id: int              # Cell ID (MCNP cell number)
cell.material: int        # Material ID (0 for void)         [read-write]
cell.density: float       # Material density                 [read-write]
cell.density_unit: str    # "g/cm3" or "atoms/b-cm"          [read-write]
cell.universe: int        # Universe this cell belongs to
cell.fill: Optional[int]  # Fill universe (None if not filled) [read-write]
cell.name: Optional[str]  # Cell name                        [read-write]
cell.bounds: Tuple[float, float, float, float, float, float]  # Bounding box
cell.is_void: bool        # True if material == 0
cell.is_filled: bool      # True if fill is set
cell.importance: float    # Particle importance               [read-write]
cell.depth: int           # Hierarchy depth (0 = root)
cell.region               # Region or _ImportedRegion
```

### Lattice properties

Only populated for cells that define a lattice (`lat_type != 0`). For non-lattice cells, these return `None` or `0`/`False`.

```python
cell.is_lattice: bool                     # True if this cell defines a lattice
cell.lattice_type: int                    # 0=none, 1=rectangular, 2=hexagonal
cell.lattice_pitch: Optional[Tuple[float, float, float]]    # Element pitch (x, y, z)
cell.lattice_fill: Optional[List[int]]    # Universe IDs filling lattice elements
cell.lattice_lower_left: Optional[Tuple[float, float, float]]  # Lower-left corner
cell.lattice_dims: Optional[Tuple[int, int, int, int, int, int]]  # (imin,imax,jmin,jmax,kmin,kmax)
```

### Methods

```python
cell.contains(x: float, y: float, z: float) -> bool
cell.fill_with(universe, transform: int = 0) -> None
cell.unfill() -> None
```

---

## CellCollection

Collection of cells with filtering.

### Basic operations

```python
len(collection)                    # Number of cells
collection[0]                      # First Cell
collection[cell_id]                # Cell by ID
for cell in collection: ...        # Iterate Cells
cell_view in collection            # Membership test
```

### Filtering

```python
collection.by_material(material_id: int) -> CellCollection
collection.by_universe(universe_id: int) -> CellCollection
collection.by_fill(fill_universe: int = None) -> CellCollection
collection.filter(predicate: Callable[[Cell], bool]) -> CellCollection
```

Filters can be chained: `model.cells.by_material(1).by_universe(0)`.

### Utility

```python
collection.get(cell_id: int) -> Optional[Cell]  # Get by ID, None if missing
collection.ids() -> List[int]
collection.materials() -> Set[int]
collection.universes() -> Set[int]
collection.to_list() -> List[Cell]
```

---

## TraceResult

Result of a ray trace through geometry.

### Iteration

```python
len(result)             # Number of segments
result[0]               # First TraceSegment
for seg in result: ...  # Iterate segments
```

### Analysis

```python
result.path_length(material: int = None) -> float
```

Total path length through a specific material, or all materials if `material=None`.

```python
result.cells_hit() -> List[Optional[Cell]]
result.materials_hit() -> Set[int]
```

---

## TraceSegment

A segment of a ray trace.

### Properties

```python
seg.cell: Optional[Cell]  # Cell (None for void)
seg.t_enter: float             # Entry distance along ray
seg.t_exit: float              # Exit distance along ray
seg.material: int              # Material ID (0 for void)
seg.density: float             # Material density
seg.length: float              # t_exit - t_enter
seg.is_void: bool              # True if cell is None
```

---

## Surface Classes

All surfaces inherit from `Surface` and support:

```python
-surface    # Negative halfspace (inside)
+surface    # Positive halfspace (outside)
surface.interior()  # Same as -surface
surface.exterior()  # Same as +surface
surface.evaluate(point)  # Signed distance function
```

### Primitives

```python
Plane(a, b, c, d, name=None, boundary='transmissive')
XPlane(x0, ...)
YPlane(y0, ...)
ZPlane(z0, ...)

Sphere(x0, y0, z0, radius, ...)
CylinderX(y0, z0, radius, ...)
CylinderY(x0, z0, radius, ...)
CylinderZ(x0, y0, radius, ...)
ConeX(x0, y0, z0, t_sq, sheet=0, ...)
ConeY(x0, y0, z0, t_sq, sheet=0, ...)
ConeZ(x0, y0, z0, t_sq, sheet=0, ...)
TorusX(x0, y0, z0, major_radius, minor_radius, ...)
TorusY(x0, y0, z0, major_radius, minor_radius, ...)
TorusZ(x0, y0, z0, major_radius, minor_radius, ...)
Quadric(A, B, C, D, E, F, G, H, J, K, ...)
```

### Macrobodies

```python
Box(xmin, xmax, ymin, ymax, zmin, zmax, ...)
RCC(base_x, base_y, base_z, height_x, height_y, height_z, radius, ...)
TRC(base_x, base_y, base_z, height_x, height_y, height_z, base_radius, top_radius, ...)
ELL(v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, major_axis_len, ...)
REC(base_x, ..., height_x, ..., axis1_x, ..., axis2_x, ..., ...)
WED(vertex_x, ..., v1_x, ..., v2_x, ..., v3_x, ..., ...)
RHP(base_x, ..., height_x, ..., r1_x, ..., r2_x, ..., r3_x, ..., ...)
GeneralBox(corner_x, ..., v1_x, ..., v2_x, ..., v3_x, ..., ...)
```

---

## Region Classes

### Halfspace

```python
Halfspace(surface, positive)   # Created by -surface or +surface
```

Properties: `surface`, `positive`.

### Intersection

```python
region_a & region_b    # Creates Intersection([a, b])
```

Property: `regions` (list of child regions). Nested intersections are flattened.

### Union

```python
region_a | region_b    # Creates Union([a, b])
```

Property: `regions` (list). Nested unions are flattened.

### Complement

```python
~region    # Creates Complement(region)
```

Double negation is eliminated: `~~region` returns the original region.

### All Regions

```python
(x, y, z) in region         # Point containment test
region.get_surfaces()        # Set of all referenced surfaces
```

---

## Dataclass: Material

```python
@dataclass
class Material:
    id: int                           # Must be positive
    name: Optional[str] = None
    density: float = 0.0
    composition: Dict[str, float] = {}
```

---

## Dataclass: _CellData

Internal representation used when building geometry programmatically. Users interact with the `Cell` class (in `collections.py`) for the public API.

```python
@dataclass
class _CellData:
    id: int                           # Must be positive
    region: Region
    material: int = 0                 # 0 = void
    density: float = 0.0
    name: Optional[str] = None
    universe: int = 0
    fill: Optional[int] = None
    fill_transform: Optional[int] = None
    importance: float = 1.0

    is_void: bool      # Property: material == 0
    is_filled: bool    # Property: fill is not None
```

---

## Dataclass: Universe

```python
@dataclass
class Universe:
    id: int
    name: Optional[str] = None
    cells: List[Cell] = []

    def add_cell(cell: Cell) -> None
```

---

## Logging

```python
ath.enable_logging()                 # Set INFO level
ath.disable_logging()                # Silence
ath.set_log_level(ath.LOG_DEBUG)     # Specific level
ath.get_log_level() -> int
```

### Log level constants

| Constant | Value | Meaning |
|----------|-------|---------|
| `ath.LOG_NONE` | 0 | No output |
| `ath.LOG_ERROR` | 1 | Errors only |
| `ath.LOG_WARN` | 2 | Warnings and above |
| `ath.LOG_INFO` | 3 | Informational |
| `ath.LOG_DEBUG` | 4 | Debug detail |
| `ath.LOG_TRACE` | 5 | Full trace |

---

## Constants

### Grid/Cell identification

| Constant | Value | Meaning |
|----------|-------|---------|
| `ath.CELL_VOID` | -1 | Void pixel in grid |
| `ath.CELL_UNDEFINED` | -2 | Undefined region |
| `ath.CELL_OVERLAP` | -3 | Overlapping cells |

### Grid error flags

| Constant | Value | Meaning |
|----------|-------|---------|
| `ath.GRID_ERROR_OK` | 0 | No error |
| `ath.GRID_ERROR_OVERLAP` | 1 | Overlap detected |
| `ath.GRID_ERROR_UNDEFINED` | 2 | Undefined region |

---

## Module: slicing

Helpers used internally by Model methods. Available for advanced use.

```python
from aleathor.slicing import (
    find_label_positions,
    find_surface_label_positions,
    check_grid_overlaps,
    _extract_slice_params,       # DRY helper for grid dimension dispatch
)
```

### _extract_slice_params

```python
_extract_slice_params(grid_result: dict)
    -> Tuple[nu, nv, origin, normal, up, u_min, u_max, v_min, v_max]
```

Extract unified slice parameters from any grid result dict (z, y, x, or arbitrary). Handles the dispatch logic for different grid dimension keys (`nx`/`ny`, `nx`/`nz`, `ny`/`nz`, `nu`/`nv`).
