<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO

SPDX-License-Identifier: MPL-2.0
-->

# Architecture

This document explains how the aleathor Python package works internally. Read this if you want to contribute, understand the design decisions, or debug issues at the boundary between Python and C.

## The Big Picture

aleathor Python is a thin layer over a C geometry engine ([libalea](https://github.com/giovanni-mariano/aleathor/tree/main/csrc/libalea)). The package structure:

```
src/aleathor/
    __init__.py        Public API surface, optional imports
    model.py           Model, Material, Universe, and C-backed model operations
    collections.py     Cell, CellCollection, TraceResult, TraceSegment
    surfaces.py        Surface classes (Sphere, Box, CylinderZ, etc.)
    geometry.py        Region classes (Halfspace, Intersection, Union, Complement)
    io.py              File I/O (read_mcnp, write_openmc, etc.)
    slicing.py         DRY helper for slice parameter extraction
    plotting.py        matplotlib visualization
    _binding/
        aleathor_binding.c   CPython extension entry point
        _bind_*.c            Binding implementation split by feature area
```

The C extension (`_alea`) is a CPython module that wraps `libalea` functions into Python-callable methods on a `System` class. The Python `Model` class owns a `System` instance (`model._sys`) and manages the translation between Python objects and C data.

## Model State

The `Model` class is C-backed. Geometry cells are pushed into `_sys` as soon as they are added, and queries operate on `_sys`. Python keeps lightweight side metadata for information that is not fully owned by the C system.

```
Model                         C side
├── _sys: _alea.System   ───► cells, CSG nodes, primitives, materials,
│                            fills, universes, and acceleration structures
├── _regions: Dict[int, Region]      Python region or imported-region view
├── _surfaces: Dict[int, Surface]    Python-created surfaces by ID
├── _materials: Dict[int, Material]  Material wrappers by material ID
├── _names: Dict[int, str]           Python-side cell names
└── _importances: Dict[int, float]   Python-side importance metadata
```

### When you build geometry in Python

```python
model.add_cell(region=-sphere, material=1, density=10.5)
```

1. `region._to_csg(model)` creates or reuses C CSG nodes.
2. The material is ensured in the C material table.
3. `_sys.add_cell(...)` creates the C-backed cell immediately.
4. Python stores metadata such as name, importance, region, and surfaces.
5. A `Cell` view is returned.

### When you query the geometry

```python
cell = model.cell_at(0, 0, 0)
```

1. `_ensure_query_caches()` prepares universe, spatial, raycast, and slice caches if needed.
2. `_sys.find_cell(x, y, z)` performs the C query.
3. The returned cell index is wrapped in a `Cell` object.

### When you load from a file

```python
model = ath.read_mcnp("geometry.inp")
```

The C parser builds the system directly. `_populate_regions()` then attaches `_ImportedRegion` wrappers to each loaded cell ID so users can inspect or reuse imported CSG regions from Python.

This is important: loaded models don't duplicate geometry in Python. The `Cell` objects returned by queries go straight to the C system.

## Query Cache Lifecycle

```
cell_at() / trace() / find_cells_grid() / plot()
        │
        ▼
    _ensure_query_caches()
        │
        ├── _sys.build_universe_index()
        ├── _sys.build_spatial_index()
        ├── _sys.prepare_raycast()
        └── _sys.prepare_slice()
```

The model no longer rebuilds all C cells from a Python `_cells` dictionary on each edit. Cell edits are applied directly to `_sys`; query preparation only builds acceleration structures.

The cache `_halfspace_node_cache` maps `(surface_id, positive)` to C node IDs. This prevents recreating the same surface primitive when multiple cells reference it.

## _get_or_create_halfspace_node

This is the bridge between Python `Surface` objects and C primitives. It's a large isinstance chain (~150 lines) that dispatches on surface type:

```python
def _get_or_create_halfspace_node(self, surface, positive):
    cache_key = (surface.id, positive)
    if cache_key in self._halfspace_node_cache:
        return self._halfspace_node_cache[cache_key]

    if isinstance(surface, Sphere):
        _, pos_node, neg_node = self._sys.sphere_surface(...)
    elif isinstance(surface, CylinderZ):
        _, pos_node, neg_node = self._sys.cylinder_z_surface(...)
    # ... etc for all ~20 surface types
```

Each `_surface()` call creates the C primitive and returns both sense nodes. Both are cached. When a cell references `-S5`, the cache returns the `neg_node` directly without creating anything.

## Cell Views

`Cell` (in `collections.py`) is the public API type. It holds a reference to the `Model`, an integer index into the C system's cell array, and optionally the original cell ID so it can resolve itself again after changes. Property access (`cell.material`, `cell.bounds`) reads from the C system. Property setters (`cell.material = 2`, `cell.density = 5.0`, `cell.fill = 10`) update the C system directly.

```python
class Cell:
    def __init__(self, model, index, _depth=0):
        self._model = model
        self._index = index
        self._depth = _depth

    @property
    def material(self):
        return self._get_info()['material_id']
```

## CellCollection

`CellCollection` wraps a list of cell indices and provides filtering:

```python
class CellCollection:
    def __init__(self, model, indices=None):
        self._model = model
        self._indices = indices  # None = all cells
```

When `indices` is `None`, iteration creates `Cell` objects for all cells (`0` to `cell_count - 1`). Filtering methods return a new `CellCollection` with a subset of indices.

The `by_material()` and `by_universe()` filters call C functions (`get_cells_by_material`, etc.) for performance rather than iterating in Python.

## TraceResult and TraceSegment

Ray tracing is delegated entirely to C. The Python side:

1. Normalizes the direction vector
2. Calls `_sys.raycast()` or `_sys.raycast_cell_aware()`
3. Gets back a list of segment dicts from C
4. Wraps each in a `TraceSegment` (which lazily creates a `Cell` for the cell)

```python
class TraceSegment:
    @property
    def cell(self):
        if self._cell_idx < 0:
            return None  # void
        return Cell(self._model, self._cell_idx)
```

## The Slicing Module

`slicing.py` exists to eliminate code triplication. Three functions — `find_label_positions`, `find_surface_label_positions`, and `check_grid_overlaps` — all need to extract the same set of parameters (dimensions, origin, normal, up, bounds) from a grid result dict. The dict format differs depending on whether it came from a Z, Y, X, or arbitrary slice.

The `_extract_slice_params()` helper handles this dispatch once:

```python
def _extract_slice_params(grid_result):
    if 'nx' in grid_result and 'ny' in grid_result:
        # Z-slice: nx, ny -> nu, nv; z -> origin; normal=(0,0,1)
    elif 'nx' in grid_result and 'nz' in grid_result:
        # Y-slice: nx, nz -> nu, nv; y -> origin; normal=(0,1,0)
    elif 'ny' in grid_result and 'nz' in grid_result:
        # X-slice: ny, nz -> nu, nv; x -> origin; normal=(1,0,0)
    elif 'nu' in grid_result and 'nv' in grid_result:
        # Arbitrary plane: direct extraction
```

Model methods that delegate to `slicing.py` are: `find_label_positions`, `find_surface_label_positions`, `check_grid_overlaps`. These are the only methods that use the delegate pattern — justified because they contain non-trivial dispatch logic.

## The Plotting Module

`plotting.py` is separated from `model.py` to keep plotting code out of the core model object. The module is imported lazily when plotting is requested.

The key functions are:

- `plot()`: high-level orchestrator. Calls `find_cells_grid()` and the appropriate `get_slice_curves_*()`, then delegates to `plot_slice_filled()` or `plot_slice_curves()`.
- `plot_views()`: creates a 1x3 subplot figure with XY, XZ, YZ slices.
- `plot_slice_filled()`: the workhorse. Takes a grid result, creates an image with cell/material coloring, overlays contour lines, handles error highlighting, and supports tally overlay.

Model delegates to these via lazy imports:

```python
def plot(self, z=None, y=None, x=None, **kwargs):
    from .plotting import plot as _plot
    return _plot(self, z=z, y=y, x=x, **kwargs)
```

## _ImportedRegion

When geometry is loaded from a file, cells reference `_ImportedRegion` objects instead of Python-created `Region` trees. An `_ImportedRegion` is a thin wrapper around a C system and CSG node ID:

```python
class _ImportedRegion(Region):
    def __init__(self, sys, node_id):
        self._sys = sys
        self._node_id = node_id

    def _to_csg(self, model):
        return self._node_id
```

This is why loaded models are efficient: they don't duplicate the CSG tree in Python. The `_to_csg()` method returns the existing C node directly.

`_ImportedRegion` also supports lazy reconstruction of a Python `Region` tree (the `.tree` property) for inspection. This walks the C CSG tree and creates `Halfspace`, `Intersection`, `Union`, and `Complement` objects.

## Surface Deduplication

The C library automatically deduplicates surfaces during loading (canonicalizing coefficients, hashing, and matching within tolerance). This is transparent to the Python side.

For programmatically-built geometry, deduplication happens implicitly through the `_halfspace_node_cache`: if you create two `Sphere` objects with the same parameters but different Python IDs, they'll get different C surface entries. This is acceptable because programmatic models rarely have duplicate surfaces.

## Interrupt Support

The C library checks a global interrupt flag during long operations. The Python binding installs a SIGINT handler:

1. User presses Ctrl+C
2. Python's default SIGINT handler raises `KeyboardInterrupt`
3. But if we're inside a C call, Python can't raise immediately
4. The binding calls `alea_interrupt()` — sets the C flag
5. The C operation returns `ALEA_ERR_INTERRUPTED`
6. The binding calls `alea_clear_interrupt()` and raises `KeyboardInterrupt`

This allows clean termination of grid queries and ray traces that might take seconds on large models.

## Module Boundaries

```
__init__.py  ─── re-exports public API
    │
    ├── model.py ────── Model, Cell, Material, Universe
    │     │
    │     ├── geometry.py ──── Region, Halfspace, Intersection, Union, Complement
    │     ├── surfaces.py ──── Surface, Sphere, Box, CylinderZ, ...
    │     ├── collections.py ── Cell, CellCollection, TraceResult, TraceSegment
    │     └── slicing.py ───── _extract_slice_params, find_label_positions, ...
    │
    ├── io.py ──────── read_mcnp, write_openmc, _ImportedRegion
    │
    └── plotting.py ── plot, plot_views, plot_slice_filled, ...
```

Dependencies flow downward. `plotting.py` depends on `model.py` (via `Model` type hints) but only through `TYPE_CHECKING` — no runtime import cycle. `slicing.py` similarly uses `TYPE_CHECKING` for the `Model` type.

The C extension (`_alea`) is imported with a try/except in modules that need it. If unavailable, `_alea` is `None` and C-backed operations raise `RuntimeError`.

## Build System

The package uses a custom `setup.py` that:

1. Compiles `libalea` via `make full` in the `csrc/libalea/` submodule
2. Builds the CPython extension (`_alea`) linking against `bin/libalea_full.a`
3. Installs the extension alongside the Python package

```bash
python3 setup.py build_ext --inplace   # Development build
pip install -e .                        # Editable install
```

The C binding uses the Python C API directly (no CFFI, no ctypes, no Cython). `aleathor_binding.c` is the extension entry point and includes feature-focused binding files such as `_bind_core.c`, `_bind_io.c`, `_bind_raycast.c`, `_bind_slice.c`, `_bind_mesh.c`, and `_bind_void.c`.

The binding includes headers for all library modules:
- `alea.h` — core system, cells, surfaces, CSG operations, simplification
- `alea_slice.h` — 2D slice curves and grid queries
- `alea_raycast.h` — ray tracing
- `alea_mesh.h` — structured mesh sampling and export (Gmsh/VTK)

Key binding sections (organized by C API area):
- **Cell info**: `get_cell`, `get_cell_by_index`, `get_cells` — includes lattice fields (`lat_type`, `lat_pitch`, `lat_fill`, etc.) when `lat_type != 0`
- **CSG simplification**: `flatten_all_cells` — full simplification pass returning stats dict
- **Primitive data**: `node_primitive_data` — returns type-specific dict for all 21 primitive types
- **BBox tightening**: `tighten_cell_bbox`, `tighten_all_bboxes`, `tighten_cell_bbox_numerical` — interval-arithmetic and numerical fallback
- **Mesh module**: `mesh_export` (one-shot export to file), `mesh_sample` (sample and return data)
