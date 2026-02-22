<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO

SPDX-License-Identifier: MPL-2.0
-->

# Architecture

This document explains how the AleaTHOR Python package works internally. Read this if you want to contribute, understand the design decisions, or debug issues at the boundary between Python and C.

## The Big Picture

AleaTHOR Python is a thin layer over a C geometry engine ([libalea](../csrc/libalea/)). The package structure:

```
src/aleathor/
    __init__.py        Public API surface, optional imports
    model.py           Model class — main container
    collections.py     Cell, CellCollection, TraceResult, TraceSegment
    surfaces.py        Surface classes (Sphere, Box, CylinderZ, etc.)
    geometry.py        Region classes (Halfspace, Intersection, Union, Complement)
    io.py              File I/O (read_mcnp, write_openmc, etc.)
    slicing.py         DRY helper for slice parameter extraction
    plotting.py        matplotlib visualization (optional dependency)
    _binding/
        aleathor_binding.c   CPython extension module
```

The C extension (`_alea`) is a CPython module that wraps `libalea` functions into Python-callable methods on a `System` class. The Python `Model` class owns a `System` instance (`model._sys`) and manages the translation between Python objects and C data.

## Dual Representation

The `Model` class maintains two parallel representations of the geometry:

```
Model
├── Python side                     C side
│   ├── _cells: Dict[int, Cell]     _sys: _alea.System
│   ├── _surfaces: Dict[int, Surface]   (contains all nodes, primitives,
│   ├── _materials: Dict[int, Material]   cells, surfaces, materials,
│   └── _universes: Dict[int, Universe]   acceleration structures)
│
└── _dirty: bool   ← links the two sides
```

### When you build geometry in Python

```python
model.add_cell(region=-sphere, material=1, density=-10.5)
```

1. A `Cell` dataclass is created and stored in `model._cells`
2. The `_dirty` flag is set to `True`
3. No C data is created yet

### When you query the geometry

```python
cell = model.cell_at(0, 0, 0)
```

1. `_rebuild_if_needed()` checks `_dirty`
2. If dirty: resets `_sys`, iterates `_cells`, calls `region._to_csg(model)` to build the C CSG tree, adds each cell to the C system, builds the universe index
3. Clears `_dirty`
4. Calls `_sys.find_cell(x, y, z)` — pure C computation
5. Returns a `Cell` wrapping the result

### When you load from a file

```python
model = ath.read_mcnp("geometry.inp")
```

The C parser builds the system directly. The Python `_cells` dict stays **empty**. Cells use `_ImportedRegion` objects — lightweight wrappers that reference C system nodes by index. The model is **not dirty** (the C system is already built).

This is important: loaded models don't duplicate geometry in Python. The `Cell` objects returned by queries go straight to the C system.

## The _rebuild_if_needed() Cycle

```
add_cell() / remove_cell() / update_cell()
        │
        ▼
    _dirty = True
        │
        ▼ (on next query)
    _rebuild_if_needed()
        │
        ├── _sys.reset()
        ├── _halfspace_node_cache.clear()
        ├── for each cell in _cells:
        │       node_id = cell.region._to_csg(model)
        │       _sys.add_cell(cell_id, node_id, material, density, universe)
        ├── _sys.build_universe_index()
        └── _dirty = False
```

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

## _CellData vs Cell

| | `_CellData` | `Cell` |
|--|--------|------------|
| Type | `@dataclass` | Regular class |
| Storage | Python `_cells` dict | Created on demand |
| Backed by | Python attributes | C system (by index) |
| Mutable | Via direct attribute access | Via property setters (marks model dirty) |
| When | Building geometry | Querying and mutating geometry |

`Cell` (in `collections.py`) is the **public API type**. Users should never interact with `_CellData` directly. A `Cell` holds a reference to the `Model` and an integer index into the C system's cell array. Property access (`cell.material`, `cell.bounds`) calls into the C system. Property setters (`cell.material = 2`) update the Python `_CellData` and mark the model dirty.

```python
class Cell:
    def __init__(self, model, index, _depth=0):
        self._model = model
        self._index = index
        self._depth = _depth

    @property
    def material(self):
        return self._model._sys.cell_material(self._index)
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

`plotting.py` is separated from `model.py` because matplotlib is an optional dependency. The module is only imported when plotting is requested.

The key functions are:

- `plot()`: high-level orchestrator. Calls the appropriate `find_cells_grid_*()` and `get_slice_curves_*()`, then delegates to `plot_slice_filled()` or `plot_slice_curves()`.
- `plot_views()`: creates a 1x3 subplot figure with XY, XZ, YZ slices.
- `plot_slice_filled()`: the workhorse. Takes a grid result, creates an image with cell/material coloring, overlays contour lines, handles error highlighting, and supports tally overlay.

Model delegates to these via lazy imports:

```python
def plot(self, z=None, y=None, x=None, **kwargs):
    from .plotting import plot as _plot
    return _plot(self, z=z, y=y, x=x, **kwargs)
```

## _ImportedRegion

When geometry is loaded from an MCNP file, cells reference `_ImportedRegion` objects instead of Python `Region` trees. An `_ImportedRegion` is a thin wrapper around a C system and cell index:

```python
class _ImportedRegion(Region):
    def __init__(self, system, cell_index):
        self._system = system
        self._cell_index = cell_index

    def _to_csg(self, model):
        # Returns the existing C node ID — no rebuild needed
        return self._system.cell_root_node(self._cell_index)
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
    └── plotting.py ── plot, plot_views, plot_slice_filled, ... (optional: matplotlib)
```

Dependencies flow downward. `plotting.py` depends on `model.py` (via `Model` type hints) but only through `TYPE_CHECKING` — no runtime import cycle. `slicing.py` similarly uses `TYPE_CHECKING` for the `Model` type.

The C extension (`_alea`) is imported at the top of `model.py` with a try/except. If unavailable, `_csg` is `None` and any geometry operation raises `RuntimeError`.

## Build System

The package uses a custom `setup.py` that:

1. Compiles `libalea` via `make lib` in the `csrc/libalea/` submodule
2. Builds the CPython extension (`_alea`) linking against the static library
3. Installs the extension alongside the Python package

```bash
python3 setup.py build_ext --inplace   # Development build
pip install -e .                        # Editable install
```

The C binding (`_binding/aleathor_binding.c`) is a single file that wraps every public `libalea` function. It uses the Python C API directly (no CFFI, no ctypes, no Cython).

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
