<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO

SPDX-License-Identifier: MPL-2.0
-->

# Tutorial

This tutorial walks through the main things you can do with AleaTHOR. Each section is self-contained: real code, real results, no hand-waving.

All examples assume:

```python
import aleathor as ath
```

## 1. Loading a Model

The most common starting point is an existing MCNP input file:

```python
model = ath.load("iter_blanket.inp")
print(model)
# Model: 4523 cells, 8901 surfaces, 12 universes
```

`ath.load()` auto-detects the format from the file extension. You can also be explicit:

```python
model = ath.read_mcnp("geometry.inp")
model = ath.read_openmc("geometry.xml")
```

Or load from a string:

```python
mcnp_input = """1 1 -10.0 -1
2 0 1

1 SO 5.0

"""
model = ath.read_mcnp_string(mcnp_input)
```

The returned `Model` is ready for queries — no separate build step required. The internal C system is built lazily on the first query.

## 2. Asking Questions About the Geometry

Once loaded, the most useful thing is asking "what's at this point?"

### Which cell contains a point?

```python
cell = model.cell_at(650.0, 0.0, 0.0)
if cell:
    print(f"Point is in cell {cell.id}, material {cell.material}")
else:
    print("Point is in void or undefined")
```

`cell_at()` returns a `Cell` — a lightweight read-only view of the cell. It traverses the full universe hierarchy: if the point is in a cell that has a FILL, it descends into the filled universe and continues until it finds a terminal cell (one with a material or void, not another fill).

Returns `None` if no cell claims the point. This means either void or a geometry error.

### What about the full hierarchy?

For debugging universe fills, you often want to see every cell the point passes through:

```python
cells = model.cells_at(650.0, 0.0, 0.0)
for cell in cells:
    print(f"  depth {cell.depth}: cell {cell.id}, "
          f"universe {cell.universe}, material {cell.material}")
```

Each `Cell` in the list has its `depth` set — 0 for root-level cells, increasing with nesting.

### Using subscript notation

If you know a cell ID:

```python
cell = model[10]      # Same as model.get_cell(10)
print(cell.material)  # Material number
print(cell.density)   # Density
print(cell.bounds)    # (xmin, xmax, ymin, ymax, zmin, zmax)
```

### Checking for overlaps

```python
overlaps = model.find_overlaps()
for c1, c2 in overlaps:
    print(f"Overlap: cell {c1.id} and cell {c2.id}")
```

This is statistical — it samples random points within cell bounding boxes and checks for multiple occupancy. Not exhaustive, but it catches the vast majority of real overlaps.

## 3. Filtering Cells

The `model.cells` property returns a `CellCollection` that supports filtering and iteration:

```python
# Iterate all cells
for cell in model.cells:
    print(f"Cell {cell.id}: material {cell.material}")

# Filter by material
tungsten = model.cells.by_material(74)
print(f"Found {len(tungsten)} tungsten cells")

# Filter by universe
blanket = model.cells.by_universe(100)

# Chain filters
u3_tungsten = model.cells.by_material(74).by_universe(3)

# Custom filter
dense = model.cells.filter(lambda c: abs(c.density) > 10.0)

# Get unique values
print(f"Materials: {model.cells.materials()}")
print(f"Universes: {model.cells.universes()}")
```

Individual cells are `Cell` objects with these properties: `id`, `material`, `density`, `universe`, `fill`, `name`, `bounds`, `is_void`, `is_filled`, `is_lattice`, `importance`, `depth`. Lattice cells also have `lattice_type`, `lattice_pitch`, `lattice_fill`, `lattice_dims`, and `lattice_lower_left`.

## 4. Visualizing the Geometry

AleaTHOR provides two complementary approaches to 2D visualization:

1. **Grid queries**: sample cell/material IDs on a pixel grid (fast, gives you colors)
2. **Analytical curves**: extract exact surface boundaries (lines, circles, ellipses — gives you contours)

The `plot()` method combines both:

```python
# XY cross-section at z=0
model.plot(z=0, bounds=(-100, 100, -100, 100))

# XZ cross-section at y=0
model.plot(y=0, bounds=(-100, 100, -100, 100))

# YZ cross-section at x=50
model.plot(x=50, bounds=(-100, 100, -100, 100))
```

### Coloring options

```python
# Color by material instead of cell ID
model.plot(z=0, bounds=(-100, 100, -100, 100), by_material=True)

# Add a colorbar
model.plot(z=0, bounds=(-100, 100, -100, 100),
           by_material=True, show_colorbar=True)

# Control contour boundaries
model.plot(z=0, bounds=(-100, 100, -100, 100),
           by_material=True, contour_by='cell')  # cell contours on material fill

# Hide contour lines
model.plot(z=0, bounds=(-100, 100, -100, 100), show_contours=False)
```

### Saving and multi-view

```python
# Save to file
model.plot(z=50, bounds=(-100, 100, -100, 100), save="midplane.png")

# Three orthogonal views at once
model.plot_views(bounds=(-100, 100, -100, 100, -100, 100),
                 save="views.png")
```

### Overlay tally data

```python
import numpy as np
flux = np.load("flux_z0.npy")

model.plot(z=0, bounds=(-100, 100, -100, 100),
           overlay=flux, overlay_norm='log',
           overlay_label='Flux [n/cm^2/s]', show_colorbar=True)
```

### Low-level access

If you need the raw data instead of a plot:

```python
# Grid: cell/material IDs at each pixel
grid = model.find_cells_grid_z(0, (-100, 100, -100, 100),
                                resolution=(800, 800))
cell_ids = grid['cell_ids']       # Flat list, 640000 elements
material_ids = grid['material_ids']

# Curves: exact surface boundaries
curves = model.get_slice_curves_z(0, (-100, 100, -100, 100))
for curve in curves['curves']:
    print(f"Surface {curve['surface_id']}: {curve['type']}")
```

The `universe_depth` parameter controls which level of the hierarchy you see:

- `-1` (default): innermost cell — what the transport code sees
- `0`: root-level cells only — useful for seeing the container structure
- `N`: cells at depth N

```python
grid = model.find_cells_grid_z(0, (-100, 100, -100, 100),
                                resolution=(400, 400),
                                universe_depth=0)
```

### Error detection

```python
grid = model.find_cells_grid_z(0, (-100, 100, -100, 100),
                                resolution=(200, 200),
                                detect_errors=True)
# grid['errors'] contains: 0=ok, 1=overlap, 2=undefined

# For thorough validation (rechecks every pixel):
errors = model.check_grid_overlaps(grid)
```

### Arbitrary planes

All slice functions support arbitrary cutting planes via `origin`, `normal`, and `up` vectors:

```python
# Diagonal slice through the origin
curves = model.get_slice_curves(
    origin=(0, 0, 0),
    normal=(1, 1, 0),     # 45-degree plane
    up=(0, 0, 1),
    bounds=(-100, 100, -100, 100)
)

# Plot on an arbitrary plane
model.plot(origin=(0, 0, 0), normal=(0, 1, 1), up=(1, 0, 0),
           bounds=(-100, 100, -100, 100))
```

## 5. Tracing Rays

Ray tracing reports every cell the ray passes through, in order:

```python
# Trace with direction
result = model.trace(origin=(0, 0, 0), direction=(1, 0, 0),
                     max_distance=500.0)

for seg in result:
    if seg.cell:
        print(f"Cell {seg.cell.id}: {seg.length:.2f} cm, "
              f"material {seg.material}")
    else:
        print(f"VOID: {seg.length:.2f} cm  <-- geometry error?")
```

### Point-to-point trace

```python
result = model.trace(start=(-500, 0, 0), end=(500, 0, 0))
```

### Path length through a material

```python
steel_path = result.path_length(material=26)
water_path = result.path_length(material=1)
total_path = result.path_length()  # all materials

print(f"Steel: {steel_path:.2f} cm, Water: {water_path:.2f} cm")
print(f"Materials hit: {result.materials_hit()}")
```

### Cell-aware tracing

For large models, cell-aware tracing is faster. Instead of testing every surface globally, it tracks through cells one at a time:

```python
result = model.trace(origin=(0, 0, 0), direction=(1, 0, 0),
                     cell_aware=True)
```

Same interface, same result, better performance on models with many surfaces.

## 6. Building Geometry from Scratch

You don't have to load from a file. You can build geometry with Python boolean operators:

```python
model = ath.Model("My Geometry")

# Create surfaces
sphere = ath.Sphere(0, 0, 0, radius=5.0)
box = ath.Box(-10, 10, -10, 10, -10, 10)

# Define regions using boolean operators
# -surface = inside (negative halfspace)
# +surface = outside (positive halfspace)
fuel = -sphere                # Inside sphere
moderator = -box & +sphere    # Inside box AND outside sphere

# Add cells
model.add_cell(fuel, material=1, density=-10.5, name="fuel")
model.add_cell(moderator, material=2, density=-1.0, name="moderator")
```

### Boolean operators

```python
region_a & region_b    # Intersection (AND)
region_a | region_b    # Union (OR)
region_a - region_b    # Difference (A AND NOT B)
~region                # Complement (NOT)
-surface               # Negative halfspace (inside)
+surface               # Positive halfspace (outside)
```

These correspond to MCNP cell definitions. For example, MCNP cell `1 1 -10.5 -1 2 -3` (intersection of -S1, +S2, -S3) becomes:

```python
region = -s1 & +s2 & -s3
```

And MCNP cell `2 0 -1:-2` (union of -S1 and -S2) becomes:

```python
region = -s1 | -s2
```

### Available surface types

| Class | Description | MCNP equivalent |
|-------|-------------|-----------------|
| `Plane(a, b, c, d)` | General plane | P |
| `XPlane(x0)` | Plane at x=x0 | PX |
| `YPlane(y0)` | Plane at y=y0 | PY |
| `ZPlane(z0)` | Plane at z=z0 | PZ |
| `Sphere(x0, y0, z0, radius)` | Sphere | S, SO |
| `CylinderX(y0, z0, radius)` | Infinite cylinder along X | CX, C/X |
| `CylinderY(x0, z0, radius)` | Infinite cylinder along Y | CY, C/Y |
| `CylinderZ(x0, y0, radius)` | Infinite cylinder along Z | CZ, C/Z |
| `ConeX(x0, y0, z0, t_sq)` | Cone along X | KX, K/X |
| `ConeY(x0, y0, z0, t_sq)` | Cone along Y | KY, K/Y |
| `ConeZ(x0, y0, z0, t_sq)` | Cone along Z | KZ, K/Z |
| `TorusX(x0, y0, z0, R, r)` | Torus with X axis | TX |
| `TorusY(x0, y0, z0, R, r)` | Torus with Y axis | TY |
| `TorusZ(x0, y0, z0, R, r)` | Torus with Z axis | TZ |
| `Box(xmin, xmax, ymin, ymax, zmin, zmax)` | Axis-aligned box | RPP |
| `RCC(bx, by, bz, hx, hy, hz, r)` | Right circular cylinder | RCC |
| `TRC(bx, by, bz, hx, hy, hz, r1, r2)` | Truncated cone | TRC |
| `ELL(v1, v2, major_axis)` | Ellipsoid | ELL |
| `REC(base, height, axis1, axis2)` | Right elliptical cylinder | REC |
| `WED(vertex, v1, v2, v3)` | Wedge | WED |
| `RHP(base, height, r1, r2, r3)` | Hexagonal prism | RHP/HEX |
| `GeneralBox(corner, v1, v2, v3)` | Oriented box | BOX |
| `Quadric(A, B, C, D, E, F, G, H, J, K)` | General quadric | GQ, SQ |

### Universe fills

```python
model = ath.Model()

# Define the fuel pin (universe 1)
fuel_surf = ath.Sphere(0, 0, 0, radius=0.5)
clad_surf = ath.CylinderZ(0, 0, radius=0.6)

model.add_cell(-fuel_surf, material=1, density=-10.0,
               universe=1, name="fuel")
model.add_cell(-clad_surf & +fuel_surf, material=2, density=-6.5,
               universe=1, name="clad")

# Container in universe 0 filled with universe 1
container = ath.Box(-5, 5, -5, 5, -5, 5)
model.add_cell(-container, fill=1, name="pin container")
```

## 7. Exporting

### To MCNP

```python
model.export_mcnp("output.inp")

# Or use the auto-detect:
model.save("output.inp")       # MCNP (default)
```

### To OpenMC

```python
model.export_openmc("geometry.xml")

# Or:
model.save("geometry.xml")     # OpenMC (from extension)
```

### To string

```python
mcnp_text = model.to_mcnp_string()
print(mcnp_text)
```

### Format conversion

Converting between MCNP and OpenMC is a two-liner:

```python
# MCNP to OpenMC
model = ath.load("input.inp")
model.save("geometry.xml")

# OpenMC to MCNP
model = ath.load("geometry.xml")
model.save("output.inp")
```

## 8. Manipulating Geometry

### Renumbering

```python
model.renumber_cells(start_id=1)       # Consecutive cell IDs from 1
model.renumber_surfaces(start_id=1)    # Consecutive surface IDs from 1
model.offset_cell_ids(10000)           # Add 10000 to all cell IDs
model.offset_surface_ids(10000)
model.offset_material_ids(100)
```

### Splitting and expanding

```python
# Split cells with top-level unions into simpler cells
n_new = model.split_union_cells()
print(f"Created {n_new} new cells")

# Expand macrobodies to primitive surfaces
n_expanded = model.expand_macrobodies()
```

### Tightening bounding boxes

```python
# Tighten all cell bounding boxes (improves query performance)
n_tightened = model.tighten_bboxes(tolerance=1.0)

# Tighten a single cell
bbox = model.tighten_cell_bbox(cell_id=10, tolerance=1.0)
# Returns (xmin, xmax, ymin, ymax, zmin, zmax)
```

### Volume estimation

```python
# Estimate cell volumes using random ray tracing
volumes = model.estimate_cell_volumes(n_rays=100000)
for i, vol in enumerate(volumes['volumes']):
    print(f"Cell index {i}: {vol:.3f} cm^3 "
          f"(+/- {volumes['rel_errors'][i]*100:.1f}%)")

# Remove tiny cells
removed = model.remove_cells_by_volume(volumes['volumes'], threshold=0.01)
print(f"Removed {removed} cells below 0.01 cm^3")
```

### Extracting sub-models

```python
# Extract a single universe
sub = model.extract_universe(universe_id=5)
sub.save("universe_5.inp")

# Extract by bounding box
region = model.extract_region((-50, 50, -50, 50, -50, 50))
region.save("region.inp")
```

### Flattening universes

```python
# Expand all fills — every cell ends up in a flat geometry
model.flatten_universe(0)
```

### CSG simplification

`simplify()` runs a full optimization pass on all cell CSG trees — eliminating complements, removing double negations, deduplicating subtrees, and more:

```python
stats = model.simplify()
print(f"Reduced {stats['nodes_before']} -> {stats['nodes_after']} CSG nodes")
print(f"Complements eliminated: {stats['complements_eliminated']}")
```

This is useful before export (smaller, cleaner CSG trees) and before analysis (faster queries).

### Numerical bounding box tightening

Interval arithmetic gives tight bounding boxes for most cells, but some complex boolean expressions produce loose bounds. The numerical fallback samples the cell boundary:

```python
# Tighten all cells with interval arithmetic first
model.tighten_bboxes()

# Then use numerical fallback on specific cells
model.tighten_bbox_numerical(cell_id=42)
```

## 9. Mesh Export

AleaTHOR can export the geometry as a structured hexahedral mesh, suitable for visualization in Gmsh or ParaView:

```python
# Export as Gmsh mesh
model.export_mesh("output.msh", nx=50, ny=50, nz=50)

# Export as VTK
model.export_mesh("output.vtk", nx=50, ny=50, nz=50, format="vtk")

# Specify bounds (otherwise auto-detected from cell bboxes)
model.export_mesh("region.msh", nx=100, ny=100, nz=100,
                  bounds=(-50, 50, -50, 50, -50, 50))
```

Each mesh element is assigned the material ID at its center. This is useful for:
- Quick 3D visualization of material assignments
- Verifying geometry before running transport
- Converting CSG geometry to a mesh representation

You can also sample without writing to file:

```python
result = model.sample_mesh(nx=50, ny=50, nz=50)
material_ids = result['material_ids']   # Flat list, Z-major order
cell_ids = result['cell_ids']
x_nodes = result['x_nodes']            # nx+1 node positions
```

## 10. Lattice Data

Lattice cells (rectangular or hexagonal) store their fill arrays and pitches. You can inspect them:

```python
for cell in model.cells:
    if cell.is_lattice:
        print(f"Cell {cell.id}: {cell.lattice_type} lattice")
        print(f"  Pitch: {cell.lattice_pitch}")
        print(f"  Dims: {cell.lattice_dims}")
        print(f"  Fill universes: {cell.lattice_fill[:10]}...")  # First 10
```

`lattice_type` is `1` for rectangular and `2` for hexagonal. `lattice_dims` gives `(imin, imax, jmin, jmax, kmin, kmax)`.

## 11. Working with Large Models

For models with deep universe hierarchies (tokamak-scale), build the spatial index before slicing:

```python
model = ath.load("iter_full.inp")
model.build_spatial_index()

# Now slicing and plotting work efficiently
model.plot(z=0, bounds=(-1000, 1000, -1000, 1000))
```

The spatial index builds a KD-tree over cell instances. Without it, queries scan cells linearly. With it, queries touch only the handful of cells whose bounding boxes contain the query point.

```python
# Check how many instances the index covers
print(f"Cell instances: {model.spatial_index_instance_count}")
```

### Configuration

System behavior is controlled through the `config` property:

```python
# Read current config
print(model.config)

# Update specific settings
model.config = {'log_level': 3, 'abs_tol': 1e-8}
```

### Logging

```python
ath.enable_logging()       # INFO level
ath.set_log_level(ath.LOG_DEBUG)  # More detail
ath.disable_logging()      # Silence

# Available levels: LOG_NONE, LOG_ERROR, LOG_WARN, LOG_INFO, LOG_DEBUG, LOG_TRACE
```

## Next Steps

- Read [Concepts](CONCEPTS.md) to understand surfaces, sense, regions, and the model lifecycle
- Read the [API Reference](API.md) for the complete class and function listing
- Read [Architecture](ARCHITECTURE.md) if you want to contribute or understand the internals
