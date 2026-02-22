<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO

SPDX-License-Identifier: MPL-2.0
-->

# Concepts

This document explains the domain concepts that AleaTHOR works with. If you know MCNP, most of this will be familiar, but there are a few places where AleaTHOR's Python representation differs from the underlying C library or from what you might expect.

## Surfaces

A surface is a mathematical object that divides all of 3D space into two halves. A plane divides space into "above" and "below." A sphere divides space into "inside" and "outside."

In AleaTHOR, a surface is always one of these types:

| Class | Equation | MCNP equivalent |
|-------|----------|-----------------|
| `Plane` | ax + by + cz = d | P, PX, PY, PZ |
| `Sphere` | (x-x0)^2 + (y-y0)^2 + (z-z0)^2 = r^2 | S, SO, SX, SY, SZ |
| `CylinderX` | (y-y0)^2 + (z-z0)^2 = r^2 | CX, C/X |
| `CylinderY` | (x-x0)^2 + (z-z0)^2 = r^2 | CY, C/Y |
| `CylinderZ` | (x-x0)^2 + (y-y0)^2 = r^2 | CZ, C/Z |
| `ConeX/Y/Z` | ... = t^2 * (axis)^2 | KX, KY, KZ |
| `TorusX/Y/Z` | fourth-degree surface of revolution | TX, TY, TZ |
| `Box` | axis-aligned box | RPP |
| `Quadric` | Ax^2 + By^2 + ... + K = 0 | GQ, SQ |

Macrobodies (`RCC`, `TRC`, `ELL`, `REC`, `WED`, `RHP`, `GeneralBox`) are also surface types. They represent composite shapes (a cylinder-with-caps, a wedge, etc.) as single objects.

Surface IDs are auto-assigned when you create a surface. If you load from MCNP, the original surface IDs are preserved.

```python
s = ath.Sphere(0, 0, 0, radius=5.0)
print(s.id)  # Auto-assigned integer
```

## Sense and Halfspaces

Every surface divides space into two sides. In AleaTHOR (and MCNP), these sides are called **positive** and **negative** sense.

- **Negative sense** (inside): for a sphere, points closer to the center than the radius. For a plane, points on the side opposite to the normal vector. Written as `-S` in MCNP.
- **Positive sense** (outside): for a sphere, points farther from the center. Written as `+S` or just `S` in MCNP.

In Python, you get halfspaces using the `-` and `+` operators on surfaces:

```python
sphere = ath.Sphere(0, 0, 0, radius=5.0)

inside = -sphere    # Negative halfspace: points where r < 5
outside = +sphere   # Positive halfspace: points where r > 5
```

These return `Halfspace` objects — a type of `Region`. You can also use the `.interior()` and `.exterior()` methods:

```python
inside = sphere.interior()   # Same as -sphere
outside = sphere.exterior()  # Same as +sphere
```

## Regions and Boolean Operations

A `Region` defines a volume of space. The simplest region is a `Halfspace` (one side of a surface). Complex regions are built by combining halfspaces with boolean operations.

AleaTHOR uses Python operators for boolean combinations:

| Operation | Python | MCNP | Meaning |
|-----------|--------|------|---------|
| Intersection | `a & b` | `-1 -2` (juxtaposition) | Inside both a and b |
| Union | `a \| b` | `-1:-2` (colon) | Inside a or b or both |
| Difference | `a - b` | `-1 2` (sense flip) | Inside a but not b |
| Complement | `~a` | `#(...)` | Everything not inside a |

```python
sphere = ath.Sphere(0, 0, 0, radius=5.0)
box = ath.Box(-10, 10, -10, 10, -10, 10)

# Hemisphere: inside sphere AND above z=0
plane = ath.ZPlane(0)
hemisphere = -sphere & +plane

# Shell: inside box minus inside sphere
shell = -box & +sphere       # Same as: -box - (-sphere)

# Union of two spheres
s1 = ath.Sphere(-3, 0, 0, radius=2.0)
s2 = ath.Sphere(3, 0, 0, radius=2.0)
dumbbell = -s1 | -s2
```

### Region types

| Class | Description |
|-------|-------------|
| `Halfspace` | One side of a surface |
| `Intersection` | AND of multiple regions (flattened automatically) |
| `Union` | OR of multiple regions (flattened automatically) |
| `Complement` | NOT of a region (double negation eliminated) |

Nested operations of the same type are flattened: `(a & b) & c` becomes `Intersection([a, b, c])`, not a tree of two intersections.

### Region methods

Every `Region` supports:

```python
point in region        # Containment test: True if (x, y, z) is inside
region.get_surfaces()  # Set of all surfaces referenced by the region
```

## Cells

A cell is a region of space with a specific material assignment. It corresponds to an MCNP cell card.

Every cell has:

- **Cell ID**: the MCNP cell number (unique identifier, positive integer)
- **Region**: a `Region` defining the shape
- **Material**: material number (0 for void)
- **Density**: material density (negative = g/cm^3, positive = atoms/barn-cm)
- **Universe**: which universe the cell belongs to (default 0)
- **Fill**: optionally, a universe to fill this cell with
- **Importance**: particle importance for transport (default 1.0)

A cell either has a material or a fill, never both. A cell with a fill is a "container" — it defines a region of space and says "look in universe N to find what's actually here."

```python
model = ath.Model()
sphere = ath.Sphere(0, 0, 0, radius=5.0)

cell = model.add_cell(
    region=-sphere,
    material=1,
    density=-10.5,
    name="fuel"
)
```

### _CellData vs Cell

AleaTHOR has two cell representations:

| Type | Purpose | Source |
|------|---------|--------|
| `_CellData` | Internal dataclass holding Python Region objects | `model.add_cell()` (internal) |
| `Cell` | Mutable view backed by the C system | `model.get_cell()`, `model.cell_at()`, iteration |

When you **build** geometry in Python, cells are stored as `_CellData` dataclass instances in `model._cells`. When you **query** geometry, the C system returns `Cell` objects.

`Cell` is the public API type. It provides access to all cell properties, with mutation support for key fields:

```python
cell = model.get_cell(10)     # Returns Cell
cell.id                        # Cell ID
cell.material                  # Material number (read-write)
cell.density                   # Density (read-write)
cell.universe                  # Universe ID
cell.fill                      # Fill universe or None (read-write)
cell.bounds                    # Bounding box tuple
cell.is_void                   # True if material == 0
cell.is_filled                 # True if fill is set
cell.is_lattice                # True if this cell defines a lattice
cell.contains(x, y, z)        # Point containment
```

To modify cell properties, assign directly or use `model.update_cell()`:

```python
cell.material = 2              # Direct mutation
cell.density = 5.0
model.update_cell(10, material=2, density=5.0)  # Also works
```

Mutations mark the model dirty so the C system is rebuilt on the next query.

## Universes

A universe is a collection of cells that share the same coordinate system. Every cell belongs to exactly one universe. By default, everything is in universe 0.

Universes become important when you have repeated geometry. Instead of defining the same fuel pin 1000 times, you define it once in universe 5, and then fill 1000 container cells with universe 5.

### Fills

A cell with `fill=N` is a container. When AleaTHOR evaluates a point query, it:

1. Finds the container cell in the current universe
2. Sees it has a fill
3. Applies the inverse of the fill transform (if any)
4. Searches for the point among the cells of universe N
5. Repeats if it hits another fill

This is recursive. AleaTHOR handles arbitrary nesting depth.

```python
# Define pin geometry in universe 1
model.add_cell(-fuel_surf, material=1, density=-10.0, universe=1)
model.add_cell(-clad_surf & +fuel_surf, material=2, density=-6.5, universe=1)

# Container in universe 0, filled with universe 1
model.add_cell(-box, fill=1)
```

### Lattices

A lattice cell is a cell that tiles space with a repeating pattern of universes. MCNP supports two lattice types:

| Type | Value | Description |
|------|-------|-------------|
| Rectangular | 1 | Elements arranged on a Cartesian grid |
| Hexagonal | 2 | Elements arranged on a hexagonal grid |

A lattice cell defines a pitch (element spacing in each dimension), a fill array (which universe goes in each element), and dimensions (the range of lattice indices).

```python
for cell in model.cells:
    if cell.is_lattice:
        print(f"Cell {cell.id}: type={cell.lattice_type}")
        print(f"  Pitch: {cell.lattice_pitch}")        # (dx, dy, dz)
        print(f"  Dims: {cell.lattice_dims}")           # (imin,imax,jmin,jmax,kmin,kmax)
        print(f"  Fill: {cell.lattice_fill[:5]}...")     # Universe IDs
        print(f"  Lower-left: {cell.lattice_lower_left}")
```

When AleaTHOR evaluates a point inside a lattice cell, it determines which lattice element the point falls in (based on pitch and lower-left corner), looks up the universe for that element in the fill array, and searches for the point among that universe's cells.

### set_fill

You can change a cell's fill after creation:

```python
model.set_fill(cell_id=10, fill_universe=5, transform=0)
```

This operates on the C system directly and takes effect immediately (no rebuild needed).

## Materials

A material is a composition of nuclides or elements. AleaTHOR stores:

- **Material ID**: the MCNP material number
- **Name**: optional human-readable name
- **Density**: stored on the cell, not the material
- **Composition**: optional nuclide fractions

Materials are preserved through loading and export. AleaTHOR doesn't do physics with materials — it carries them along so queries return the correct material number and exports produce valid input files.

```python
model.add_material(1, name="UO2", density=-10.5)
model.add_material(2, name="Steel")

# Create a mixture
mix_id = model.create_mixture([1, 2], [0.7, 0.3], name="fuel-steel-mix")
```

## Model Lifecycle

The `Model` class manages a dual representation:

1. **Python objects**: `Cell` dataclass instances in `model._cells`, `Surface` objects in `model._surfaces`
2. **C system**: the compiled geometry in `model._sys`, used for all queries

### The dirty flag

When you modify cells (add, remove, update), the Python objects are changed immediately but the C system becomes **stale**. The `_dirty` flag tracks this:

```python
model.add_cell(...)      # Sets _dirty = True
model.cell_at(0, 0, 0)  # Triggers _rebuild_if_needed() — rebuilds C system
```

The rebuild is automatic and transparent. You never need to call it manually.

### Loaded models

When you load a model from an MCNP file, the C system is populated directly by the parser. Cells use `_ImportedRegion` objects instead of Python `Region` trees — these are lightweight wrappers that reference the C system's CSG nodes. The Python `_cells` dict is empty for loaded models.

This means loaded models are efficient: they don't duplicate the geometry in Python objects.

## Macrobodies

Macrobodies are composite surfaces that MCNP defines as single surface cards but that are really combinations of simpler surfaces:

| Macrobody | Description | Decomposition |
|-----------|-------------|---------------|
| `RCC` | Right circular cylinder | Cylinder + 2 planes |
| `Box` / RPP | Axis-aligned box | 6 planes (or single primitive) |
| `TRC` | Truncated right cone | Cone + 2 planes |
| `ELL` | Ellipsoid | Quadric surface |
| `REC` | Right elliptical cylinder | Quadric + 2 planes |
| `WED` | Wedge | 5 planes |
| `RHP` / HEX | Hexagonal prism | 8 planes |
| `GeneralBox` | Oriented box | 6 planes |

You can expand macrobodies to their primitive surfaces:

```python
n = model.expand_macrobodies()
print(f"Expanded {n} macrobodies")
```

## Boundary Conditions

Surfaces can have boundary conditions that affect particle behavior in transport:

| Type | Meaning |
|------|---------|
| `'transmissive'` | Particles cross normally (default) |
| `'reflective'` | Particles reflect at the surface |
| `'vacuum'` | Particles are killed |

Boundary conditions are preserved through loading and export. They don't affect point queries or ray tracing (those are purely geometric operations).

```python
s = ath.Sphere(0, 0, 0, radius=100.0, boundary='reflective')
```

## Configuration

System behavior is controlled through the `config` property:

```python
cfg = model.config       # Returns dict
model.config = {         # Update specific keys
    'abs_tol': 1e-8,
    'log_level': 3,
}
```

Key configuration fields:

| Key | Default | Meaning |
|-----|---------|---------|
| `abs_tol` | 1e-6 | Absolute tolerance for surface matching |
| `rel_tol` | 1e-9 | Relative tolerance |
| `log_level` | 0 | 0=none, 1=error, 2=warn, 3=info, 4=debug, 5=trace |
| `dedup` | true | Enable surface deduplication |
| `export_materials` | true | Include material cards in export |

## Error Handling

AleaTHOR functions raise standard Python exceptions:

| Exception | When |
|-----------|------|
| `KeyError` | Cell or material not found |
| `ValueError` | Invalid arguments (negative radius, missing bounds, etc.) |
| `RuntimeError` | C extension not available, internal error |
| `FileNotFoundError` | Input file not found |
| `KeyboardInterrupt` | Ctrl+C during long operations (ray tracing, grid queries) |

The C library supports cooperative interruption: long-running operations check a flag periodically. When you press Ctrl+C during a grid query or ray trace, the operation terminates cleanly and Python raises `KeyboardInterrupt`.

## Plotting Constants

For working with raw grid data:

| Constant | Value | Meaning |
|----------|-------|---------|
| `CELL_VOID` | -1 | Grid pixel is in void |
| `CELL_UNDEFINED` | -2 | No cell claims this pixel |
| `CELL_OVERLAP` | -3 | Multiple cells claim this pixel |
| `GRID_ERROR_OK` | 0 | No error at this pixel |
| `GRID_ERROR_OVERLAP` | 1 | Overlap detected |
| `GRID_ERROR_UNDEFINED` | 2 | Undefined region |
