# aleathor

<p align="center">
  <img src="assets/aleathor_logo.png" alt="aleathor Logo" width="400">
</p>

aleathor is a Python package for inspecting, debugging, and converting CSG
geometry models used in particle transport workflows.

The library is aimed at existing MCNP/OpenMC-style models: load a geometry,
query cells and materials, trace rays, inspect universe paths, plot slices, and
export or sample the model when needed.

aleathor is currently alpha software. Use it with independent checks before
relying on conversion or analysis results in production work.

## Status

Current capabilities include:

- MCNP input loading
- OpenMC XML loading, currently alpha
- point queries with `cell_at()`
- nested universe path queries with `cell_path_at()`
- ray tracing through cells and materials
- 2-D slice plotting
- structured mesh sampling and export
- MCNP, OpenMC, and Serpent export

See [Current Status](docs/STATUS.md) for the detailed state of the codebase.

## Installation

From a checkout:

```bash
git clone --recurse-submodules https://github.com/giovanni-mariano/aleathor.git
cd aleathor
pip install -e .
```

Requirements:

- Python >= 3.9
- C compiler, such as gcc or clang
- Make

`matplotlib` and `numpy` are installed as package dependencies.

### Build Options

| Variable | Default | Effect |
|----------|---------|--------|
| `PORTABLE` | `1` | When `0`, compile with `-march=native`. Do not distribute wheels built this way. |
| `USE_OPENMP` | `0` | When `1`, enable OpenMP if the compiler supports it. |

Example:

```bash
PORTABLE=0 USE_OPENMP=1 pip install -e .
```

For a source build from PyPI:

```bash
PORTABLE=0 USE_OPENMP=1 pip install --no-binary aleathor aleathor
```

## Quick Start

```python
import aleathor as ath

model = ath.load("model.inp")

print(model)
print(f"cells: {len(model.cells)}")
print(f"surfaces: {len(model.surfaces)}")

cell = model.cell_at(0.0, 0.0, 0.0)
if cell is not None:
    print(cell)
```

Trace a ray:

```python
trace = model.trace(start=(-100.0, 0.0, 0.0), end=(100.0, 0.0, 0.0))

for segment in trace:
    print(segment.cell, segment.length, segment.material)
```

Plot a slice:

```python
model.plot(z=0.0, bounds=(-100, 100, -100, 100), by_material=True)
```

Get raw slice data:

```python
grid = model.slice.grid(
    axis="z",
    value=0.0,
    bounds=(-100, 100, -100, 100),
    resolution=(300, 300),
)

curves = model.slice.curves(
    axis="z",
    value=0.0,
    bounds=(-100, 100, -100, 100),
)
```

## API Shape

`Model` owns the geometry. The common API stays on `Model`:

```python
model.cells
model.materials
model.surfaces

model.add_cell(...)
model.add_material(...)
model.cell_at(x, y, z)
model.cell_path_at(x, y, z)
model.trace(...)
model.plot(...)
model.save("out.inp")
```

Advanced operations are grouped under explicit namespaces:

```python
model.slice.grid(...)
model.slice.curves(...)
model.slice.labels(grid)

model.mesh.sample(...)
model.mesh.export("mesh.vtk", format="vtk")

model.analysis.find_overlaps()
model.analysis.estimate_cell_volumes()

model.repair.simplify()
model.repair.tighten_bboxes()

model.void.generate(...)
model.void.add(voids)

model.backend.config
```

Cells and materials are live views into the model. Mutating them updates the
backend model immediately:

```python
cell = model.cells[10]
cell.material = 2
cell.fill = 5
cell.fill = None

mat = model.get_material(1)
mat.density = 10.5
mat.add_nuclide(92235, 0.04)
```

Surfaces created in Python are immutable geometry definitions.

## Documentation

The documentation site is:

https://giovanni-mariano.github.io/aleathor/

Local documentation files:

- [Tutorial](docs/TUTORIAL.md)
- [Concepts](docs/CONCEPTS.md)
- [Architecture](docs/ARCHITECTURE.md)
- [API Reference](docs/API.md)
- [Current Status](docs/STATUS.md)

## Examples

See `examples/`:

- `basic_usage.py`
- `advanced_surfaces.py`
- `plotting_example.py`
- `plot_geometry.py`

## Development Note

This package was developed with support from AI tools.

## License

MPL-2.0
