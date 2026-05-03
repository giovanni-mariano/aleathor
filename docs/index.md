<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO

SPDX-License-Identifier: MPL-2.0
-->

# aleathor

aleathor is a Python library for debugging and analyzing large Constructive Solid Geometry (CSG) models used in neutron and gamma transport simulations.

It focuses on inspecting existing MCNP and OpenMC-style geometry models, finding geometry issues, tracing rays, querying cells, and generating cross-section plots for transport workflows.

## Start Here

- [Current Status](STATUS.md): what is available now and what is still experimental.
- [Tutorial](TUTORIAL.md): load a model, query geometry, plot slices, and export results.
- [Concepts](CONCEPTS.md): user-facing geometry and API mental model.
- [API Reference](API.md): find public classes, functions, and methods grouped by task.
- [Architecture](ARCHITECTURE.md): contributor-facing internals and C/Python design.

## Installation

```bash
git clone --recurse-submodules https://github.com/giovanni-mariano/aleathor.git
cd aleathor
pip install -e .
```

## Quick Example

```python
import aleathor as ath

model = ath.load("geometry.inp")

cell = model.cell_at(100.0, 0.0, 0.0)
if cell:
    print(f"Point is in cell {cell.id}, material {cell.material}")

model.plot(z=0, bounds=(-100, 100, -100, 100))
```

## Project Status

aleathor is under active development. APIs and behavior may change while the project is still in its alpha stage.
