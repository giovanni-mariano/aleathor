<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO

SPDX-License-Identifier: MPL-2.0
-->

# Current Status

aleathor is in alpha. The codebase is usable, but the public API is still being refined and some workflows need more validation on large production models.

## Available Now

- MCNP input parsing from files and strings.
- OpenMC XML parsing from files and strings.
- Programmatic CSG model construction with surfaces, regions, cells, materials, and universes.
- Point queries with `model.cell_at()` and hierarchy inspection with `model.cell_path_at()`.
- Ray tracing through geometry with `model.trace()`.
- Cell collection filtering by material, universe, fill, bounding box, and custom predicates.
- Axis-aligned and arbitrary-plane slice sampling.
- Matplotlib/numpy slice plotting, contours, labels, multi-view plots, and optional tally overlays. Plotting dependencies are installed with aleathor.
- MCNP, OpenMC XML, and Serpent export through `model.save()` and explicit export methods.
- Structured mesh sampling and export to Gmsh or VTK.
- Volume estimation, bounding-box tightening, ID renumbering, and void-region generation helpers.

## Experimental Or Still Evolving

- OpenMC and Serpent conversion paths need broader real-world validation.
- Overlap detection is sampling-based unless you use grid-level checks on a selected slice.
- Volume estimates are Monte Carlo estimates; results depend on the number of rays and bounding volume.
- Imported cells expose a reconstructed or C-backed region view, but not every piece of original input formatting is preserved.
- The Python API may still change before a stable release.

## Point Query Names

Use `model.cell_at(x, y, z)` for normal "what material/cell is here?" queries. It returns the terminal cell reached after following universe fills.

Use `model.cell_path_at(x, y, z)` when debugging nested universes. It returns the full containment path, ordered from the outermost cell to the innermost cell, with `cell.depth` set on each result.
