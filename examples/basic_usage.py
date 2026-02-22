#!/usr/bin/env python3
"""
Basic usage example for AleaTHOR CSG library.

This example demonstrates how to:
1. Create geometry using Python surface objects
2. Query the geometry (point queries, ray tracing)
3. Filter and inspect cells
4. Export to MCNP format
"""

import aleathor as ath
from aleathor.surfaces import Sphere, CylinderZ, Box, Plane

# =============================================================================
# 1. Create a simple reactor pin cell geometry
# =============================================================================

print("=" * 60)
print("Creating a simple PWR fuel pin cell geometry")
print("=" * 60)

model = ath.Model("PWR Fuel Pin Cell")

# Define surfaces
# Fuel pellet (UO2) - radius 0.4096 cm
fuel_surf = Sphere(0, 0, 0, radius=0.4096)

# Fuel-clad gap (Helium) - radius 0.418 cm
gap_surf = CylinderZ(0, 0, radius=0.418)

# Cladding outer (Zircaloy) - radius 0.475 cm
clad_surf = CylinderZ(0, 0, radius=0.475)

# Water box (half-pitch 0.63 cm)
water_box = Box(-0.63, 0.63, -0.63, 0.63, -1.0, 1.0)

# Define cells using CSG region operators
# -surface means "inside" (negative side)
# +surface means "outside" (positive side)
# & means intersection, | means union

model.add_cell(
    region=-fuel_surf,
    material=1,
    density=10.5,
    name="fuel"
)

model.add_cell(
    region=-gap_surf & +fuel_surf,
    material=2,
    density=0.0001,
    name="gap"
)

model.add_cell(
    region=-clad_surf & +gap_surf,
    material=3,
    density=6.56,
    name="cladding"
)

model.add_cell(
    region=-water_box & +clad_surf,
    material=4,
    density=0.7,
    name="moderator"
)

# Print model summary
print(f"\n{model}")
print()

# =============================================================================
# 2. Inspect the geometry
# =============================================================================

print("=" * 60)
print("Inspecting the geometry")
print("=" * 60)

# List all cells
print("\nAll cells:")
for cell in model.cells:
    print(f"  Cell {cell.id}: {cell.name}")
    print(f"    Material: {cell.material}, Density: {cell.density} g/cm³")
    print(f"    Bounds: {cell.bounds}")

# Filter by material
print("\nCells with material 1 (fuel):")
for cell in model.cells.by_material(1):
    print(f"  Cell {cell.id}: {cell.name}")

# Get unique materials
print(f"\nUnique materials: {model.cells.materials()}")

# =============================================================================
# 3. Point queries
# =============================================================================

print("\n" + "=" * 60)
print("Point queries")
print("=" * 60)

test_points = [
    (0.0, 0.0, 0.0),    # Center (fuel)
    (0.41, 0.0, 0.0),   # In gap
    (0.45, 0.0, 0.0),   # In cladding
    (0.5, 0.0, 0.0),    # In moderator
    (1.0, 0.0, 0.0),    # Outside geometry
]

print("\nPoint-in-cell queries:")
for pt in test_points:
    cell = model.cell_at(*pt)
    if cell:
        print(f"  {pt} -> Cell {cell.id} ({cell.name}), material {cell.material}")
    else:
        print(f"  {pt} -> void (outside geometry)")

# Check if specific cell contains a point
fuel_cell = model.cells.get(1)
print(f"\nfuel_cell.contains(0, 0, 0) = {fuel_cell.contains(0, 0, 0)}")
print(f"fuel_cell.contains(0.5, 0, 0) = {fuel_cell.contains(0.5, 0, 0)}")

# =============================================================================
# 4. Ray tracing
# =============================================================================

print("\n" + "=" * 60)
print("Ray tracing")
print("=" * 60)

# Trace a ray from left to right through the pin
print("\nRay trace from (-1, 0, 0) to (1, 0, 0):")
trace = model.trace(start=(-1, 0, 0), end=(1, 0, 0))

print(f"  Total segments: {len(trace)}")
for i, seg in enumerate(trace):
    cell_name = seg.cell.name if seg.cell else "void"
    if seg.length < 1e10:  # Skip infinite segments
        print(f"  Segment {i+1}: {cell_name}, length={seg.length:.4f} cm, material={seg.material}")

# Calculate path length through fuel
fuel_path = trace.path_length(material=1)
print(f"\nTotal path length through fuel: {fuel_path:.4f} cm")

# Trace using origin + direction
print("\nRay trace from origin along +X direction:")
trace2 = model.trace(origin=(0, 0, 0), direction=(1, 0, 0), max_distance=1.0)
for seg in trace2:
    if seg.length < 1e10:
        print(f"  {seg.cell.name if seg.cell else 'void'}: {seg.length:.4f} cm")

# =============================================================================
# 5. Export to MCNP
# =============================================================================

print("\n" + "=" * 60)
print("Export to MCNP")
print("=" * 60)

output_file = "/tmp/pin_cell.inp"
model.save(output_file)
print(f"\nSaved to: {output_file}")

# Show the generated MCNP input
print("\nGenerated MCNP input:")
print("-" * 40)
with open(output_file) as f:
    content = f.read()
    print(content)

# =============================================================================
# 6. Working with surfaces
# =============================================================================

print("=" * 60)
print("Surface information")
print("=" * 60)

print(f"\nModel has {len(model.surfaces)} surfaces:")
for surf_id, surf in model.surfaces.items():
    print(f"  Surface {surf_id}: {type(surf).__name__}")

print("\n" + "=" * 60)
print("Example complete!")
print("=" * 60)
