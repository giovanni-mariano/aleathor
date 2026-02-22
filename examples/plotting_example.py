#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""
Plotting example for AleaTHOR CSG library.

This example demonstrates the visualization capabilities:
- 2D slice plots (XY, XZ, YZ planes)
- Arbitrary plane slices (diagonal, tilted)
- Color by cell or material
- Multi-view layouts
- Saving plots to files
"""

import aleathor as ath
from aleathor.surfaces import Sphere, CylinderZ, Box, TorusZ

# =============================================================================
# Create a tokamak-like geometry for visualization
# =============================================================================

print("=" * 60)
print("Creating tokamak-like geometry for plotting demo")
print("=" * 60)

model = ath.Model("Tokamak Demo")

# Plasma chamber (torus)
plasma_torus = TorusZ(0, 0, 0, major_radius=6.0, minor_radius=2.0)

# First wall (slightly larger torus)
first_wall = TorusZ(0, 0, 0, major_radius=6.0, minor_radius=2.3)

# Vacuum vessel (cylinder)
vessel_inner = CylinderZ(0, 0, radius=3.5)
vessel_outer = CylinderZ(0, 0, radius=9.0)

# Outer boundary
boundary = Box(-12, 12, -12, 12, -5, 5)

# Define cells
model.add_cell(
    region=-plasma_torus,
    material=0,  # Void (plasma region)
    density=0,
    name="plasma"
)

model.add_cell(
    region=-first_wall & +plasma_torus,
    material=1,
    density=7.9,  # Steel first wall
    name="first_wall"
)

model.add_cell(
    region=-vessel_outer & +vessel_inner & +first_wall,
    material=2,
    density=2.5,  # Concrete shielding
    name="shielding"
)

model.add_cell(
    region=-boundary & +vessel_outer,
    material=3,
    density=0.001,  # Air
    name="air"
)

model.add_cell(
    region=-vessel_inner & +first_wall,
    material=3,
    density=0.001,
    name="inner_air"
)

print(f"\nCreated model with {len(model.cells)} cells")

# =============================================================================
# Basic slice plot
# =============================================================================

print("\n" + "=" * 60)
print("Generating plots...")
print("=" * 60)

# Plot XY slice at z=0 (midplane)
print("\n1. XY slice at z=0 (midplane view)")
model.plot(z=0, bounds=(-12, 12, -12, 12))

# =============================================================================
# Plot colored by material
# =============================================================================

print("\n2. XY slice colored by material with colorbar")
model.plot(z=0, bounds=(-12, 12, -12, 12), by_material=True, show_colorbar=True)

# =============================================================================
# Different slice planes
# =============================================================================

print("\n3. XZ slice at y=0 (vertical cross-section)")
model.plot(y=0, bounds=(-12, 12, -5, 5))

print("\n4. YZ slice at x=6 (through torus)")
model.plot(x=6, bounds=(-12, 12, -5, 5))

# =============================================================================
# Save plots to files
# =============================================================================

print("\n5. Saving plots to files...")

model.plot(z=0, bounds=(-12, 12, -12, 12), save="/tmp/tokamak_xy.png")
print("   Saved: /tmp/tokamak_xy.png")

model.plot(y=0, bounds=(-12, 12, -5, 5), save="/tmp/tokamak_xz.png")
print("   Saved: /tmp/tokamak_xz.png")

model.plot(z=0, bounds=(-12, 12, -12, 12), by_material=True, save="/tmp/tokamak_materials.png")
print("   Saved: /tmp/tokamak_materials.png")

# =============================================================================
# Three orthogonal views
# =============================================================================

print("\n6. Three orthogonal views in one figure")
model.plot_views(bounds=(-12, 12, -12, 12, -5, 5), save="/tmp/tokamak_views.png")
print("   Saved: /tmp/tokamak_views.png")

# =============================================================================
# Arbitrary plane slices
# =============================================================================

print("\n7. Arbitrary plane slices")

# Diagonal plane through the tokamak
# origin: point on the plane
# normal: perpendicular to the plane
# up: defines the vertical direction in the plot

# 45-degree diagonal slice (XY plane rotated around Z)
import math
angle = math.radians(45)
print(f"\n   45-degree diagonal slice:")
model.plot(
    origin=(0, 0, 0),
    normal=(math.sin(angle), math.cos(angle), 0),  # Normal at 45 deg in XY
    up=(0, 0, 1),  # Z is up
    bounds=(-12, 12, -5, 5),
    save="/tmp/tokamak_diagonal.png"
)
print("   Saved: /tmp/tokamak_diagonal.png")

# Tilted plane (looking at tokamak from an angle)
print(f"\n   Tilted plane (30-degree tilt from horizontal):")
tilt = math.radians(30)
model.plot(
    origin=(0, 0, 0),
    normal=(0, math.sin(tilt), math.cos(tilt)),  # Tilted in YZ plane
    up=(0, math.cos(tilt), -math.sin(tilt)),  # Perpendicular up vector
    bounds=(-12, 12, -12, 12),
    save="/tmp/tokamak_tilted.png"
)
print("   Saved: /tmp/tokamak_tilted.png")

# Slice through a specific point at an angle
print(f"\n   Slice through point (6, 0, 0) at 30-degree angle:")
model.plot(
    origin=(6, 0, 0),  # Point on the torus
    normal=(math.cos(math.radians(30)), math.sin(math.radians(30)), 0),
    up=(0, 0, 1),
    bounds=(-5, 5, -5, 5),
    save="/tmp/tokamak_local.png"
)
print("   Saved: /tmp/tokamak_local.png")

# =============================================================================
# Counting cells, surfaces, materials in a slice
# =============================================================================

print("\n" + "=" * 60)
print("Slice statistics (cells, surfaces, materials)")
print("=" * 60)

# Get analytical data for a slice
bounds = (-12, 12, -12, 12)
curves = model.get_slice_curves_z(0, bounds)
grid = model.find_cells_grid_z(0, bounds, resolution=(100, 100))

# Count using the helper functions
from aleathor import count_cells, count_surfaces, count_materials, get_slice_stats

print(f"\nXY slice at z=0:")
print(f"  Unique cells: {count_cells(grid)}")
print(f"  Unique surfaces: {count_surfaces(curves)}")
print(f"  Unique materials: {count_materials(grid)}")

# Or get all stats at once
stats = get_slice_stats(curves, grid)
print(f"\n  All stats: {stats}")

# =============================================================================
# Ray trace visualization
# =============================================================================

print("\n" + "=" * 60)
print("Ray trace analysis")
print("=" * 60)

# Trace through the tokamak
trace = model.trace(start=(-12, 0, 0), end=(12, 0, 0))

print("\nRay trace from (-12, 0, 0) to (12, 0, 0):")
print("-" * 50)
for seg in trace:
    if seg.length < 1e10:
        name = seg.cell.name if seg.cell else "void"
        mat = seg.material
        print(f"  {name:15s}: {seg.length:8.3f} cm  (material {mat})")

# Path lengths by material
print("\nPath lengths by material:")
for mat in sorted(trace.materials_hit()):
    path = trace.path_length(material=mat)
    print(f"  Material {mat}: {path:.3f} cm")

print("\n" + "=" * 60)
print("Plotting example complete!")
print("=" * 60)
