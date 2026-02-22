#!/usr/bin/env python3
"""
Advanced surfaces example for AleaTHOR CSG library.

This example demonstrates the various surface types available:
- Basic surfaces: Plane, Sphere, Cylinder, Cone
- Torus surfaces
- Macrobodies: RCC, TRC, Box
- General quadric surface
"""

import aleathor as ath
from aleathor.surfaces import (
    # Basic surfaces
    Plane, XPlane, YPlane, ZPlane,
    Sphere,
    CylinderX, CylinderY, CylinderZ,
    ConeX, ConeY, ConeZ,
    # Torus
    TorusX, TorusY, TorusZ,
    # Macrobodies
    Box, RCC, TRC,
    # General
    Quadric,
)

print("=" * 60)
print("AleaTHOR Surface Types Demo")
print("=" * 60)

# =============================================================================
# 1. Planes
# =============================================================================

print("\n--- Planes ---")

model = ath.Model("Plane Demo")

# General plane: ax + by + cz = d
plane = Plane(a=1, b=1, c=0, d=5)  # x + y = 5

# Axis-aligned planes (convenient shortcuts)
px = XPlane(x=10)   # Plane at x=10
py = YPlane(y=-5)   # Plane at y=-5
pz = ZPlane(z=0)    # Plane at z=0

print(f"  General plane: {plane}")
print(f"  X-plane at x=10: {px}")
print(f"  Y-plane at y=-5: {py}")
print(f"  Z-plane at z=0: {pz}")

# =============================================================================
# 2. Spheres
# =============================================================================

print("\n--- Spheres ---")

# Sphere centered at origin
s1 = Sphere(0, 0, 0, radius=5)

# Sphere at arbitrary position
s2 = Sphere(10, 5, 3, radius=2.5)

print(f"  Sphere at origin, r=5: {s1}")
print(f"  Sphere at (10,5,3), r=2.5: {s2}")

# =============================================================================
# 3. Cylinders
# =============================================================================

print("\n--- Cylinders ---")

# Cylinders along each axis
cyl_x = CylinderX(y0=0, z0=0, radius=3)  # Along X axis
cyl_y = CylinderY(x0=0, z0=0, radius=3)  # Along Y axis
cyl_z = CylinderZ(x0=0, y0=0, radius=3)  # Along Z axis

# Off-center cylinder
cyl_off = CylinderZ(x0=5, y0=5, radius=2)

print(f"  Cylinder along X: {cyl_x}")
print(f"  Cylinder along Y: {cyl_y}")
print(f"  Cylinder along Z: {cyl_z}")
print(f"  Off-center cylinder: {cyl_off}")

# =============================================================================
# 4. Cones
# =============================================================================

print("\n--- Cones ---")

# Cones along each axis
# t_sq is tan²(half-angle)
cone_x = ConeX(x0=0, y0=0, z0=0, t_sq=0.25)  # 26.57° half-angle
cone_y = ConeY(x0=0, y0=0, z0=0, t_sq=0.25)
cone_z = ConeZ(x0=0, y0=0, z0=0, t_sq=1.0)   # 45° half-angle

print(f"  Cone along X (26.57° half-angle): {cone_x}")
print(f"  Cone along Y (26.57° half-angle): {cone_y}")
print(f"  Cone along Z (45° half-angle): {cone_z}")

# =============================================================================
# 5. Torus
# =============================================================================

print("\n--- Torus ---")

# Torus (donut shape) along each axis
# major_radius = distance from center to tube center
# minor_radius = radius of the tube
torus_z = TorusZ(0, 0, 0, major_radius=5, minor_radius=1)
torus_x = TorusX(0, 0, 0, major_radius=5, minor_radius=1)
torus_y = TorusY(0, 0, 0, major_radius=5, minor_radius=1)

print(f"  Torus along Z (R=5, r=1): {torus_z}")
print(f"  Torus along X (R=5, r=1): {torus_x}")
print(f"  Torus along Y (R=5, r=1): {torus_y}")

# =============================================================================
# 6. Macrobodies
# =============================================================================

print("\n--- Macrobodies ---")

# Box (RPP - Rectangular Parallelepiped)
box = Box(xmin=-5, xmax=5, ymin=-3, ymax=3, zmin=-10, zmax=10)

# RCC (Right Circular Cylinder)
# Defined by base point, height vector, and radius
rcc = RCC(
    base_x=0, base_y=0, base_z=-5,
    height_x=0, height_y=0, height_z=10,
    radius=2
)

# TRC (Truncated Right Cone)
# Cone frustum with different radii at base and top
trc = TRC(
    base_x=0, base_y=0, base_z=0,
    height_x=0, height_y=0, height_z=10,
    base_radius=3, top_radius=1
)

print(f"  Box: {box}")
print(f"  RCC (cylinder): {rcc}")
print(f"  TRC (truncated cone): {trc}")

# =============================================================================
# 7. General Quadric
# =============================================================================

print("\n--- General Quadric ---")

# General quadric: Ax² + By² + Cz² + Dxy + Eyz + Fzx + Gx + Hy + Jz + K = 0
# Example: Ellipsoid x²/4 + y²/9 + z²/16 = 1
# Rewrite: x²/4 + y²/9 + z²/16 - 1 = 0
# Multiply by 36: 9x² + 4y² + 2.25z² - 36 = 0
ellipsoid = Quadric(
    A=9, B=4, C=2.25,   # x², y², z² coefficients
    D=0, E=0, F=0,       # xy, yz, zx coefficients
    G=0, H=0, J=0,       # x, y, z coefficients
    K=-36                # constant
)

print(f"  Ellipsoid via quadric: {ellipsoid}")

# =============================================================================
# 8. Build a geometry with various surfaces
# =============================================================================

print("\n" + "=" * 60)
print("Building geometry with torus")
print("=" * 60)

model = ath.Model("Torus Example")

# Create a torus inside a box
torus = TorusZ(0, 0, 0, major_radius=5, minor_radius=1.5)
outer = Box(-10, 10, -10, 10, -5, 5)

model.add_cell(region=-torus, material=1, density=8.0, name="torus_interior")
model.add_cell(region=-outer & +torus, material=2, density=1.0, name="surrounding")

print(f"\n{model}")

# Query points
print("\nPoint queries:")
test_points = [(5, 0, 0), (0, 0, 0), (8, 0, 0)]
for pt in test_points:
    cell = model.cell_at(*pt)
    if cell:
        print(f"  {pt} -> {cell.name}")
    else:
        print(f"  {pt} -> void")

# Ray trace through torus
print("\nRay trace through torus (y=0, z=0):")
trace = model.trace(start=(-10, 0, 0), end=(10, 0, 0))
for seg in trace:
    if seg.length < 1e10:
        name = seg.cell.name if seg.cell else "void"
        print(f"  {name}: {seg.length:.3f} cm")

# Export
output = "/tmp/torus_example.inp"
model.save(output)
print(f"\nExported to: {output}")

print("\n" + "=" * 60)
print("Surface demo complete!")
print("=" * 60)
