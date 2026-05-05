#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""
Material definition example.

Demonstrates how to define materials with nuclide/element composition,
create mixtures, and use them in a geometry model.
"""

import aleathor as ath

# =============================================================================
# Create a model
# =============================================================================

model = ath.Model("PWR Fuel Pin")

# =============================================================================
# Define materials
# =============================================================================

# UO2 fuel — enriched to 4.5%
uo2 = model.add_material(1, name="UO2 fuel", density=10.97)
uo2.weight_fractions = False  # atom fractions
uo2.add_nuclide(92235, 0.045, "80c")   # U-235
uo2.add_nuclide(92238, 0.955, "80c")   # U-238
uo2.add_element(8, 2.0)                  # Oxygen (natural)

print(uo2)
print(f"  Nuclides: {uo2.nuclides}")
print(f"  Elements: {uo2.elements}")

# Zircaloy-4 cladding — by weight fraction
zirc = model.add_material(2, name="Zircaloy-4", density=6.56)
zirc.weight_fractions = True
zirc.add_element(40, 0.9823)   # Zr
zirc.add_element(50, 0.0145)   # Sn
zirc.add_element(26, 0.0021)   # Fe
zirc.add_element(24, 0.0011)   # Cr

print(f"\n{zirc}")
print(f"  Elements: {zirc.elements}")

# Expand Zircaloy elements to explicit nuclides
zirc.expand_elements()
print(f"  After expand: {len(zirc.nuclides)} nuclides")

# Light water — using element-based definition
water = model.add_material(3, name="Light water", density=0.997)
water.add_element(1, 2.0)    # Hydrogen
water.add_element(8, 1.0)    # Oxygen

print(f"\n{water}")

# Stainless steel 304
ss304 = model.add_material(4, name="SS304", density=7.94)
ss304.weight_fractions = True
ss304.add_element(26, 0.695)   # Fe
ss304.add_element(24, 0.190)   # Cr
ss304.add_element(28, 0.095)   # Ni
ss304.add_element(25, 0.020)   # Mn

print(f"\n{ss304}")

# =============================================================================
# Mixtures
# =============================================================================

# Homogenized fuel-clad mixture (e.g. for smeared pin regions)
# 85% UO2 + 15% Zircaloy-4 by weight
fuel_clad_id = model.create_mixture([1, 2], [0.85, 0.15],
                                     new_id=10, name="UO2-Zirc mix")
print(f"\nMixture M{fuel_clad_id}: 85% UO2 + 15% Zircaloy-4")

# Water + SS304 mix (just as demo)
water_ss_id = model.create_mixture([3, 4], [0.5, 0.5],
                                    new_id=20, name="Water-SS mix")
print(f"Mixture M{water_ss_id}: 50% Water + 50% SS304")

# =============================================================================
# Build geometry: simple fuel pin
# =============================================================================

# Surfaces
fuel_r  = ath.ZCylinder(0, 0, radius=0.4096)
clad_ir = ath.ZCylinder(0, 0, radius=0.418)
clad_or = ath.ZCylinder(0, 0, radius=0.475)
pitch   = ath.Box(-0.63, 0.63, -0.63, 0.63, -10, 10)

# Cells
model.add_cell(region=-fuel_r, material=1, density=10.97,
               name="fuel")
model.add_cell(region=+fuel_r & -clad_ir, material=0, density=0,
               name="gap")
model.add_cell(region=+clad_ir & -clad_or, material=2, density=6.56,
               name="clad")
model.add_cell(region=+clad_or & -pitch, material=3, density=0.997,
               name="water")

# =============================================================================
# Summary
# =============================================================================

print("\n" + "=" * 50)
print("Model summary")
print("=" * 50)
print(model)

print(f"\nMaterials ({len(model.materials)}):")
for mat in model.materials:
    nn = len(mat.nuclides)
    ne = len(mat.elements)
    d = mat.density
    print(f"  M{mat.id} ({mat.name}): density={d}, "
          f"{nn} nuclides, {ne} elements")

print(f"\nCells ({len(model.cells)}):")
for cell in model.cells:
    print(f"  Cell {cell.id} ({cell.name}): mat={cell.material}, "
          f"rho={cell.density}")

# =============================================================================
# Export
# =============================================================================

import tempfile, os
out = os.path.join(tempfile.gettempdir(), "fuel_pin.inp")
model.export_mcnp(out)
print(f"\nMCNP input written to: {out}")

# Show excerpt
with open(out) as f:
    lines = f.readlines()
print("\n" + "=" * 50)
print("Generated MCNP input (excerpt)")
print("=" * 50)
for line in lines[:40]:
    print(line, end="")
if len(lines) > 40:
    print(f"  ... ({len(lines) - 40} more lines)")
