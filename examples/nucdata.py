#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""
Nuclear data example.

Demonstrates how to:
1. Load a cross-section directory (xsdir)
2. Query nuclide properties and cross sections
3. Build a material and compute macroscopic quantities
4. Collapse to multigroup constants
5. Doppler-broaden cross sections

Requires ACE-format nuclear data files and an xsdir index.
"""

import aleathor as ath
from aleathor.nucdata import XsDir, Nuclide, NucMaterial, Multigroup, parse_zaid

# =============================================================================
# 1. Load cross-section directory
# =============================================================================

print("=" * 60)
print("Loading cross-section directory")
print("=" * 60)

# Point this to your xsdir file
xsdir = XsDir("/path/to/xsdir")
print(f"\n{xsdir}")

# =============================================================================
# 2. Load nuclides and inspect properties
# =============================================================================

print("\n" + "=" * 60)
print("Nuclide properties")
print("=" * 60)

u235 = Nuclide(xsdir, "92235.80c")
u238 = Nuclide(xsdir, "92238.80c")
o16 = Nuclide(xsdir, "8016.80c")
h1 = Nuclide(xsdir, "1001.80c")

for nuc in [u235, u238, o16, h1]:
    print(f"\n  {nuc}")
    print(f"    AWR: {nuc.awr:.4f}")
    print(f"    Temperature: {nuc.temperature:.2e} MeV")
    print(f"    Fissile: {nuc.is_fissile}")
    print(f"    Energy grid points: {nuc.n_energies}")
    print(f"    Reactions: {nuc.n_reactions}")

# =============================================================================
# 3. Cross-section lookup
# =============================================================================

print("\n" + "=" * 60)
print("Cross-section lookup")
print("=" * 60)

energies = [1e-8, 0.0253e-6, 1e-3, 1.0, 14.0]  # MeV

print(f"\nU-235 cross sections (barns):")
print(f"  {'Energy (MeV)':>14s}  {'Total':>10s}  {'Absorption':>10s}  {'Elastic':>10s}")
for e in energies:
    st = u235.xs_total(e)
    sa = u235.xs_absorption(e)
    se = u235.xs_elastic(e)
    print(f"  {e:14.4e}  {st:10.3f}  {sa:10.3f}  {se:10.3f}")

# Fission-specific data
print(f"\nU-235 nu-bar (avg. neutrons per fission):")
for e in [0.0253e-6, 1.0, 14.0]:
    print(f"  {e:.4e} MeV -> nu = {u235.nu_bar(e):.4f}")

# Reaction list
print(f"\nU-235 reactions:")
for rxn in u235.reactions()[:10]:
    print(f"  MT={rxn['mt']:3d}  Q={rxn['q_value']:+8.4f} MeV")

# =============================================================================
# 4. Materials and macroscopic cross sections
# =============================================================================

print("\n" + "=" * 60)
print("Material (UO2 fuel, 4.5% enriched)")
print("=" * 60)

# Number densities for UO2 at 10.97 g/cm3, 4.5% enriched
mat = NucMaterial()
mat.add(u235, 0.001105)   # atoms/barn-cm
mat.add(u238, 0.02234)
mat.add(o16, 0.04689)

print(f"\nMacroscopic cross sections at 1 MeV:")
print(f"  Sigma_t   = {mat.xs_total(1.0):.4f} cm^-1")
print(f"  Sigma_a   = {mat.xs_absorption(1.0):.4f} cm^-1")
print(f"  Sigma_el  = {mat.xs_elastic(1.0):.4f} cm^-1")
print(f"  MFP       = {mat.mean_free_path(1.0):.4f} cm")

print(f"\nMacroscopic cross sections at thermal (0.0253 eV):")
e_th = 0.0253e-6
print(f"  Sigma_t   = {mat.xs_total(e_th):.4f} cm^-1")
print(f"  Sigma_a   = {mat.xs_absorption(e_th):.4f} cm^-1")
print(f"  MFP       = {mat.mean_free_path(e_th):.4f} cm")

# Monte Carlo sampling demo
print(f"\nMonte Carlo sampling (5 random collisions at 1 MeV):")
import random
for i in range(5):
    dist = mat.sample_distance(1.0, random.random())
    idx, zaid = mat.sample_nuclide(1.0, random.random())
    print(f"  Collision {i+1}: distance={dist:.4f} cm, hit {zaid}")

# =============================================================================
# 5. Multigroup collapse
# =============================================================================

print("\n" + "=" * 60)
print("Multigroup collapse (U-235)")
print("=" * 60)

# CASMO-4 style energy group boundaries (MeV, descending)
bounds = [20.0, 6.065, 0.821, 0.00553, 6.25e-7, 1e-11]
mg = Multigroup(bounds)
mg.collapse(u235)

data = mg.get_data()
print(f"\n{mg}")
print(f"\n  {'Group':>5s}  {'E_high (MeV)':>12s}  {'E_low (MeV)':>12s}  "
      f"{'sigma_t':>10s}  {'sigma_a':>10s}  {'nu*sigma_f':>10s}  {'chi':>8s}")
for g in range(data['n_groups']):
    print(f"  {g+1:5d}  {bounds[g]:12.4e}  {bounds[g+1]:12.4e}  "
          f"{data['sigma_t'][g]:10.3f}  {data['sigma_a'][g]:10.3f}  "
          f"{data['nu_sigma_f'][g]:10.3f}  {data['chi'][g]:8.5f}")

# Scattering matrix
print(f"\n  Scattering matrix (sigma_s g->g'):")
ng = data['n_groups']
header = "  " + " " * 8 + "".join(f"  g'={g+1:2d}  " for g in range(ng))
print(header)
for g_from in range(ng):
    row = f"  g={g_from+1:2d}  "
    for g_to in range(ng):
        row += f"  {mg.scatter(g_from, g_to):8.3f}"
    print(row)

# =============================================================================
# 6. Doppler broadening
# =============================================================================

print("\n" + "=" * 60)
print("Doppler broadening")
print("=" * 60)

# Load a fresh (non-cached) copy for broadening
u235_cold = Nuclide(xsdir, "92235.80c", cached=False)
print(f"\n  Original temperature: {u235_cold.temperature:.4e} MeV")

# Broaden to 600 K (kT = 5.17e-8 MeV)
kT_600K = 5.17e-8
u235_cold.doppler_broaden(kT_600K)
print(f"  After broadening:    {u235_cold.temperature:.4e} MeV")

# Compare XS at a resonance energy
e_res = 6.67e-6  # 6.67 eV resonance of U-235
xs_orig = u235.xs_total(e_res)
xs_broad = u235_cold.xs_total(e_res)
print(f"\n  Total XS at {e_res:.2e} MeV (near 6.67 eV resonance):")
print(f"    Original:  {xs_orig:.1f} barns")
print(f"    Broadened: {xs_broad:.1f} barns")

# =============================================================================
# 7. ZAID parsing utility
# =============================================================================

print("\n" + "=" * 60)
print("ZAID parsing")
print("=" * 60)

for zaid in ["92235.80c", "94239.70c", "1001.80c", "8016.80c"]:
    info = parse_zaid(zaid)
    print(f"  {zaid:>12s} -> Z={info['Z']:3d}, A={info['A']:3d}, type={info['type']}")

print("\n" + "=" * 60)
print("Example complete!")
print("=" * 60)
