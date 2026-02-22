# MCNP Test Geometries

This directory contains MCNP input files for testing the aleathor library.

## Simple Geometries

| File | Description |
|------|-------------|
| `sphere_simple.inp` | Single sphere at origin |
| `sphere_offcenter.inp` | Sphere not at origin (tests S surface) |
| `concentric_spheres.inp` | Multiple concentric spherical shells |
| `hemisphere.inp` | Sphere cut by plane |
| `void_only.inp` | Minimal void-only geometry |

## Cylinders

| File | Description |
|------|-------------|
| `cylinder_z.inp` | Vertical cylinder with caps |
| `cylinder_x.inp` | Horizontal cylinder along X-axis |
| `multi_cylinder.inp` | Multiple parallel cylinders |

## Other Surfaces

| File | Description |
|------|-------------|
| `box_rpp.inp` | RPP macrobody (rectangular box) |
| `cone_z.inp` | Conical surface (KZ) |
| `torus.inp` | Torus surface (TZ) |
| `ellipsoid.inp` | Ellipsoid using GQ surface |

## Macrobodies

| File | Description |
|------|-------------|
| `rcc_macrobody.inp` | Right Circular Cylinder |
| `trc_macrobody.inp` | Truncated Right Cone |

## Boolean Operations

| File | Description |
|------|-------------|
| `boolean_complex.inp` | Complex CSG (intersection, union, complement) |

## Universe/FILL

| File | Description |
|------|-------------|
| `universe_fill.inp` | Simple universe with FILL |
| `nested_universes.inp` | Three-level universe hierarchy |


## Realistic Configurations

| File | Description |
|------|-------------|
| `pin_cell.inp` | PWR fuel pin cell |
| `reactor_core.inp` | Simplified PWR core cross-section |
| `shielding.inp` | Multi-layer shielding configuration |

| `spallation_target.inp` | SNS-like spallation target |

## Usage

```python
import aleathor as ath

# Load any test file
model = ath.load_mcnp("tests/data/concentric_spheres.inp")

# Plot
model.plot(z=0, bounds=(-12, 12, -12, 12))
```

Or from command line:
```bash
python examples/plot_geometry.py tests/data/concentric_spheres.inp
```
