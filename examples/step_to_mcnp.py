#!/usr/bin/env python3
"""
STEP to MCNP conversion example using OCP (OpenCascade) and aleathor.

This example demonstrates how to:
1. Load a STEP file using OCP's STEPControl_Reader
2. Analyze the B-Rep geometry to extract primitive surfaces
3. Map recognized shapes to aleathor CSG surfaces
4. Build an aleathor model and export to MCNP format

Requirements:
    pip install cadquery-ocp aleathor   (or any OCP distribution)

Note: Full B-Rep to CSG conversion is a hard problem. This example handles
common primitives (planes, cylinders, spheres, cones, tori) that are
axis-aligned or have analytically recognizable parameters. Complex freeform
surfaces (B-splines, NURBS) are not convertible to MCNP CSG and will be
skipped with a warning.
"""

import math
import sys
import tempfile
from pathlib import Path

from OCP.BRep import BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.GeomAbs import (
    GeomAbs_Plane,
    GeomAbs_Cylinder,
    GeomAbs_Sphere,
    GeomAbs_Cone,
    GeomAbs_Torus,
)
from OCP.STEPControl import STEPControl_Reader
from OCP.TopAbs import (
    TopAbs_FACE, TopAbs_SOLID, TopAbs_SHELL,
    TopAbs_FORWARD, TopAbs_REVERSED,
)
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS
from OCP.GProp import GProp_GProps
from OCP.BRepGProp import BRepGProp
from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib
from OCP.IFSelect import IFSelect_RetDone

import aleathor as ath
from aleathor.surfaces import (
    Plane,
    Sphere,
    CylinderX,
    CylinderY,
    CylinderZ,
    ConeX,
    ConeY,
    ConeZ,
    TorusX,
    TorusY,
    TorusZ,
    RCC,
    Box,
)


# ---------------------------------------------------------------------------
# Surface recognition helpers
# ---------------------------------------------------------------------------

AXIS_TOL = 1e-6  # tolerance for axis-alignment checks


def _is_axis(direction, ref):
    """Check if a direction (tuple or gp_Dir) is aligned with a reference axis."""
    if isinstance(direction, tuple):
        dx, dy, dz = direction
    else:
        dx, dy, dz = direction.X(), direction.Y(), direction.Z()
    return (abs(dx - ref[0]) < AXIS_TOL
            and abs(dy - ref[1]) < AXIS_TOL
            and abs(dz - ref[2]) < AXIS_TOL)


def _surface_info(face):
    """Extract analytic surface parameters from an OCC face.

    Returns a dict with 'type' and relevant geometric parameters,
    or None if the surface type is not supported.
    """
    adaptor = BRepAdaptor_Surface(TopoDS.Face_s(face))
    stype = adaptor.GetType()

    if stype == GeomAbs_Plane:
        pln = adaptor.Plane()
        loc = pln.Location()
        axis = pln.Axis().Direction()
        # Plane: a*x + b*y + c*z = d
        a, b, c = axis.X(), axis.Y(), axis.Z()
        d = a * loc.X() + b * loc.Y() + c * loc.Z()
        return {"type": "plane", "a": a, "b": b, "c": c, "d": d}

    elif stype == GeomAbs_Cylinder:
        cyl = adaptor.Cylinder()
        loc = cyl.Location()
        axis = cyl.Axis().Direction()
        r = cyl.Radius()
        return {
            "type": "cylinder",
            "origin": (loc.X(), loc.Y(), loc.Z()),
            "axis": (axis.X(), axis.Y(), axis.Z()),
            "radius": r,
        }

    elif stype == GeomAbs_Sphere:
        sph = adaptor.Sphere()
        loc = sph.Location()
        r = sph.Radius()
        return {
            "type": "sphere",
            "center": (loc.X(), loc.Y(), loc.Z()),
            "radius": r,
        }

    elif stype == GeomAbs_Cone:
        cone = adaptor.Cone()
        loc = cone.Location()
        axis = cone.Axis().Direction()
        half_angle = cone.SemiAngle()
        ref_radius = cone.RefRadius()
        return {
            "type": "cone",
            "apex": (loc.X(), loc.Y(), loc.Z()),
            "axis": (axis.X(), axis.Y(), axis.Z()),
            "half_angle": half_angle,
            "ref_radius": ref_radius,
        }

    elif stype == GeomAbs_Torus:
        tor = adaptor.Torus()
        loc = tor.Location()
        axis = tor.Axis().Direction()
        R = tor.MajorRadius()
        r = tor.MinorRadius()
        return {
            "type": "torus",
            "center": (loc.X(), loc.Y(), loc.Z()),
            "axis": (axis.X(), axis.Y(), axis.Z()),
            "major_radius": R,
            "minor_radius": r,
        }

    return None  # B-spline, NURBS, etc.


def _make_aleathor_surface(info, name=None):
    """Convert a surface info dict to an aleathor Surface object.

    Returns None if the surface geometry cannot be represented in MCNP CSG.
    """
    t = info["type"]

    if t == "plane":
        return Plane(info["a"], info["b"], info["c"], info["d"], name=name)

    elif t == "sphere":
        cx, cy, cz = info["center"]
        return Sphere(cx, cy, cz, info["radius"], name=name)

    elif t == "cylinder":
        ax = info["axis"]
        ox, oy, oz = info["origin"]
        r = info["radius"]
        if abs(abs(ax[2]) - 1) < AXIS_TOL:
            return CylinderZ(ox, oy, r, name=name)
        elif abs(abs(ax[1]) - 1) < AXIS_TOL:
            return CylinderY(ox, oz, r, name=name)
        elif abs(abs(ax[0]) - 1) < AXIS_TOL:
            return CylinderX(oy, oz, r, name=name)
        else:
            print(f"  WARNING: Off-axis cylinder (axis={ax}) not yet supported, skipping.")
            return None

    elif t == "cone":
        ax = info["axis"]
        apex = info["apex"]
        ha = info["half_angle"]
        t_sq = math.tan(ha) ** 2
        if abs(abs(ax[2]) - 1) < AXIS_TOL:
            return ConeZ(*apex, t_sq, name=name)
        elif abs(abs(ax[1]) - 1) < AXIS_TOL:
            return ConeY(*apex, t_sq, name=name)
        elif abs(abs(ax[0]) - 1) < AXIS_TOL:
            return ConeX(*apex, t_sq, name=name)
        else:
            print(f"  WARNING: Off-axis cone (axis={ax}) not yet supported, skipping.")
            return None

    elif t == "torus":
        ax = info["axis"]
        c = info["center"]
        R, r = info["major_radius"], info["minor_radius"]
        if abs(abs(ax[2]) - 1) < AXIS_TOL:
            return TorusZ(*c, R, r, name=name)
        elif abs(abs(ax[0]) - 1) < AXIS_TOL:
            return TorusX(*c, R, r, name=name)
        elif abs(abs(ax[1]) - 1) < AXIS_TOL:
            return TorusY(*c, R, r, name=name)
        else:
            print(f"  WARNING: Off-axis torus not yet supported, skipping.")
            return None

    return None


# ---------------------------------------------------------------------------
# Solid-level analysis
# ---------------------------------------------------------------------------

def _axis_aligned(ax):
    """Check if axis tuple is aligned with X, Y, or Z."""
    for ref in [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]:
        if all(abs(a - r) < AXIS_TOL for a, r in zip(ax, ref)):
            return True
    return False


def _classify_solid(solid):
    """Try to classify an OCC solid as a single MCNP macrobody.

    Recognizes:
    - Axis-aligned boxes (RPP)
    - Axis-aligned cylinders (RCC)
    - Spheres

    Returns (surface, info_str) or (None, reason).
    """
    faces = []
    explorer = TopExp_Explorer(solid, TopAbs_FACE)
    while explorer.More():
        info = _surface_info(explorer.Current())
        if info is not None:
            faces.append(info)
        else:
            return None, "contains non-analytic (B-spline) faces"
        explorer.Next()

    if not faces:
        return None, "no faces found"

    types = set(f["type"] for f in faces)

    # Pure sphere
    if types == {"sphere"} and len(faces) >= 1:
        s = faces[0]
        return Sphere(*s["center"], s["radius"]), "sphere"

    # Box: 6 planes
    if types == {"plane"} and len(faces) == 6:
        # Collect plane offsets per axis
        xs, ys, zs = [], [], []
        for f in faces:
            a, b, c, d = f["a"], f["b"], f["c"], f["d"]
            if abs(abs(a) - 1) < AXIS_TOL and abs(b) < AXIS_TOL and abs(c) < AXIS_TOL:
                xs.append(d * (1 if a > 0 else -1))
            elif abs(a) < AXIS_TOL and abs(abs(b) - 1) < AXIS_TOL and abs(c) < AXIS_TOL:
                ys.append(d * (1 if b > 0 else -1))
            elif abs(a) < AXIS_TOL and abs(b) < AXIS_TOL and abs(abs(c) - 1) < AXIS_TOL:
                zs.append(d * (1 if c > 0 else -1))
        if len(xs) == 2 and len(ys) == 2 and len(zs) == 2:
            return Box(min(xs), max(xs), min(ys), max(ys), min(zs), max(zs)), "box"

    # Cylinder: 2 planes + N cylinder faces (all same radius/axis)
    planes = [f for f in faces if f["type"] == "plane"]
    cyls = [f for f in faces if f["type"] == "cylinder"]
    if len(planes) == 2 and len(cyls) >= 1 and types <= {"plane", "cylinder"}:
        c0 = cyls[0]
        # All cylinder faces should share the same axis and radius
        if all(
            abs(c["radius"] - c0["radius"]) < AXIS_TOL
            and all(abs(a - b) < AXIS_TOL for a, b in zip(c["axis"], c0["axis"]))
            for c in cyls
        ):
            ax = c0["axis"]
            if _axis_aligned(ax):
                # Determine base and top from the two planes
                ox, oy, oz = c0["origin"]
                r = c0["radius"]
                # Project plane offsets onto cylinder axis
                dists = []
                for p in planes:
                    # d = a*x0 + b*y0 + c*z0 for point on plane along axis
                    dists.append(p["d"])
                # Base and height along the axis
                ax_abs = tuple(abs(v) for v in ax)
                if abs(ax_abs[2] - 1) < AXIS_TOL:  # Z-axis
                    z0, z1 = sorted(dists)
                    return RCC(ox, oy, z0, 0, 0, z1 - z0, r), "rcc_z"
                elif abs(ax_abs[1] - 1) < AXIS_TOL:  # Y-axis
                    y0, y1 = sorted(dists)
                    return RCC(ox, y0, oz, 0, y1 - y0, 0, r), "rcc_y"
                elif abs(ax_abs[0] - 1) < AXIS_TOL:  # X-axis
                    x0, x1 = sorted(dists)
                    return RCC(x0, oy, oz, x1 - x0, 0, 0, r), "rcc_x"

    return None, f"complex shape with face types: {types}"


# ---------------------------------------------------------------------------
# High-level conversion
# ---------------------------------------------------------------------------

def step_to_model(step_path, title="STEP Import", default_material=1,
                  default_density=1.0):
    """Convert a STEP file to an aleathor Model.

    Each solid in the STEP file becomes a cell in the model. Recognized
    primitive shapes (boxes, cylinders, spheres) are converted to MCNP
    macrobodies. More complex solids are decomposed face-by-face into
    bounding surfaces.

    Args:
        step_path: Path to the STEP file.
        title: Model title for the MCNP card.
        default_material: Material number to assign to all cells.
        default_density: Density (g/cm3) to assign to all cells.

    Returns:
        (aleathor.Model, bounds) where bounds is
        (xmin, xmax, ymin, ymax, zmin, zmax) with margin.
    """
    # Load the STEP file using OCP directly
    print(f"Loading STEP file: {step_path}")
    reader = STEPControl_Reader()
    status = reader.ReadFile(str(step_path))
    if status != IFSelect_RetDone:
        raise RuntimeError(f"Failed to read STEP file: {step_path}")
    reader.TransferRoots()
    shape = reader.OneShape()

    # Extract all solids
    solids = []
    explorer = TopExp_Explorer(shape, TopAbs_SOLID)
    while explorer.More():
        solids.append(explorer.Current())
        explorer.Next()
    print(f"Found {len(solids)} solid(s) in STEP file\n")

    # Compute overall bounding box with margin for void generation
    bbox = Bnd_Box()
    BRepBndLib.Add_s(shape, bbox)
    xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
    margin = max(xmax - xmin, ymax - ymin, zmax - zmin) * 0.1
    model_bounds = (
        xmin - margin, xmax + margin,
        ymin - margin, ymax + margin,
        zmin - margin, zmax + margin,
    )

    model = ath.Model(title)

    for i, solid in enumerate(solids):
        solid_name = f"solid_{i+1}"

        # Compute volume for info
        props = GProp_GProps()
        BRepGProp.VolumeProperties_s(solid, props)
        volume = props.Mass()

        print(f"Solid {i+1}: volume = {volume:.4f} cm3")

        # Try macrobody classification first
        macro_surf, classification = _classify_solid(solid)
        if macro_surf is not None:
            print(f"  -> Recognized as: {classification}")
            model.add_cell(
                region=-macro_surf,
                material=default_material,
                density=default_density,
                name=solid_name,
            )
            continue

        # Fall back to shell-by-shell decomposition.
        # A solid is bounded by shells:
        #   FORWARD shell  -> outer boundary (material is inside)
        #   REVERSED shell -> cavity/hole (material is outside)
        # Each face within a shell has an orientation:
        #   FORWARD  -> normal outward -> material on negative side
        #   REVERSED -> normal inward  -> material on positive side
        # The solid region = (outer shell interior) AND NOT (cavity interiors)
        print(f"  -> Decomposing by shells ({classification})")
        outer_regions = []
        cavity_regions = []
        unsupported = 0
        face_idx = 0

        shell_exp = TopExp_Explorer(solid, TopAbs_SHELL)
        while shell_exp.More():
            shell = shell_exp.Current()
            is_cavity = (shell.Orientation() == TopAbs_REVERSED)

            halfspaces = []
            face_exp = TopExp_Explorer(shell, TopAbs_FACE)
            while face_exp.More():
                face = face_exp.Current()
                orientation = face.Orientation()
                info = _surface_info(face)
                if info is not None:
                    surf = _make_aleathor_surface(
                        info, name=f"{solid_name}_f{face_idx}"
                    )
                    if surf is not None:
                        # FORWARD face -> normal outward -> interior on negative side
                        # REVERSED face -> normal inward -> interior on positive side
                        # For a REVERSED shell (cavity), the sense is flipped:
                        # the shell defines the cavity interior, not the solid.
                        use_negative = (orientation == TopAbs_FORWARD)
                        if is_cavity:
                            use_negative = not use_negative
                        halfspaces.append((surf, use_negative))
                    else:
                        unsupported += 1
                else:
                    unsupported += 1
                    print(f"  WARNING: Face {face_idx} is non-analytic, skipping")
                face_idx += 1
                face_exp.Next()

            if halfspaces:
                # Build shell interior as intersection of halfspaces
                surf, neg = halfspaces[0]
                shell_region = -surf if neg else +surf
                for surf, neg in halfspaces[1:]:
                    shell_region = shell_region & (-surf if neg else +surf)

                if is_cavity:
                    cavity_regions.append(shell_region)
                else:
                    outer_regions.append(shell_region)

            shell_exp.Next()

        if unsupported > 0:
            print(f"  WARNING: {unsupported} face(s) could not be converted")

        if outer_regions:
            # Start with intersection of all outer shell regions
            region = outer_regions[0]
            for r in outer_regions[1:]:
                region = region & r
            # Subtract each cavity
            for cavity in cavity_regions:
                region = region & ~cavity
            model.add_cell(
                region=region,
                material=default_material,
                density=default_density,
                name=solid_name,
            )
            print(f"  -> Created cell ({len(outer_regions)} outer shell(s), "
                  f"{len(cavity_regions)} cavity/ies)")
        else:
            print(f"  WARNING: No convertible surfaces for solid {i+1}")

    print(f"\nModel created: {len(model.cells)} cell(s)")
    return model, model_bounds


# ---------------------------------------------------------------------------
# Demo: create a sample STEP file and convert it
# ---------------------------------------------------------------------------

def create_sample_step(path):
    """Create a sample STEP file with a shielding assembly.

    The geometry is a simplified neutron shielding plug:
    - A steel cylindrical shell (hollow cylinder)
    - A concrete filling block inside
    - A lead sphere embedded in the concrete
    """
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakeSphere
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
    from OCP.gp import gp_Pnt, gp_Dir, gp_Ax2
    from OCP.TopoDS import TopoDS_Compound, TopoDS_Builder
    from OCP.STEPControl import STEPControl_Writer, STEPControl_AsIs

    ax = gp_Ax2(gp_Pnt(0, 0, -10), gp_Dir(0, 0, 1))

    # Outer steel cylinder (h=20, r=8)
    outer_cyl = BRepPrimAPI_MakeCylinder(ax, 8, 20).Shape()
    # Inner cut (h=18, r=6, offset 1cm from base)
    ax_inner = gp_Ax2(gp_Pnt(0, 0, -9), gp_Dir(0, 0, 1))
    inner_cyl = BRepPrimAPI_MakeCylinder(ax_inner, 6, 18).Shape()
    steel_shell = BRepAlgoAPI_Cut(outer_cyl, inner_cyl).Shape()

    # Concrete filling (h=18, r=5.9)
    concrete_cyl = BRepPrimAPI_MakeCylinder(ax_inner, 5.9, 18).Shape()

    # Lead sphere at center (r=3)
    lead_sphere = BRepPrimAPI_MakeSphere(gp_Pnt(0, 0, 0), 3).Shape()

    # Concrete with sphere cavity
    concrete_with_cavity = BRepAlgoAPI_Cut(concrete_cyl, lead_sphere).Shape()

    # Assemble into a compound
    compound = TopoDS_Compound()
    builder = TopoDS_Builder()
    builder.MakeCompound(compound)
    builder.Add(compound, steel_shell)
    builder.Add(compound, concrete_with_cavity)
    builder.Add(compound, lead_sphere)

    # Write STEP
    writer = STEPControl_Writer()
    writer.Transfer(compound, STEPControl_AsIs)
    writer.Write(str(path))
    print(f"Sample STEP file written to: {path}\n")


def main():
    if len(sys.argv) > 1:
        step_path = Path(sys.argv[1])
        if not step_path.exists():
            print(f"Error: file not found: {step_path}")
            sys.exit(1)
    else:
        # Create a sample STEP file for demonstration
        step_path = Path(tempfile.gettempdir()) / "shielding_sample.step"
        create_sample_step(step_path)

    # ---- Convert STEP -> aleathor Model ----
    model, bounds = step_to_model(
        step_path,
        title="Shielding Assembly from STEP",
        default_material=1,
        default_density=7.8,  # steel density as default
    )

    # Assign proper materials and densities
    for cell in model.cells:
        if "steel" in cell.name:
            cell.material = 1   # steel
            cell.density = 7.85
        elif "concrete" in cell.name:
            cell.material = 2   # concrete
            cell.density = 2.35
        elif "lead" in cell.name:
            cell.material = 3   # lead
            cell.density = 11.34

    # ---- Generate void cells and graveyard ----
    print("\n" + "=" * 50)
    print("Void generation")
    print("=" * 50)
    voids = model.generate_void(bounds=bounds, max_depth=3)
    print(f"Found {len(voids)} void boxes")
    voids.merge()
    print(f"After merge: {len(voids)} boxes")
    n_added = voids.add_cells()
    print(f"Added {n_added} void cell(s)")
    voids.add_graveyard()
    print("Added graveyard cell")

    # ---- Simplify CSG (eliminate #, double negation, etc.) ----
    stats = model.simplify()
    print(f"\nCSG simplification: {stats['nodes_before']} -> {stats['nodes_after']} nodes")

    # ---- Inspect the model ----
    print("\n" + "=" * 50)
    print("Model summary")
    print("=" * 50)
    for cell in model.cells:
        mat_str = f"mat={cell.material}, rho={cell.density} g/cm3" if cell.material else "void"
        print(f"  Cell {cell.id} ({cell.name}): {mat_str}")

    # ---- Export to MCNP ----
    output = Path(tempfile.gettempdir()) / "shielding_from_step.inp"
    model.save(str(output))
    print(f"\nMCNP input written to: {output}")

    # Show the first part of the generated file
    print("\n" + "=" * 50)
    print("Generated MCNP input (excerpt)")
    print("=" * 50)
    with open(output) as f:
        lines = f.readlines()
    for line in lines[:60]:
        print(line, end="")
    if len(lines) > 60:
        print(f"  ... ({len(lines) - 60} more lines)")


if __name__ == "__main__":
    main()
