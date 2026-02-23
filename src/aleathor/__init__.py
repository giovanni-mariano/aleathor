# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""
AleaTHOR: Constructive Solid Geometry for Nuclear Simulations

A Python package for building, manipulating, and exporting CSG geometries
compatible with MCNP and OpenMC.

Example:
    import aleathor as ath

    # Create geometry
    fuel = ath.Sphere(0, 0, 0, radius=5.0).interior()
    clad = ath.CylinderZ(0, 0, radius=5.5).interior() - fuel
    water = ath.Box(-10, 10, -10, 10, -10, 10).interior() - fuel - clad

    # Create model
    model = ath.Model()
    model.add_cell(fuel, material=1, density=-10.5, name="fuel")
    model.add_cell(clad, material=2, density=-6.5, name="clad")
    model.add_cell(water, material=3, density=-1.0, name="water")

    # Export
    model.export_mcnp("geometry.inp")
    model.export_openmc("geometry.xml")
"""

__version__ = "0.1.0"

from .geometry import (
    Region,
    Halfspace,
    Intersection,
    Union,
    Complement,
)

from .surfaces import (
    Surface,
    Plane,
    XPlane,
    YPlane,
    ZPlane,
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
    Box,
    RCC,
    TRC,
    ELL,
    REC,
    WED,
    RHP,
    GeneralBox,
    Quadric,
)

from .model import (
    Model,
    Material,
    Universe,
)

from .io import (
    read_mcnp,
    read_mcnp_string,
    read_openmc,
    write_mcnp,
    write_openmc,
)

from .collections import (
    Cell,
    CellCollection,
    TraceResult,
    TraceSegment,
)

# Design document API aliases
load_mcnp = read_mcnp
load_openmc = read_openmc


def load(filename: str) -> Model:
    """Load a geometry file, auto-detecting format.

    Supported formats:
    - MCNP input files (.inp, .i, .mcnp, or no extension)
    - OpenMC XML files (.xml)

    Args:
        filename: Path to geometry file

    Returns:
        Model with loaded geometry

    Example:
        model = ath.load("input.inp")     # MCNP
        model = ath.load("geometry.xml")  # OpenMC
    """
    from pathlib import Path
    path = Path(filename)
    ext = path.suffix.lower()

    if ext == '.xml':
        return read_openmc(filename)
    else:
        # Default to MCNP
        return read_mcnp(filename)

# Logging integration (routes C library logs to Python logging module)
try:
    import _alea
    # Log level constants
    LOG_NONE = _alea.LOG_NONE
    LOG_ERROR = _alea.LOG_ERROR
    LOG_WARN = _alea.LOG_WARN
    LOG_INFO = _alea.LOG_INFO
    LOG_DEBUG = _alea.LOG_DEBUG
    LOG_TRACE = _alea.LOG_TRACE
    # Functions
    set_log_level = _alea.set_log_level
    get_log_level = _alea.get_log_level

    def enable_logging():
        """Enable logging at INFO level."""
        _alea.set_log_level(_alea.LOG_INFO)

    def disable_logging():
        """Disable logging."""
        _alea.set_log_level(_alea.LOG_NONE)

    HAS_NATIVE = True
except ImportError:
    HAS_NATIVE = False
    LOG_NONE = LOG_ERROR = LOG_WARN = LOG_INFO = LOG_DEBUG = LOG_TRACE = 0
    def set_log_level(level): pass
    def get_log_level(): return 0
    def enable_logging(): pass
    def disable_logging(): pass

# Optional plotting (requires matplotlib)
try:
    from .plotting import (
        plot_slice_curves,
        plot_slice_filled,
        plot_curve,
        plot_ray_path,
        count_surfaces,
        count_cells,
        count_materials,
        get_slice_stats,
        filter_boundary_curves,
        filter_boundary_curves_arbitrary,
        CELL_VOID,
        CELL_UNDEFINED,
        CELL_OVERLAP,
        GRID_ERROR_OK,
        GRID_ERROR_OVERLAP,
        GRID_ERROR_UNDEFINED,
    )
    HAS_PLOTTING = True
except ImportError:
    HAS_PLOTTING = False

__all__ = [
    # Geometry
    'Region',
    'Halfspace',
    'Intersection',
    'Union',
    'Complement',
    # Surfaces
    'Surface',
    'Plane',
    'XPlane',
    'YPlane',
    'ZPlane',
    'Sphere',
    'CylinderX',
    'CylinderY',
    'CylinderZ',
    'ConeX',
    'ConeY',
    'ConeZ',
    'TorusX',
    'TorusY',
    'TorusZ',
    'Box',
    'RCC',
    'TRC',
    'ELL',
    'REC',
    'WED',
    'RHP',
    'GeneralBox',
    'Quadric',
    # Model
    'Model',
    'Cell',
    'Material',
    'Universe',
    # I/O
    'load',
    'load_mcnp',
    'load_openmc',
    'read_mcnp',
    'read_mcnp_string',
    'read_openmc',
    'write_mcnp',
    'write_openmc',
    # Collections
    'Cell',
    'CellCollection',
    'TraceResult',
    'TraceSegment',
    # Logging
    'LOG_NONE',
    'LOG_ERROR',
    'LOG_WARN',
    'LOG_INFO',
    'LOG_DEBUG',
    'LOG_TRACE',
    'set_log_level',
    'get_log_level',
    'enable_logging',
    'disable_logging',
    # Plotting (optional, analytical)
    'plot_slice_curves',
    'plot_slice_filled',
    'plot_curve',
    'plot_ray_path',
    # Counting functions
    'count_surfaces',
    'count_cells',
    'count_materials',
    'get_slice_stats',
    # Curve filtering
    'filter_boundary_curves',
    'filter_boundary_curves_arbitrary',
    'CELL_VOID',
    'CELL_UNDEFINED',
    'CELL_OVERLAP',
    # Grid error constants
    'GRID_ERROR_OK',
    'GRID_ERROR_OVERLAP',
    'GRID_ERROR_UNDEFINED',
]
