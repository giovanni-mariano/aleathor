# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""
Surface definitions for CSG geometry.

Surfaces divide space into two halfspaces (positive and negative).
Use surface.positive() / surface.negative() or +surface / -surface
to get regions.
"""

from __future__ import annotations
from abc import ABC, abstractmethod
from typing import Tuple, Optional
import math


class Surface(ABC):
    """Abstract base class for surfaces.

    A surface divides 3D space into two regions (halfspaces).
    The sign convention is:
        - Positive halfspace: f(x,y,z) > 0
        - Negative halfspace: f(x,y,z) < 0

    For closed surfaces (sphere, cylinder, cone):
        - Negative = inside
        - Positive = outside

    Attributes:
        id: Surface ID (assigned when added to model).
        name: Optional human-readable name.
        boundary: Boundary condition ('transmissive', 'reflective', 'vacuum').
    """

    _next_id = 1  # Auto-increment ID counter

    def __init__(self, name: Optional[str] = None,
                 boundary: str = 'transmissive',
                 surface_id: Optional[int] = None):
        """
        Args:
            name: Optional name for the surface.
            boundary: Boundary condition type.
            surface_id: Optional explicit surface ID.
        """
        self.name = name
        self.boundary = boundary

        if surface_id is not None:
            self.id = surface_id
        else:
            self.id = Surface._next_id
            Surface._next_id += 1

    @abstractmethod
    def evaluate(self, point: Tuple[float, float, float]) -> float:
        """Evaluate surface equation at point.

        Returns positive value if point is on positive side,
        negative if on negative side, zero if on surface.
        """
        pass

    @abstractmethod
    def _get_type(self) -> str:
        """Get surface type string for export."""
        pass

    @abstractmethod
    def _get_params(self) -> Tuple:
        """Get surface parameters for export."""
        pass

    def positive(self) -> 'Halfspace':
        """Get positive halfspace (f(x,y,z) > 0)."""
        from .geometry import Halfspace
        return Halfspace(self, positive=True)

    def negative(self) -> 'Halfspace':
        """Get negative halfspace (f(x,y,z) < 0)."""
        from .geometry import Halfspace
        return Halfspace(self, positive=False)

    def interior(self) -> 'Halfspace':
        """Get interior (negative halfspace). Alias for negative()."""
        return self.negative()

    def exterior(self) -> 'Halfspace':
        """Get exterior (positive halfspace). Alias for positive()."""
        return self.positive()

    def __pos__(self) -> 'Halfspace':
        """Positive halfspace: +surface"""
        return self.positive()

    def __neg__(self) -> 'Halfspace':
        """Negative halfspace: -surface"""
        return self.negative()

    def __repr__(self) -> str:
        if self.name:
            return f"{self.__class__.__name__}({self.id}, name='{self.name}')"
        return f"{self.__class__.__name__}({self.id})"


class Plane(Surface):
    """General plane: ax + by + cz = d

    The sign convention is: f(x,y,z) = ax + by + cz - d
        - Positive halfspace: ax + by + cz > d (in direction of normal)
        - Negative halfspace: ax + by + cz < d (opposite to normal)
    """

    def __init__(self, a: float, b: float, c: float, d: float,
                 name: Optional[str] = None, **kwargs):
        """
        Args:
            a, b, c: Normal vector components.
            d: Distance from origin (ax + by + cz = d).
            name: Optional surface name.
        """
        super().__init__(name=name, **kwargs)
        # Normalize
        norm = math.sqrt(a*a + b*b + c*c)
        if norm < 1e-10:
            raise ValueError("Plane normal vector cannot be zero")
        self.a = a / norm
        self.b = b / norm
        self.c = c / norm
        self.d = d / norm

    def evaluate(self, point: Tuple[float, float, float]) -> float:
        x, y, z = point
        return self.a * x + self.b * y + self.c * z - self.d

    def _get_type(self) -> str:
        return 'P'

    def _get_params(self) -> Tuple:
        return (self.a, self.b, self.c, self.d)


class XPlane(Plane):
    """Plane perpendicular to X axis: x = x0"""

    def __init__(self, x0: float, name: Optional[str] = None, **kwargs):
        """
        Args:
            x0: X coordinate of plane.
        """
        super().__init__(1, 0, 0, x0, name=name, **kwargs)
        self.x0 = x0

    def _get_type(self) -> str:
        return 'PX'

    def _get_params(self) -> Tuple:
        return (self.x0,)


class YPlane(Plane):
    """Plane perpendicular to Y axis: y = y0"""

    def __init__(self, y0: float, name: Optional[str] = None, **kwargs):
        """
        Args:
            y0: Y coordinate of plane.
        """
        super().__init__(0, 1, 0, y0, name=name, **kwargs)
        self.y0 = y0

    def _get_type(self) -> str:
        return 'PY'

    def _get_params(self) -> Tuple:
        return (self.y0,)


class ZPlane(Plane):
    """Plane perpendicular to Z axis: z = z0"""

    def __init__(self, z0: float, name: Optional[str] = None, **kwargs):
        """
        Args:
            z0: Z coordinate of plane.
        """
        super().__init__(0, 0, 1, z0, name=name, **kwargs)
        self.z0 = z0

    def _get_type(self) -> str:
        return 'PZ'

    def _get_params(self) -> Tuple:
        return (self.z0,)


class Sphere(Surface):
    """Sphere: (x-x0)² + (y-y0)² + (z-z0)² = R²

    Sign convention:
        - Negative: inside sphere (r < R)
        - Positive: outside sphere (r > R)
    """

    def __init__(self, x0: float, y0: float, z0: float, radius: float,
                 name: Optional[str] = None, **kwargs):
        """
        Args:
            x0, y0, z0: Center coordinates.
            radius: Sphere radius.
        """
        super().__init__(name=name, **kwargs)
        if radius <= 0:
            raise ValueError("Sphere radius must be positive")
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.radius = radius

    def evaluate(self, point: Tuple[float, float, float]) -> float:
        x, y, z = point
        dx, dy, dz = x - self.x0, y - self.y0, z - self.z0
        return dx*dx + dy*dy + dz*dz - self.radius*self.radius

    def _get_type(self) -> str:
        if abs(self.x0) < 1e-10 and abs(self.y0) < 1e-10 and abs(self.z0) < 1e-10:
            return 'SO'  # Sphere centered at origin
        return 'S'

    def _get_params(self) -> Tuple:
        if self._get_type() == 'SO':
            return (self.radius,)
        return (self.x0, self.y0, self.z0, self.radius)


class CylinderX(Surface):
    """Infinite cylinder along X axis: (y-y0)² + (z-z0)² = R²"""

    def __init__(self, y0: float, z0: float, radius: float,
                 name: Optional[str] = None, **kwargs):
        """
        Args:
            y0, z0: Center in YZ plane.
            radius: Cylinder radius.
        """
        super().__init__(name=name, **kwargs)
        if radius <= 0:
            raise ValueError("Cylinder radius must be positive")
        self.y0 = y0
        self.z0 = z0
        self.radius = radius

    def evaluate(self, point: Tuple[float, float, float]) -> float:
        x, y, z = point
        dy, dz = y - self.y0, z - self.z0
        return dy*dy + dz*dz - self.radius*self.radius

    def _get_type(self) -> str:
        if abs(self.y0) < 1e-10 and abs(self.z0) < 1e-10:
            return 'CX'
        return 'C/X'

    def _get_params(self) -> Tuple:
        if self._get_type() == 'CX':
            return (self.radius,)
        return (self.y0, self.z0, self.radius)


class CylinderY(Surface):
    """Infinite cylinder along Y axis: (x-x0)² + (z-z0)² = R²"""

    def __init__(self, x0: float, z0: float, radius: float,
                 name: Optional[str] = None, **kwargs):
        """
        Args:
            x0, z0: Center in XZ plane.
            radius: Cylinder radius.
        """
        super().__init__(name=name, **kwargs)
        if radius <= 0:
            raise ValueError("Cylinder radius must be positive")
        self.x0 = x0
        self.z0 = z0
        self.radius = radius

    def evaluate(self, point: Tuple[float, float, float]) -> float:
        x, y, z = point
        dx, dz = x - self.x0, z - self.z0
        return dx*dx + dz*dz - self.radius*self.radius

    def _get_type(self) -> str:
        if abs(self.x0) < 1e-10 and abs(self.z0) < 1e-10:
            return 'CY'
        return 'C/Y'

    def _get_params(self) -> Tuple:
        if self._get_type() == 'CY':
            return (self.radius,)
        return (self.x0, self.z0, self.radius)


class CylinderZ(Surface):
    """Infinite cylinder along Z axis: (x-x0)² + (y-y0)² = R²"""

    def __init__(self, x0: float, y0: float, radius: float,
                 name: Optional[str] = None, **kwargs):
        """
        Args:
            x0, y0: Center in XY plane.
            radius: Cylinder radius.
        """
        super().__init__(name=name, **kwargs)
        if radius <= 0:
            raise ValueError("Cylinder radius must be positive")
        self.x0 = x0
        self.y0 = y0
        self.radius = radius

    def evaluate(self, point: Tuple[float, float, float]) -> float:
        x, y, z = point
        dx, dy = x - self.x0, y - self.y0
        return dx*dx + dy*dy - self.radius*self.radius

    def _get_type(self) -> str:
        if abs(self.x0) < 1e-10 and abs(self.y0) < 1e-10:
            return 'CZ'
        return 'C/Z'

    def _get_params(self) -> Tuple:
        if self._get_type() == 'CZ':
            return (self.radius,)
        return (self.x0, self.y0, self.radius)


class ConeX(Surface):
    """Cone along X axis: (y-y0)²/t² + (z-z0)²/t² = (x-x0)²"""

    def __init__(self, x0: float, y0: float, z0: float, t_sq: float,
                 sheet: int = 0, name: Optional[str] = None, **kwargs):
        """
        Args:
            x0, y0, z0: Apex coordinates.
            t_sq: tan²(half-angle).
            sheet: Sheet selection (+1, -1, or 0 for both).
        """
        super().__init__(name=name, **kwargs)
        self.x0, self.y0, self.z0 = x0, y0, z0
        self.t_sq = t_sq
        self.sheet = sheet

    def evaluate(self, point: Tuple[float, float, float]) -> float:
        x, y, z = point
        dx, dy, dz = x - self.x0, y - self.y0, z - self.z0
        return dy*dy + dz*dz - self.t_sq * dx*dx

    def _get_type(self) -> str:
        return 'KX'

    def _get_params(self) -> Tuple:
        return (self.x0, self.t_sq, self.sheet)


class ConeY(Surface):
    """Cone along Y axis."""

    def __init__(self, x0: float, y0: float, z0: float, t_sq: float,
                 sheet: int = 0, name: Optional[str] = None, **kwargs):
        super().__init__(name=name, **kwargs)
        self.x0, self.y0, self.z0 = x0, y0, z0
        self.t_sq = t_sq
        self.sheet = sheet

    def evaluate(self, point: Tuple[float, float, float]) -> float:
        x, y, z = point
        dx, dy, dz = x - self.x0, y - self.y0, z - self.z0
        return dx*dx + dz*dz - self.t_sq * dy*dy

    def _get_type(self) -> str:
        return 'KY'

    def _get_params(self) -> Tuple:
        return (self.y0, self.t_sq, self.sheet)


class ConeZ(Surface):
    """Cone along Z axis."""

    def __init__(self, x0: float, y0: float, z0: float, t_sq: float,
                 sheet: int = 0, name: Optional[str] = None, **kwargs):
        super().__init__(name=name, **kwargs)
        self.x0, self.y0, self.z0 = x0, y0, z0
        self.t_sq = t_sq
        self.sheet = sheet

    def evaluate(self, point: Tuple[float, float, float]) -> float:
        x, y, z = point
        dx, dy, dz = x - self.x0, y - self.y0, z - self.z0
        return dx*dx + dy*dy - self.t_sq * dz*dz

    def _get_type(self) -> str:
        return 'KZ'

    def _get_params(self) -> Tuple:
        return (self.z0, self.t_sq, self.sheet)


class Box(Surface):
    """Axis-aligned rectangular parallelepiped (RPP macrobody).

    A box is defined by its extent in each dimension.
    Unlike other surfaces, Box is a closed region.

    Use box.interior() to get the inside region.
    """

    def __init__(self, xmin: float, xmax: float,
                 ymin: float, ymax: float,
                 zmin: float, zmax: float,
                 name: Optional[str] = None, **kwargs):
        """
        Args:
            xmin, xmax: X extent.
            ymin, ymax: Y extent.
            zmin, zmax: Z extent.
        """
        super().__init__(name=name, **kwargs)
        if xmin >= xmax or ymin >= ymax or zmin >= zmax:
            raise ValueError("Box min values must be less than max values")
        self.xmin, self.xmax = xmin, xmax
        self.ymin, self.ymax = ymin, ymax
        self.zmin, self.zmax = zmin, zmax

    def evaluate(self, point: Tuple[float, float, float]) -> float:
        """For Box, negative = inside, positive = outside."""
        x, y, z = point
        # Distance to nearest face (negative if inside)
        dx = max(self.xmin - x, x - self.xmax, 0)
        dy = max(self.ymin - y, y - self.ymax, 0)
        dz = max(self.zmin - z, z - self.zmax, 0)

        if dx == 0 and dy == 0 and dz == 0:
            # Inside: return negative of distance to nearest face
            return -min(x - self.xmin, self.xmax - x,
                       y - self.ymin, self.ymax - y,
                       z - self.zmin, self.zmax - z)
        # Outside: return positive distance
        return math.sqrt(dx*dx + dy*dy + dz*dz)

    def _get_type(self) -> str:
        return 'RPP'

    def _get_params(self) -> Tuple:
        return (self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax)

    @property
    def center(self) -> Tuple[float, float, float]:
        """Get box center."""
        return (
            (self.xmin + self.xmax) / 2,
            (self.ymin + self.ymax) / 2,
            (self.zmin + self.zmax) / 2
        )

    @property
    def dimensions(self) -> Tuple[float, float, float]:
        """Get box dimensions (dx, dy, dz)."""
        return (
            self.xmax - self.xmin,
            self.ymax - self.ymin,
            self.zmax - self.zmin
        )


class RCC(Surface):
    """Right Circular Cylinder macrobody.

    Defined by base point, height vector, and radius.
    """

    def __init__(self, base_x: float, base_y: float, base_z: float,
                 height_x: float, height_y: float, height_z: float,
                 radius: float, name: Optional[str] = None, **kwargs):
        """
        Args:
            base_x, base_y, base_z: Base center point.
            height_x, height_y, height_z: Height vector (direction * length).
            radius: Cylinder radius.
        """
        super().__init__(name=name, **kwargs)
        if radius <= 0:
            raise ValueError("RCC radius must be positive")
        self.base_x, self.base_y, self.base_z = base_x, base_y, base_z
        self.height_x, self.height_y, self.height_z = height_x, height_y, height_z
        self.radius = radius

        # Pre-compute normalized axis
        h = math.sqrt(height_x**2 + height_y**2 + height_z**2)
        if h < 1e-10:
            raise ValueError("RCC height vector cannot be zero")
        self._axis = (height_x/h, height_y/h, height_z/h)
        self._height = h

    def evaluate(self, point: Tuple[float, float, float]) -> float:
        """Evaluate RCC: negative inside, positive outside."""
        x, y, z = point
        # Vector from base to point
        px, py, pz = x - self.base_x, y - self.base_y, z - self.base_z

        # Project onto axis
        ax, ay, az = self._axis
        proj = px*ax + py*ay + pz*az

        # Check axial bounds
        if proj < 0:
            # Below base
            return abs(proj)
        elif proj > self._height:
            # Above top
            return proj - self._height

        # Radial distance
        rx = px - proj * ax
        ry = py - proj * ay
        rz = pz - proj * az
        r = math.sqrt(rx*rx + ry*ry + rz*rz)

        return r - self.radius

    def _get_type(self) -> str:
        return 'RCC'

    def _get_params(self) -> Tuple:
        return (self.base_x, self.base_y, self.base_z,
                self.height_x, self.height_y, self.height_z,
                self.radius)


class TorusZ(Surface):
    """Torus with axis along Z: ((√(x²+y²) - R)² + z²) = r²

    A torus (donut shape) centered at a point with major radius R
    and minor radius r.
    """

    def __init__(self, x0: float, y0: float, z0: float,
                 major_radius: float, minor_radius: float,
                 name: Optional[str] = None, **kwargs):
        """
        Args:
            x0, y0, z0: Center of torus.
            major_radius: Distance from center to tube center (R).
            minor_radius: Radius of the tube (r).
        """
        super().__init__(name=name, **kwargs)
        if major_radius <= 0 or minor_radius <= 0:
            raise ValueError("Torus radii must be positive")
        self.x0, self.y0, self.z0 = x0, y0, z0
        self.major_radius = major_radius
        self.minor_radius = minor_radius

    def evaluate(self, point: Tuple[float, float, float]) -> float:
        x, y, z = point
        dx, dy, dz = x - self.x0, y - self.y0, z - self.z0
        rho = math.sqrt(dx*dx + dy*dy)
        return (rho - self.major_radius)**2 + dz*dz - self.minor_radius**2

    def _get_type(self) -> str:
        return 'TZ'

    def _get_params(self) -> Tuple:
        return (self.x0, self.y0, self.z0, self.major_radius, self.minor_radius)


class TorusX(Surface):
    """Torus with axis along X."""

    def __init__(self, x0: float, y0: float, z0: float,
                 major_radius: float, minor_radius: float,
                 name: Optional[str] = None, **kwargs):
        super().__init__(name=name, **kwargs)
        if major_radius <= 0 or minor_radius <= 0:
            raise ValueError("Torus radii must be positive")
        self.x0, self.y0, self.z0 = x0, y0, z0
        self.major_radius = major_radius
        self.minor_radius = minor_radius

    def evaluate(self, point: Tuple[float, float, float]) -> float:
        x, y, z = point
        dx, dy, dz = x - self.x0, y - self.y0, z - self.z0
        rho = math.sqrt(dy*dy + dz*dz)
        return (rho - self.major_radius)**2 + dx*dx - self.minor_radius**2

    def _get_type(self) -> str:
        return 'TX'

    def _get_params(self) -> Tuple:
        return (self.x0, self.y0, self.z0, self.major_radius, self.minor_radius)


class TorusY(Surface):
    """Torus with axis along Y."""

    def __init__(self, x0: float, y0: float, z0: float,
                 major_radius: float, minor_radius: float,
                 name: Optional[str] = None, **kwargs):
        super().__init__(name=name, **kwargs)
        if major_radius <= 0 or minor_radius <= 0:
            raise ValueError("Torus radii must be positive")
        self.x0, self.y0, self.z0 = x0, y0, z0
        self.major_radius = major_radius
        self.minor_radius = minor_radius

    def evaluate(self, point: Tuple[float, float, float]) -> float:
        x, y, z = point
        dx, dy, dz = x - self.x0, y - self.y0, z - self.z0
        rho = math.sqrt(dx*dx + dz*dz)
        return (rho - self.major_radius)**2 + dy*dy - self.minor_radius**2

    def _get_type(self) -> str:
        return 'TY'

    def _get_params(self) -> Tuple:
        return (self.x0, self.y0, self.z0, self.major_radius, self.minor_radius)


class Quadric(Surface):
    """General quadric surface (GQ).

    Ax² + By² + Cz² + Dxy + Eyz + Fzx + Gx + Hy + Jz + K = 0

    This is the most general second-degree surface and can represent
    any conic section: ellipsoid, hyperboloid, paraboloid, cone, cylinder, etc.
    """

    def __init__(self, A: float, B: float, C: float,
                 D: float, E: float, F: float,
                 G: float, H: float, J: float, K: float,
                 name: Optional[str] = None, **kwargs):
        """
        Args:
            A, B, C: Coefficients of x², y², z²
            D, E, F: Coefficients of xy, yz, zx
            G, H, J: Coefficients of x, y, z
            K: Constant term
        """
        super().__init__(name=name, **kwargs)
        self.A, self.B, self.C = A, B, C
        self.D, self.E, self.F = D, E, F
        self.G, self.H, self.J = G, H, J
        self.K = K

    def evaluate(self, point: Tuple[float, float, float]) -> float:
        x, y, z = point
        return (self.A*x*x + self.B*y*y + self.C*z*z +
                self.D*x*y + self.E*y*z + self.F*z*x +
                self.G*x + self.H*y + self.J*z + self.K)

    def _get_type(self) -> str:
        return 'GQ'

    def _get_params(self) -> Tuple:
        return (self.A, self.B, self.C, self.D, self.E, self.F,
                self.G, self.H, self.J, self.K)


class TRC(Surface):
    """Truncated Right Cone macrobody.

    A cone frustum defined by base center, height vector, and two radii.
    """

    def __init__(self, base_x: float, base_y: float, base_z: float,
                 height_x: float, height_y: float, height_z: float,
                 base_radius: float, top_radius: float,
                 name: Optional[str] = None, **kwargs):
        """
        Args:
            base_x, base_y, base_z: Base center point.
            height_x, height_y, height_z: Height vector to top center.
            base_radius: Radius at base.
            top_radius: Radius at top.
        """
        super().__init__(name=name, **kwargs)
        if base_radius <= 0 or top_radius < 0:
            raise ValueError("TRC radii must be non-negative (base > 0)")
        self.base_x, self.base_y, self.base_z = base_x, base_y, base_z
        self.height_x, self.height_y, self.height_z = height_x, height_y, height_z
        self.base_radius = base_radius
        self.top_radius = top_radius

        h = math.sqrt(height_x**2 + height_y**2 + height_z**2)
        if h < 1e-10:
            raise ValueError("TRC height vector cannot be zero")
        self._height = h
        self._axis = (height_x/h, height_y/h, height_z/h)

    def evaluate(self, point: Tuple[float, float, float]) -> float:
        x, y, z = point
        px, py, pz = x - self.base_x, y - self.base_y, z - self.base_z
        ax, ay, az = self._axis
        proj = px*ax + py*ay + pz*az

        if proj < 0:
            return abs(proj)
        elif proj > self._height:
            return proj - self._height

        # Interpolate radius
        t = proj / self._height
        r_at_proj = self.base_radius + t * (self.top_radius - self.base_radius)

        rx = px - proj * ax
        ry = py - proj * ay
        rz = pz - proj * az
        r = math.sqrt(rx*rx + ry*ry + rz*rz)

        return r - r_at_proj

    def _get_type(self) -> str:
        return 'TRC'

    def _get_params(self) -> Tuple:
        return (self.base_x, self.base_y, self.base_z,
                self.height_x, self.height_y, self.height_z,
                self.base_radius, self.top_radius)


class ELL(Surface):
    """Ellipsoid macrobody.

    Defined by two focus points and the major axis length.
    """

    def __init__(self, v1_x: float, v1_y: float, v1_z: float,
                 v2_x: float, v2_y: float, v2_z: float,
                 major_axis_len: float,
                 name: Optional[str] = None, **kwargs):
        """
        Args:
            v1_x, v1_y, v1_z: First focus point.
            v2_x, v2_y, v2_z: Second focus point.
            major_axis_len: Length of major axis (2a).
        """
        super().__init__(name=name, **kwargs)
        self.v1_x, self.v1_y, self.v1_z = v1_x, v1_y, v1_z
        self.v2_x, self.v2_y, self.v2_z = v2_x, v2_y, v2_z
        self.major_axis_len = major_axis_len

        # Distance between foci
        dx = v2_x - v1_x
        dy = v2_y - v1_y
        dz = v2_z - v1_z
        self._focal_dist = math.sqrt(dx*dx + dy*dy + dz*dz)

        if major_axis_len <= self._focal_dist:
            raise ValueError("Major axis must be longer than distance between foci")

    def evaluate(self, point: Tuple[float, float, float]) -> float:
        x, y, z = point
        # Distance to first focus
        d1 = math.sqrt((x-self.v1_x)**2 + (y-self.v1_y)**2 + (z-self.v1_z)**2)
        # Distance to second focus
        d2 = math.sqrt((x-self.v2_x)**2 + (y-self.v2_y)**2 + (z-self.v2_z)**2)
        # Ellipsoid: d1 + d2 = major_axis_len
        return d1 + d2 - self.major_axis_len

    def _get_type(self) -> str:
        return 'ELL'

    def _get_params(self) -> Tuple:
        return (self.v1_x, self.v1_y, self.v1_z,
                self.v2_x, self.v2_y, self.v2_z,
                self.major_axis_len)


class REC(Surface):
    """Right Elliptical Cylinder macrobody.

    A cylinder with elliptical cross-section.
    """

    def __init__(self, base_x: float, base_y: float, base_z: float,
                 height_x: float, height_y: float, height_z: float,
                 axis1_x: float, axis1_y: float, axis1_z: float,
                 axis2_x: float, axis2_y: float, axis2_z: float,
                 name: Optional[str] = None, **kwargs):
        """
        Args:
            base_x, base_y, base_z: Base center point.
            height_x, height_y, height_z: Height vector.
            axis1_x, axis1_y, axis1_z: First ellipse semi-axis vector.
            axis2_x, axis2_y, axis2_z: Second ellipse semi-axis vector.
        """
        super().__init__(name=name, **kwargs)
        self.base_x, self.base_y, self.base_z = base_x, base_y, base_z
        self.height_x, self.height_y, self.height_z = height_x, height_y, height_z
        self.axis1_x, self.axis1_y, self.axis1_z = axis1_x, axis1_y, axis1_z
        self.axis2_x, self.axis2_y, self.axis2_z = axis2_x, axis2_y, axis2_z

        self._height = math.sqrt(height_x**2 + height_y**2 + height_z**2)
        self._a = math.sqrt(axis1_x**2 + axis1_y**2 + axis1_z**2)
        self._b = math.sqrt(axis2_x**2 + axis2_y**2 + axis2_z**2)

    def evaluate(self, point: Tuple[float, float, float]) -> float:
        # Simplified evaluation - full implementation would project onto ellipse plane
        x, y, z = point
        px, py, pz = x - self.base_x, y - self.base_y, z - self.base_z

        # Normalize height axis
        h = self._height
        if h < 1e-10:
            return 1.0  # Degenerate
        hx, hy, hz = self.height_x/h, self.height_y/h, self.height_z/h

        # Project onto axis
        proj = px*hx + py*hy + pz*hz
        if proj < 0 or proj > h:
            return abs(proj - h/2) - h/2

        # This is a simplified radial check
        rx = px - proj * hx
        ry = py - proj * hy
        rz = pz - proj * hz
        r = math.sqrt(rx*rx + ry*ry + rz*rz)
        return r - max(self._a, self._b)

    def _get_type(self) -> str:
        return 'REC'

    def _get_params(self) -> Tuple:
        return (self.base_x, self.base_y, self.base_z,
                self.height_x, self.height_y, self.height_z,
                self.axis1_x, self.axis1_y, self.axis1_z,
                self.axis2_x, self.axis2_y, self.axis2_z)


class WED(Surface):
    """Wedge macrobody.

    A triangular prism (wedge) defined by a vertex and three edge vectors.
    """

    def __init__(self, vertex_x: float, vertex_y: float, vertex_z: float,
                 v1_x: float, v1_y: float, v1_z: float,
                 v2_x: float, v2_y: float, v2_z: float,
                 v3_x: float, v3_y: float, v3_z: float,
                 name: Optional[str] = None, **kwargs):
        """
        Args:
            vertex_x, vertex_y, vertex_z: Base vertex of wedge.
            v1_x, v1_y, v1_z: First edge vector (forms base triangle).
            v2_x, v2_y, v2_z: Second edge vector (forms base triangle).
            v3_x, v3_y, v3_z: Height vector (prism extrusion direction).
        """
        super().__init__(name=name, **kwargs)
        self.vertex_x, self.vertex_y, self.vertex_z = vertex_x, vertex_y, vertex_z
        self.v1_x, self.v1_y, self.v1_z = v1_x, v1_y, v1_z
        self.v2_x, self.v2_y, self.v2_z = v2_x, v2_y, v2_z
        self.v3_x, self.v3_y, self.v3_z = v3_x, v3_y, v3_z

    def evaluate(self, point: Tuple[float, float, float]) -> float:
        # Simplified - returns approximate signed distance
        return 0.0  # Proper evaluation requires half-plane checks

    def _get_type(self) -> str:
        return 'WED'

    def _get_params(self) -> Tuple:
        return (self.vertex_x, self.vertex_y, self.vertex_z,
                self.v1_x, self.v1_y, self.v1_z,
                self.v2_x, self.v2_y, self.v2_z,
                self.v3_x, self.v3_y, self.v3_z)


class RHP(Surface):
    """Right Hexagonal Prism macrobody.

    A hexagonal prism defined by base center, height vector, and
    three vectors from center to vertices (every other vertex).
    """

    def __init__(self, base_x: float, base_y: float, base_z: float,
                 height_x: float, height_y: float, height_z: float,
                 r1_x: float, r1_y: float, r1_z: float,
                 r2_x: float, r2_y: float, r2_z: float,
                 r3_x: float, r3_y: float, r3_z: float,
                 name: Optional[str] = None, **kwargs):
        """
        Args:
            base_x, base_y, base_z: Base center point.
            height_x, height_y, height_z: Height vector.
            r1, r2, r3: Vectors from center to alternating vertices.
        """
        super().__init__(name=name, **kwargs)
        self.base_x, self.base_y, self.base_z = base_x, base_y, base_z
        self.height_x, self.height_y, self.height_z = height_x, height_y, height_z
        self.r1_x, self.r1_y, self.r1_z = r1_x, r1_y, r1_z
        self.r2_x, self.r2_y, self.r2_z = r2_x, r2_y, r2_z
        self.r3_x, self.r3_y, self.r3_z = r3_x, r3_y, r3_z

    def evaluate(self, point: Tuple[float, float, float]) -> float:
        # Simplified - returns approximate signed distance
        return 0.0  # Proper evaluation requires half-plane checks

    def _get_type(self) -> str:
        return 'RHP'

    def _get_params(self) -> Tuple:
        return (self.base_x, self.base_y, self.base_z,
                self.height_x, self.height_y, self.height_z,
                self.r1_x, self.r1_y, self.r1_z,
                self.r2_x, self.r2_y, self.r2_z,
                self.r3_x, self.r3_y, self.r3_z)


class GeneralBox(Surface):
    """General box macrobody (arbitrarily oriented).

    A parallelepiped defined by a corner and three edge vectors.
    Unlike Box (RPP), this can be arbitrarily oriented.
    """

    def __init__(self, corner_x: float, corner_y: float, corner_z: float,
                 v1_x: float, v1_y: float, v1_z: float,
                 v2_x: float, v2_y: float, v2_z: float,
                 v3_x: float, v3_y: float, v3_z: float,
                 name: Optional[str] = None, **kwargs):
        """
        Args:
            corner_x, corner_y, corner_z: Corner point.
            v1, v2, v3: Three edge vectors from corner.
        """
        super().__init__(name=name, **kwargs)
        self.corner_x, self.corner_y, self.corner_z = corner_x, corner_y, corner_z
        self.v1_x, self.v1_y, self.v1_z = v1_x, v1_y, v1_z
        self.v2_x, self.v2_y, self.v2_z = v2_x, v2_y, v2_z
        self.v3_x, self.v3_y, self.v3_z = v3_x, v3_y, v3_z

    def evaluate(self, point: Tuple[float, float, float]) -> float:
        # Simplified - returns approximate signed distance
        return 0.0  # Proper evaluation requires half-plane checks

    def _get_type(self) -> str:
        return 'BOX'

    def _get_params(self) -> Tuple:
        return (self.corner_x, self.corner_y, self.corner_z,
                self.v1_x, self.v1_y, self.v1_z,
                self.v2_x, self.v2_y, self.v2_z,
                self.v3_x, self.v3_y, self.v3_z)


# Import Halfspace for type annotations
from .geometry import Halfspace
