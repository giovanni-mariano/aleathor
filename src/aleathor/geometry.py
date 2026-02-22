# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""
CSG region definitions supporting boolean operations via Python operators.

Regions represent volumes in space defined by CSG operations on surfaces.
Use Python operators to combine regions:
    - `a & b` or `a * b`: intersection (AND)
    - `a | b`: union (OR)
    - `~a` or `-a`: complement (NOT)
    - `a - b`: difference (a AND NOT b)
"""

from __future__ import annotations
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, List, Optional, Tuple, Set
import copy

if TYPE_CHECKING:
    from .surfaces import Surface


class Region(ABC):
    """Abstract base class for CSG regions.

    A region represents a volume in 3D space defined by CSG operations.
    Regions can be combined using Python operators:

        intersection = region_a & region_b
        union = region_a | region_b
        complement = ~region_a
        difference = region_a - region_b

    Regions are immutable - operations create new Region objects.
    """

    @abstractmethod
    def __contains__(self, point: Tuple[float, float, float]) -> bool:
        """Test if point is inside region."""
        pass

    @abstractmethod
    def get_surfaces(self) -> Set['Surface']:
        """Get all surfaces used in this region."""
        pass

    @abstractmethod
    def _to_csg(self, model: 'Model') -> int:
        """Convert to CSG node in the model's internal system."""
        pass

    def __and__(self, other: 'Region') -> 'Intersection':
        """Intersection (AND): region_a & region_b"""
        return Intersection(self, other)

    def __mul__(self, other: 'Region') -> 'Intersection':
        """Intersection (AND): region_a * region_b (alternative syntax)"""
        return Intersection(self, other)

    def __or__(self, other: 'Region') -> 'Union':
        """Union (OR): region_a | region_b"""
        return Union(self, other)

    def __invert__(self) -> 'Complement':
        """Complement (NOT): ~region"""
        return Complement(self)

    def __neg__(self) -> 'Complement':
        """Complement (NOT): -region (alternative syntax)"""
        return Complement(self)

    def __sub__(self, other: 'Region') -> 'Intersection':
        """Difference: region_a - region_b = region_a & ~region_b"""
        return Intersection(self, Complement(other))

    def __pos__(self) -> 'Region':
        """Positive: +region returns self"""
        return self

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(...)"


class Halfspace(Region):
    """Region defined by one side of a surface.

    A halfspace is the fundamental building block - the region on one side
    of a surface. For a plane ax+by+cz+d=0:
        - positive halfspace: ax+by+cz+d > 0
        - negative halfspace: ax+by+cz+d < 0

    For closed surfaces like spheres:
        - negative halfspace: inside the sphere
        - positive halfspace: outside the sphere

    Created by calling surface.positive() or surface.negative(), or using
    +surface and -surface syntax.
    """

    def __init__(self, surface: 'Surface', positive: bool = True):
        """
        Args:
            surface: The surface defining the halfspace.
            positive: True for positive side, False for negative side.
        """
        self.surface = surface
        self.positive = positive

    def __contains__(self, point: Tuple[float, float, float]) -> bool:
        """Test if point is in this halfspace."""
        value = self.surface.evaluate(point)
        if self.positive:
            return value > 0
        else:
            return value < 0

    def get_surfaces(self) -> Set['Surface']:
        """Get the single surface defining this halfspace."""
        return {self.surface}

    def _to_csg(self, model: 'Model') -> int:
        """Convert to CSG node."""
        return model._get_or_create_halfspace_node(self.surface, self.positive)

    def __repr__(self) -> str:
        sign = '+' if self.positive else '-'
        return f"{sign}{self.surface}"

    def __invert__(self) -> 'Halfspace':
        """Complement of halfspace is the opposite halfspace."""
        return Halfspace(self.surface, not self.positive)

    def __neg__(self) -> 'Halfspace':
        """Complement using -halfspace syntax."""
        return Halfspace(self.surface, not self.positive)


class Intersection(Region):
    """Region formed by intersection (AND) of multiple regions.

    Points are inside an intersection only if they are inside ALL
    constituent regions.

    Can contain other Intersections, which are flattened.
    """

    def __init__(self, *regions: Region):
        """
        Args:
            *regions: Two or more regions to intersect.
        """
        # Flatten nested intersections
        self.regions: List[Region] = []
        for r in regions:
            if isinstance(r, Intersection):
                self.regions.extend(r.regions)
            else:
                self.regions.append(r)

        if len(self.regions) < 2:
            raise ValueError("Intersection requires at least 2 regions")

    def __contains__(self, point: Tuple[float, float, float]) -> bool:
        """Point is inside only if inside ALL regions."""
        return all(point in r for r in self.regions)

    def get_surfaces(self) -> Set['Surface']:
        """Get all surfaces from all constituent regions."""
        surfaces = set()
        for r in self.regions:
            surfaces.update(r.get_surfaces())
        return surfaces

    def _to_csg(self, model: 'Model') -> int:
        """Convert to CSG intersection node."""
        nodes = [r._to_csg(model) for r in self.regions]
        return model._sys.create_intersection_many(nodes)

    def __repr__(self) -> str:
        return " & ".join(repr(r) for r in self.regions)

    def __and__(self, other: Region) -> 'Intersection':
        """Extend intersection with another region."""
        if isinstance(other, Intersection):
            return Intersection(*self.regions, *other.regions)
        return Intersection(*self.regions, other)


class Union(Region):
    """Region formed by union (OR) of multiple regions.

    Points are inside a union if they are inside ANY of the
    constituent regions.

    Can contain other Unions, which are flattened.
    """

    def __init__(self, *regions: Region):
        """
        Args:
            *regions: Two or more regions to unite.
        """
        # Flatten nested unions
        self.regions: List[Region] = []
        for r in regions:
            if isinstance(r, Union):
                self.regions.extend(r.regions)
            else:
                self.regions.append(r)

        if len(self.regions) < 2:
            raise ValueError("Union requires at least 2 regions")

    def __contains__(self, point: Tuple[float, float, float]) -> bool:
        """Point is inside if inside ANY region."""
        return any(point in r for r in self.regions)

    def get_surfaces(self) -> Set['Surface']:
        """Get all surfaces from all constituent regions."""
        surfaces = set()
        for r in self.regions:
            surfaces.update(r.get_surfaces())
        return surfaces

    def _to_csg(self, model: 'Model') -> int:
        """Convert to CSG union node."""
        nodes = [r._to_csg(model) for r in self.regions]
        return model._sys.create_union_many(nodes)

    def __repr__(self) -> str:
        return " | ".join(repr(r) for r in self.regions)

    def __or__(self, other: Region) -> 'Union':
        """Extend union with another region."""
        if isinstance(other, Union):
            return Union(*self.regions, *other.regions)
        return Union(*self.regions, other)


class Complement(Region):
    """Region formed by complement (NOT) of another region.

    Points are inside a complement if they are OUTSIDE the original region.
    """

    def __init__(self, region: Region):
        """
        Args:
            region: Region to complement.
        """
        # Double negation elimination
        if isinstance(region, Complement):
            # ~(~A) = A, but we still need to wrap it
            self._inner = region.region
            self._double_neg = True
        else:
            self._inner = region
            self._double_neg = False

        self.region = region

    def __contains__(self, point: Tuple[float, float, float]) -> bool:
        """Point is inside complement if OUTSIDE original region."""
        return point not in self.region

    def get_surfaces(self) -> Set['Surface']:
        """Get surfaces from the complemented region."""
        return self.region.get_surfaces()

    def _to_csg(self, model: 'Model') -> int:
        """Convert to CSG complement node."""
        inner_node = self.region._to_csg(model)
        return model._sys.create_complement(inner_node)

    def __repr__(self) -> str:
        return f"~({self.region})"

    def __invert__(self) -> Region:
        """Double negation: ~(~A) = A"""
        return self.region


# Avoid circular import
from .model import Model
