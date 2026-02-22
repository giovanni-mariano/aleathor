# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""
Collection classes for cells, surfaces, and universes.

These provide Pythonic iteration, filtering, and selection over model contents.
Collections are views into the model - modifying the model updates the collection.
"""

from __future__ import annotations
from typing import (
    TYPE_CHECKING, Iterator, List, Optional, Tuple, Union,
    Callable, Any, Dict, Set
)
from dataclasses import dataclass, field
import math

if TYPE_CHECKING:
    from .model import Model, Universe


@dataclass
class Cell:
    """Mutable view of a cell in the model.

    Provides access to cell properties and supports mutation via property
    setters.  Mutations update the Python ``Cell`` object, mark the model
    dirty, and invalidate the cached C info so the next read reflects the
    change.

    Attributes:
        id: Cell ID (MCNP cell number)
        material: Material ID (0 for void)
        density: Material density
        universe: Universe this cell belongs to
        fill: Universe ID filling this cell (None if not filled)
        region: Python Region object (None for loaded-only models)
        importance: Particle importance (1.0 if not set)
        depth: Hierarchy depth (0 = root universe)
    """
    _model: 'Model'
    _index: int
    _depth: int = field(default=0, repr=False)
    _info: Optional[Dict[str, Any]] = field(default=None, repr=False, compare=False)
    _cell_id: Optional[int] = field(default=None, repr=False, compare=False)

    @property
    def id(self) -> int:
        """Cell ID (MCNP cell number)."""
        return self._get_info()['cell_id']

    @property
    def material(self) -> int:
        """Material ID (0 for void)."""
        return self._get_info()['material_id']

    @material.setter
    def material(self, value: int) -> None:
        self._require_py_cell().material = value
        self._model._dirty = True
        self._invalidate_cache()

    @property
    def density(self) -> float:
        """Material density (always positive)."""
        return self._get_info()['density']

    @density.setter
    def density(self, value: float) -> None:
        self._require_py_cell().density = abs(value)
        self._model._dirty = True
        self._invalidate_cache()

    @property
    def density_unit(self) -> str:
        """Density unit — ``"g/cm3"`` or ``"atoms/b-cm"``."""
        py_cell = self._get_py_cell()
        if py_cell is not None and hasattr(py_cell, 'density_unit'):
            return py_cell.density_unit
        # Fall back to C info for loaded-only models
        info = self._get_info()
        if info.get('is_mass_density', True):
            return "g/cm3"
        return "atoms/b-cm"

    @density_unit.setter
    def density_unit(self, value: str) -> None:
        if value not in ("g/cm3", "atoms/b-cm"):
            raise ValueError(
                f"density_unit must be 'g/cm3' or 'atoms/b-cm', got {value!r}"
            )
        self._require_py_cell().density_unit = value
        self._model._dirty = True
        self._invalidate_cache()

    @property
    def universe(self) -> int:
        """Universe this cell belongs to."""
        return self._get_info()['universe_id']

    @property
    def fill(self) -> Optional[int]:
        """Universe ID filling this cell, or None."""
        fill = self._get_info().get('fill_universe', -1)
        return fill if fill >= 0 else None

    @fill.setter
    def fill(self, value: Optional[int]) -> None:
        py_cell = self._require_py_cell()
        py_cell.fill = value
        py_cell.fill_transform = 0 if value is not None else None
        if not self._model._dirty and self._model._sys is not None:
            fill_id = value if value is not None else -1
            self._model._sys.set_fill(self._index, fill_id, 0)
        else:
            self._model._dirty = True
        self._invalidate_cache()

    @property
    def name(self) -> Optional[str]:
        """Cell name (if assigned during construction)."""
        py_cell = self._get_py_cell()
        if py_cell is not None and hasattr(py_cell, 'name'):
            return py_cell.name
        return None

    @name.setter
    def name(self, value: Optional[str]) -> None:
        self._require_py_cell().name = value

    @property
    def bounds(self) -> Tuple[float, float, float, float, float, float]:
        """Bounding box (xmin, xmax, ymin, ymax, zmin, zmax)."""
        return self._get_info().get('bbox', (0, 0, 0, 0, 0, 0))

    @property
    def is_void(self) -> bool:
        """True if cell has no material."""
        return self.material == 0

    @property
    def is_filled(self) -> bool:
        """True if cell is filled with another universe."""
        return self.fill is not None

    @property
    def region(self) -> Any:
        """Python Region object for the cell.

        For models constructed via the Python API, returns the Region
        used to define the cell. For models loaded from file (e.g. MCNP),
        returns an ``_ImportedRegion`` whose ``repr()`` shows the
        reconstructed CSG expression (e.g. ``-1 & +2``) and whose
        ``get_surfaces()`` returns the referenced surface IDs.
        """
        py_cell = self._get_py_cell()
        if py_cell is not None and hasattr(py_cell, 'region'):
            return py_cell.region
        return None

    @property
    def importance(self) -> float:
        """Particle importance (defaults to 1.0 if not set)."""
        py_cell = self._get_py_cell()
        if py_cell is not None and hasattr(py_cell, 'importance'):
            return py_cell.importance
        return 1.0

    @importance.setter
    def importance(self, value: float) -> None:
        self._require_py_cell().importance = value

    @property
    def depth(self) -> int:
        """Hierarchy depth (0 = root universe).

        Set by cells_at() for hierarchical queries. For cells obtained
        via iteration or get_cell(), this is always 0.
        """
        return self._depth

    @property
    def is_lattice(self) -> bool:
        """True if this cell defines a lattice."""
        return self._get_info().get('lat_type', 0) != 0

    @property
    def lattice_type(self) -> int:
        """Lattice type: 0=none, 1=rectangular, 2=hexagonal."""
        return self._get_info().get('lat_type', 0)

    @property
    def lattice_pitch(self) -> Optional[Tuple[float, float, float]]:
        """Lattice element pitch (x, y, z), or None if not a lattice."""
        return self._get_info().get('lat_pitch')

    @property
    def lattice_fill(self) -> Optional[List[int]]:
        """Lattice fill universe IDs, or None if not a lattice."""
        return self._get_info().get('lat_fill')

    @property
    def lattice_lower_left(self) -> Optional[Tuple[float, float, float]]:
        """Lattice lower-left corner, or None if not a lattice."""
        return self._get_info().get('lat_lower_left')

    @property
    def lattice_dims(self) -> Optional[Tuple[int, int, int, int, int, int]]:
        """Lattice dimensions (imin,imax,jmin,jmax,kmin,kmax), or None if not a lattice."""
        return self._get_info().get('lat_fill_dims')

    def _get_py_cell(self) -> Any:
        """Get the Python Cell object if available."""
        if hasattr(self._model, '_cells'):
            return self._model._cells.get(self.id)
        return None

    def _require_py_cell(self):
        """Return the Python Cell backing this view, or raise."""
        py_cell = self._get_py_cell()
        if py_cell is None:
            raise RuntimeError(
                f"Cell {self.id} cannot be modified (no Python backing object)"
            )
        return py_cell

    def _get_info(self) -> Dict[str, Any]:
        """Get cell info from C system (cached)."""
        if self._info is None:
            self._model._rebuild_if_needed()
            # Re-resolve index after a possible rebuild
            if self._cell_id is not None:
                resolved = self._model._sys.cell_find(self._cell_id)
                if resolved is not None:
                    self._index = resolved
            self._info = self._model._sys.get_cell_by_index(self._index)
            if self._cell_id is None:
                self._cell_id = self._info['cell_id']
        return self._info

    def _invalidate_cache(self) -> None:
        """Clear cached info (called after cell mutation)."""
        self._info = None

    def fill_with(self, universe, transform: int = 0) -> None:
        """Fill this cell with a universe.

        Args:
            universe: Universe object or universe ID (int).
            transform: Transform ID to apply (default 0 = identity).
        """
        from .model import Universe
        if isinstance(universe, Universe):
            fill_id = universe.id
        elif isinstance(universe, int):
            fill_id = universe
        else:
            raise TypeError(
                f"Expected Universe or int, got {type(universe).__name__}"
            )

        py_cell = self._require_py_cell()
        py_cell.fill = fill_id
        py_cell.fill_transform = transform
        if not self._model._dirty and self._model._sys is not None:
            self._model._sys.set_fill(self._index, fill_id, transform)
        else:
            self._model._dirty = True
        self._invalidate_cache()

    def unfill(self) -> None:
        """Remove the fill from this cell."""
        self.fill = None

    def contains(self, x: float, y: float, z: float) -> bool:
        """Check if point is inside this cell.

        Args:
            x, y, z: Point coordinates

        Returns:
            True if point is inside this cell
        """
        self._model._rebuild_if_needed()
        result = self._model._sys.find_cell(x, y, z)
        if result is None:
            return False
        cell_idx, _ = result
        return cell_idx == self._index

    def __repr__(self) -> str:
        return f"Cell({self.id}, material={self.material}, universe={self.universe})"

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Cell):
            return self._model is other._model and self._index == other._index
        return False

    def __hash__(self) -> int:
        return hash((id(self._model), self._index))


class CellCollection:
    """Collection of cells with filtering capabilities.

    Supports iteration, indexing, and filtering operations.
    Filter methods return new CellCollection instances that can be chained.

    Example:
        # Iterate all cells
        for cell in model.cells:
            print(cell.id, cell.material)

        # Filter by material
        fuel_cells = model.cells.by_material(92)
        for cell in fuel_cells:
            print(cell.id)

        # Chain filters
        u0_fuel = model.cells.by_material(92).by_universe(0)

        # Count
        print(len(fuel_cells))
    """

    def __init__(self, model: 'Model', indices: Optional[List[int]] = None):
        """
        Args:
            model: Parent model
            indices: Cell indices to include (None = all cells)
        """
        self._model = model
        self._indices = indices  # None means all cells

    def _get_indices(self) -> List[int]:
        """Get list of cell indices in this collection."""
        if self._indices is not None:
            return self._indices
        # All cells
        self._model._rebuild_if_needed()
        return list(range(self._model._sys.cell_count))

    def __len__(self) -> int:
        return len(self._get_indices())

    def __iter__(self) -> Iterator[Cell]:
        for idx in self._get_indices():
            yield Cell(self._model, idx)

    def __getitem__(self, key: Union[int, slice]) -> Union[Cell, 'CellCollection']:
        indices = self._get_indices()
        if isinstance(key, int):
            if key < 0:
                key = len(indices) + key
            if key < 0 or key >= len(indices):
                raise IndexError(f"Cell index {key} out of range")
            return Cell(self._model, indices[key])
        elif isinstance(key, slice):
            return CellCollection(self._model, indices[key])
        else:
            raise TypeError(f"Invalid index type: {type(key)}")

    def __contains__(self, cell: Cell) -> bool:
        if not isinstance(cell, Cell):
            return False
        return cell._index in self._get_indices()

    def __repr__(self) -> str:
        return f"CellCollection({len(self)} cells)"

    # =========================================================================
    # Filtering Methods
    # =========================================================================

    def by_material(self, material_id: int) -> 'CellCollection':
        """Filter cells by material ID.

        Args:
            material_id: Material ID to filter by

        Returns:
            New CellCollection with matching cells
        """
        self._model._rebuild_if_needed()

        # Use C API filter if available and we're filtering all cells
        if self._indices is None:
            filtered = self._model._sys.get_cells_by_material(material_id)
            return CellCollection(self._model, filtered)

        # Filter existing indices
        filtered = []
        for idx in self._get_indices():
            info = self._model._sys.get_cell_by_index(idx)
            if info['material_id'] == material_id:
                filtered.append(idx)
        return CellCollection(self._model, filtered)

    def by_universe(self, universe_id: int) -> 'CellCollection':
        """Filter cells by universe ID.

        Args:
            universe_id: Universe ID to filter by

        Returns:
            New CellCollection with matching cells
        """
        self._model._rebuild_if_needed()

        # Use C API filter if available and we're filtering all cells
        if self._indices is None:
            filtered = self._model._sys.get_cells_by_universe(universe_id)
            return CellCollection(self._model, filtered)

        # Filter existing indices
        filtered = []
        for idx in self._get_indices():
            info = self._model._sys.get_cell_by_index(idx)
            if info['universe_id'] == universe_id:
                filtered.append(idx)
        return CellCollection(self._model, filtered)

    def by_fill(self, fill_universe: Optional[int] = None) -> 'CellCollection':
        """Filter cells by fill universe.

        Args:
            fill_universe: Universe ID that fills the cell.
                          If None, returns all cells with any fill.

        Returns:
            New CellCollection with matching cells
        """
        self._model._rebuild_if_needed()

        if fill_universe is not None and self._indices is None:
            # Use C API filter
            filtered = self._model._sys.get_cells_filling_universe(fill_universe)
            return CellCollection(self._model, filtered)

        # Filter existing indices
        filtered = []
        for idx in self._get_indices():
            info = self._model._sys.get_cell_by_index(idx)
            fill = info.get('fill_universe', -1)
            if fill_universe is None:
                if fill >= 0:
                    filtered.append(idx)
            else:
                if fill == fill_universe:
                    filtered.append(idx)
        return CellCollection(self._model, filtered)

    def filter(self, predicate: Callable[[Cell], bool]) -> 'CellCollection':
        """Filter cells using a custom predicate function.

        Args:
            predicate: Function that takes a Cell and returns True to include

        Returns:
            New CellCollection with matching cells

        Example:
            # Get cells with density > 5
            heavy = model.cells.filter(lambda c: abs(c.density) > 5)
        """
        filtered = []
        for idx in self._get_indices():
            cell = Cell(self._model, idx)
            if predicate(cell):
                filtered.append(idx)
        return CellCollection(self._model, filtered)

    # =========================================================================
    # Utility Methods
    # =========================================================================

    def get(self, cell_id: int) -> Optional[Cell]:
        """Get a cell by its ID.

        Uses O(1) lookup via the C library when searching all cells.

        Args:
            cell_id: Cell ID to find

        Returns:
            Cell if found, None otherwise
        """
        self._model._rebuild_if_needed()

        if self._indices is None:
            # O(1) lookup via C API
            idx = self._model._sys.cell_find(cell_id)
            if idx is None:
                return None
            return Cell(self._model, idx)

        # Filtered collection: linear scan
        for idx in self._indices:
            info = self._model._sys.get_cell_by_index(idx)
            if info['cell_id'] == cell_id:
                return Cell(self._model, idx)
        return None

    def ids(self) -> List[int]:
        """Get list of cell IDs in this collection."""
        return [Cell(self._model, idx).id for idx in self._get_indices()]

    def materials(self) -> Set[int]:
        """Get set of unique material IDs in this collection."""
        return {Cell(self._model, idx).material for idx in self._get_indices()}

    def universes(self) -> Set[int]:
        """Get set of unique universe IDs in this collection."""
        return {Cell(self._model, idx).universe for idx in self._get_indices()}

    def to_list(self) -> List[Cell]:
        """Convert to a list of Cell objects."""
        return list(self)


@dataclass
class TraceSegment:
    """A segment of a ray trace through geometry.

    Attributes:
        cell: The cell this segment passes through (None for void)
        t_enter: Entry distance along ray
        t_exit: Exit distance along ray
        material: Material ID (0 for void)
        density: Material density
    """
    cell: Optional[Cell]
    t_enter: float
    t_exit: float
    material: int
    density: float

    @property
    def length(self) -> float:
        """Length of this segment."""
        if self.t_exit > 1e30:
            return float('inf')
        return self.t_exit - self.t_enter

    @property
    def is_void(self) -> bool:
        """True if segment is in void."""
        return self.cell is None

    def __repr__(self) -> str:
        cell_str = f"cell={self.cell.id}" if self.cell else "void"
        return f"TraceSegment({cell_str}, length={self.length:.4f})"


class TraceResult:
    """Result of a ray trace through geometry.

    Supports iteration over segments and provides convenience methods
    for analyzing the trace.

    Example:
        trace = model.trace(origin=(0,0,0), direction=(1,0,0))

        for seg in trace:
            print(f"{seg.cell.id}: {seg.length} cm")

        # Total path length through material
        fuel_length = trace.path_length(material=92)
    """

    def __init__(self, model: 'Model', segments: List[Dict[str, Any]]):
        self._model = model
        self._raw_segments = segments
        self._segments: Optional[List[TraceSegment]] = None

    def _build_segments(self) -> List[TraceSegment]:
        """Convert raw segments to TraceSegment objects.

        Always returns Cell for consistency.
        """
        if self._segments is not None:
            return self._segments

        self._segments = []
        for seg in self._raw_segments:
            cell_id = seg['cell_id']
            cell = None
            if cell_id >= 0:
                # Always use Cell for consistent return type
                idx = self._model._sys.cell_find(cell_id)
                if idx is not None:
                    cell = Cell(self._model, idx)

            self._segments.append(TraceSegment(
                cell=cell,
                t_enter=seg['t_enter'],
                t_exit=seg['t_exit'],
                material=seg['material_id'],
                density=seg['density']
            ))

        return self._segments

    def __len__(self) -> int:
        return len(self._raw_segments)

    def __iter__(self) -> Iterator[TraceSegment]:
        for seg in self._build_segments():
            yield seg

    def __getitem__(self, key: int) -> TraceSegment:
        return self._build_segments()[key]

    def path_length(self, material: Optional[int] = None) -> float:
        """Calculate total path length through material.

        Args:
            material: Material ID to sum (None = all materials)

        Returns:
            Total path length in cm
        """
        total = 0.0
        for seg in self._build_segments():
            if seg.t_exit > 1e30:
                continue
            if material is None or seg.material == material:
                total += seg.length
        return total

    def cells_hit(self) -> List[Optional[Cell]]:
        """Get list of cells hit by the ray (in order)."""
        return [seg.cell for seg in self._build_segments()]

    def materials_hit(self) -> Set[int]:
        """Get set of materials encountered."""
        return {seg.material for seg in self._build_segments()}

    def __repr__(self) -> str:
        return f"TraceResult({len(self)} segments)"
