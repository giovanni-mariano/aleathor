# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""Tests for CSG tree inspection and region reconstruction."""

import pytest
from pathlib import Path

DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture
def sphere_model():
    """Load sphere_simple.inp: cell 1 = -1, cell 2 = +1."""
    import aleathor as ath
    return ath.read_mcnp(DATA_DIR / "sphere_simple.inp")


@pytest.fixture
def pin_cell_model():
    """Load pin_cell.inp with multiple cells."""
    import aleathor as ath
    return ath.read_mcnp(DATA_DIR / "pin_cell.inp")


@pytest.fixture
def boolean_model():
    """Load boolean_complex.inp with unions and complements."""
    import aleathor as ath
    return ath.read_mcnp(DATA_DIR / "boolean_complex.inp")


# ---------------------------------------------------------------------------
# _SurfaceRef
# ---------------------------------------------------------------------------

class TestSurfaceRef:
    def test_repr(self):
        from aleathor.io import _SurfaceRef
        s = _SurfaceRef(42)
        assert repr(s) == "42"

    def test_equality(self):
        from aleathor.io import _SurfaceRef
        a = _SurfaceRef(5)
        b = _SurfaceRef(5)
        c = _SurfaceRef(6)
        assert a == b
        assert a != c

    def test_hashable(self):
        from aleathor.io import _SurfaceRef
        a = _SurfaceRef(5)
        b = _SurfaceRef(5)
        assert hash(a) == hash(b)
        assert len({a, b}) == 1


# ---------------------------------------------------------------------------
# node_tree() C binding
# ---------------------------------------------------------------------------

class TestNodeTree:
    def test_primitive_negative(self, sphere_model):
        """Cell 1 = -1 (negative halfspace of surface 1)."""
        cell = sphere_model._cells[1]
        raw = sphere_model._sys.node_tree(cell.region._node_id)
        # Should be a primitive: (0, surface_id, sense)
        assert raw[0] == 0  # CSG_OP_PRIMITIVE
        assert raw[1] == 1  # surface ID
        assert raw[2] == -1  # negative sense

    def test_primitive_positive(self, sphere_model):
        """Cell 2 = +1 (positive halfspace of surface 1)."""
        cell = sphere_model._cells[2]
        raw = sphere_model._sys.node_tree(cell.region._node_id)
        assert raw[0] == 0
        assert raw[1] == 1
        assert raw[2] == 1  # positive sense

    def test_intersection(self, pin_cell_model):
        """Cell 2 = +1 -2 (intersection of two halfspaces)."""
        cell = pin_cell_model._cells[2]
        raw = pin_cell_model._sys.node_tree(cell.region._node_id)
        # Should be an intersection: (2, left, right)
        assert raw[0] == 2  # CSG_OP_INTERSECTION

    def test_nested_structure(self, pin_cell_model):
        """Cell 4 = +3 -4 +5 -6 (intersection of four halfspaces)."""
        cell = pin_cell_model._cells[4]
        raw = pin_cell_model._sys.node_tree(cell.region._node_id)
        # Should be a tree of intersections
        assert raw[0] == 2  # root is intersection


# ---------------------------------------------------------------------------
# _ImportedRegion.tree property
# ---------------------------------------------------------------------------

class TestImportedRegionTree:
    def test_tree_returns_region(self, sphere_model):
        """tree property should return a Python Region."""
        from aleathor.geometry import Region
        cell = sphere_model._cells[1]
        tree = cell.region.tree
        assert isinstance(tree, Region)

    def test_tree_caching(self, sphere_model):
        """tree property should cache the result."""
        cell = sphere_model._cells[1]
        tree1 = cell.region.tree
        tree2 = cell.region.tree
        assert tree1 is tree2

    def test_simple_halfspace(self, sphere_model):
        """Cell 1 = -1 should produce a Halfspace."""
        from aleathor.geometry import Halfspace
        cell = sphere_model._cells[1]
        tree = cell.region.tree
        assert isinstance(tree, Halfspace)
        assert tree.positive is False  # negative sense

    def test_intersection_type(self, pin_cell_model):
        """Cell 2 = +1 -2 should produce an Intersection."""
        from aleathor.geometry import Intersection
        cell = pin_cell_model._cells[2]
        tree = cell.region.tree
        assert isinstance(tree, Intersection)


# ---------------------------------------------------------------------------
# _ImportedRegion.__repr__
# ---------------------------------------------------------------------------

class TestImportedRegionRepr:
    def test_simple_negative(self, sphere_model):
        """Cell 1 = -1 should repr as '-1'."""
        cell = sphere_model._cells[1]
        assert repr(cell.region) == "-1"

    def test_simple_positive(self, sphere_model):
        """Cell 2 = +1 should repr as '+1'."""
        cell = sphere_model._cells[2]
        assert repr(cell.region) == "+1"

    def test_intersection_repr(self, pin_cell_model):
        """Cell 2 = +1 -2 should repr with & operator."""
        cell = pin_cell_model._cells[2]
        r = repr(cell.region)
        # Should contain both surface references and & operator
        assert "&" in r
        assert "1" in r
        assert "2" in r


# ---------------------------------------------------------------------------
# _ImportedRegion.get_surfaces()
# ---------------------------------------------------------------------------

class TestImportedRegionGetSurfaces:
    def test_single_surface(self, sphere_model):
        """Cell 1 = -1 should reference surface 1."""
        cell = sphere_model._cells[1]
        surfs = cell.region.get_surfaces()
        ids = {s.id for s in surfs}
        assert ids == {1}

    def test_two_surfaces(self, pin_cell_model):
        """Cell 2 = +1 -2 should reference surfaces 1 and 2."""
        cell = pin_cell_model._cells[2]
        surfs = cell.region.get_surfaces()
        ids = {s.id for s in surfs}
        assert ids == {1, 2}

    def test_four_surfaces(self, pin_cell_model):
        """Cell 4 = +3 -4 +5 -6 should reference surfaces 3, 4, 5, 6."""
        cell = pin_cell_model._cells[4]
        surfs = cell.region.get_surfaces()
        ids = {s.id for s in surfs}
        assert ids == {3, 4, 5, 6}


# ---------------------------------------------------------------------------
# __contains__ still works
# ---------------------------------------------------------------------------

class TestImportedRegionContains:
    def test_point_inside_sphere(self, sphere_model):
        """Origin should be inside cell 1 (sphere r=5)."""
        cell = sphere_model._cells[1]
        assert (0, 0, 0) in cell.region

    def test_point_outside_sphere(self, sphere_model):
        """Far point should be inside cell 2 (void outside)."""
        cell = sphere_model._cells[2]
        assert (100, 0, 0) in cell.region

    def test_point_not_inside(self, sphere_model):
        """Far point should NOT be inside cell 1."""
        cell = sphere_model._cells[1]
        assert (100, 0, 0) not in cell.region
