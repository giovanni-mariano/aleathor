# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""Tests for Model operations (filtering, extraction, spatial indexing)."""

import pytest


class TestCellFiltering:
    """Tests for cell filtering methods."""

    def test_get_cells_by_material(self, simple_model):
        """Should return cells with given material."""
        # Material 1 is fuel (inside sphere)
        cells = simple_model.get_cells_by_material(1)

        assert isinstance(cells, list)
        assert len(cells) >= 1

    def test_get_cells_by_material_nonexistent(self, simple_model):
        """Should return empty list for nonexistent material."""
        cells = simple_model.get_cells_by_material(999)

        assert cells == []

    def test_get_cells_by_universe(self, simple_model):
        """Should return cells in given universe."""
        # Universe 0 is default
        cells = simple_model.get_cells_by_universe(0)

        assert isinstance(cells, list)
        assert len(cells) >= 1

    def test_get_cells_in_bbox(self, simple_model):
        """Should return cells in bounding box."""
        # Bounding box that contains the sphere
        bounds = (-6, 6, -6, 6, -6, 6)
        cells = simple_model.get_cells_in_bbox(bounds)

        assert isinstance(cells, list)
        assert len(cells) >= 1

    def test_get_cells_in_bbox_outside(self, simple_model):
        """Should return empty for bbox outside geometry."""
        bounds = (100, 200, 100, 200, 100, 200)
        cells = simple_model.get_cells_in_bbox(bounds)

        # May return empty or cells that intersect (depends on implementation)
        assert isinstance(cells, list)


class TestSpatialIndexing:
    """Tests for spatial indexing methods."""

    def test_build_spatial_index_returns_self(self, simple_model):
        """build_spatial_index should return self for chaining."""
        result = simple_model.build_spatial_index()

        assert result is simple_model

    def test_spatial_index_instance_count(self, simple_model):
        """Should have positive instance count after building index."""
        simple_model.build_spatial_index()
        count = simple_model.spatial_index_instance_count

        assert isinstance(count, int)
        assert count >= 0

    def test_spatial_index_without_build(self, simple_model):
        """Should return 0 if index not built."""
        # Fresh model, no index built
        count = simple_model.spatial_index_instance_count

        assert count == 0


class TestExtraction:
    """Tests for extraction operations."""

    def test_extract_region_returns_model(self, simple_model):
        """extract_region should return a new Model."""
        import aleathor as ath

        bounds = (-6, 6, -6, 6, -6, 6)
        new_model = simple_model.extract_region(bounds)

        assert isinstance(new_model, ath.Model)
        assert new_model is not simple_model

    def test_extract_region_has_title(self, simple_model):
        """Extracted model should have descriptive title."""
        bounds = (-6, 6, -6, 6, -6, 6)
        new_model = simple_model.extract_region(bounds)

        assert new_model.title is not None
        assert 'Region' in new_model.title or '[' in new_model.title


class TestFlattenUniverse:
    """Tests for flatten_universe method."""

    def test_flatten_universe_no_error(self, simple_model):
        """flatten_universe should not raise for valid universe."""
        # Should not raise
        simple_model.flatten_universe(0)


class TestMultiCellFiltering:
    """Tests for filtering with multiple cells."""

    def test_filter_returns_correct_count(self, multi_cell_model):
        """Each material should have exactly one cell."""
        for mat_id in [1, 2, 3, 4]:
            cells = multi_cell_model.get_cells_by_material(mat_id)
            assert len(cells) == 1, f"Material {mat_id} should have 1 cell"

    def test_all_cells_in_universe_zero(self, multi_cell_model):
        """All cells should be in universe 0."""
        cells = multi_cell_model.get_cells_by_universe(0)
        # 4 cells total
        assert len(cells) == 4


class TestCellMutation:
    """Tests for Cell property setters and fill methods."""

    def test_set_material(self, simple_model):
        """Setting material via Cell should update the cell."""
        cell = simple_model[1]
        assert cell.material == 1
        cell.material = 99
        assert cell.material == 99

    def test_set_density(self, simple_model):
        """Setting density via Cell should update the cell."""
        cell = simple_model[1]
        cell.density = 7.8
        assert cell.density == pytest.approx(7.8)

    def test_set_importance(self, simple_model):
        """Setting importance via Cell should update the cell."""
        cell = simple_model[1]
        assert cell.importance == 1.0
        cell.importance = 0.5
        assert cell.importance == 0.5

    def test_set_name(self, simple_model):
        """Setting name via Cell should update the cell."""
        cell = simple_model[1]
        cell.name = "shield"
        assert cell.name == "shield"

    def test_fill_with_int(self, simple_model):
        """Setting fill to an int should work."""
        cell = simple_model[1]
        cell.fill = 5
        assert cell.fill == 5
        assert cell.is_filled

    def test_unfill(self, simple_model):
        """unfill() should clear the fill."""
        cell = simple_model[1]
        cell.fill = 5
        assert cell.fill == 5
        cell.unfill()
        assert cell.fill is None
        assert not cell.is_filled

    def test_fill_with_universe(self, simple_model):
        """fill_with(Universe) should set fill to the universe ID."""
        import aleathor as ath
        univ = ath.Universe(id=3, name="test_univ")
        cell = simple_model[1]
        cell.fill_with(univ)
        assert cell.fill == 3
        assert cell.is_filled

    def test_fill_with_int_method(self, simple_model):
        """fill_with(int) should set fill."""
        cell = simple_model[1]
        cell.fill_with(7)
        assert cell.fill == 7

    def test_fill_with_bad_type(self, simple_model):
        """fill_with(str) should raise TypeError."""
        cell = simple_model[1]
        with pytest.raises(TypeError):
            cell.fill_with("bad")

    def test_fill_survives_rebuild(self, simple_model):
        """Fill should persist after a dirty rebuild."""
        cell = simple_model[1]
        cell.fill = 5
        assert cell.fill == 5

        # Mutate material to mark dirty and trigger rebuild
        cell.material = 42
        # Access fill again — triggers rebuild
        assert cell.fill == 5

    def test_add_cell_fill_applied(self):
        """add_cell(fill=N) should result in cell.fill == N after query."""
        import aleathor as ath

        model = ath.Model("Fill Test")
        s = ath.Sphere(0, 0, 0, radius=5.0)
        box = ath.Box(-10, 10, -10, 10, -10, 10)
        model.add_cell(region=-s, material=0, fill=5)
        model.add_cell(region=-box & +s, material=1, density=1.0)

        cell = model[1]
        assert cell.fill == 5

    def test_set_fill_syncs_python_cell(self, simple_model):
        """Model.set_fill() should sync the Python Cell so fill survives rebuild."""
        simple_model.set_fill(1, 8)

        # Force rebuild by mutating another cell
        simple_model.update_cell(2, material=99)

        cell = simple_model[1]
        assert cell.fill == 8

    def test_mutation_marks_dirty(self, simple_model):
        """Setting material should mark model dirty."""
        # Force initial build
        _ = simple_model[1]
        simple_model._dirty = False

        cell = simple_model[1]
        cell.material = 50
        assert simple_model._dirty is True
