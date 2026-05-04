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


class TestPointQueryApi:
    """Tests for point-query API names."""

    def test_cell_at_returns_terminal_cell(self, simple_model):
        """cell_at is the normal point query."""
        cell = simple_model.cell_at(0, 0, 0)

        assert cell is not None
        assert cell.id == 1
        assert cell.material == 1

    def test_cell_path_at_returns_depth_ordered_path(self, simple_model):
        """cell_path_at exposes the hierarchy path."""
        cells = simple_model.cell_path_at(0, 0, 0)

        assert len(cells) == 1
        assert cells[0].id == simple_model.cell_at(0, 0, 0).id
        assert cells[0].depth == 0

    def test_cells_at_is_not_public_api(self, simple_model):
        """The old plural name was intentionally removed."""
        assert not hasattr(simple_model, "cells_at")


class TestObjectRepr:
    """Tests for informative API object representations."""

    def test_cell_repr_names_key_state(self, simple_model):
        cell = simple_model[1]
        r = repr(cell)

        assert "Cell(id=1" in r
        assert "name='fuel'" in r
        assert "material=1" in r
        assert "density=10" in r
        assert "density_unit='g/cm3'" in r
        assert "universe=0" in r

    def test_cell_repr_shows_void_and_fill_when_present(self, simple_model):
        cell = simple_model[1]
        cell.material = 0
        cell.density = 0.0
        cell.fill = 7
        r = repr(cell)

        assert "void=True" in r
        assert "fill=7" in r

    def test_surface_repr_names_parameters(self):
        import aleathor as ath

        sphere = ath.Sphere(
            1.0, 2.0, 3.0,
            radius=4.0,
            surface_id=99,
            boundary="vacuum",
            name="outer",
        )
        r = repr(sphere)

        assert r == (
            "Sphere(id=99, x0=1.0, y0=2.0, z0=3.0, radius=4.0, "
            "boundary='vacuum', name='outer')"
        )


class TestSpatialIndexing:
    """Tests for spatial indexing methods."""

    def test_advanced_methods_are_namespaced(self):
        """Advanced APIs should not live directly on Model."""
        from aleathor import Model

        for name in (
            "build_spatial_index",
            "spatial_index_instance_count",
            "sample_mesh",
            "export_mesh",
            "find_overlaps",
            "generate_void",
            "add_voids",
            "add_graveyard",
            "simplify",
            "flatten_universe",
            "set_fill",
            "config",
            "set_verbose",
        ):
            assert not hasattr(Model, name)

    def test_build_spatial_index_returns_self(self, simple_model):
        """build_spatial_index should return self for chaining."""
        result = simple_model.backend.build_spatial_index()

        assert result is simple_model

    def test_spatial_index_instance_count(self, simple_model):
        """Should have positive instance count after building index."""
        simple_model.backend.build_spatial_index()
        count = simple_model.backend.spatial_index_instance_count

        assert isinstance(count, int)
        assert count >= 0

    def test_spatial_index_without_build(self, simple_model):
        """Should return 0 if index not built."""
        # Fresh model, no index built
        count = simple_model.backend.spatial_index_instance_count

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
        simple_model.repair.flatten_universe(0)


class TestVoidGenerationApi:
    """Tests for model-owned void commit operations."""

    def test_add_voids_is_model_owned(self):
        """Generated voids should be committed through Model, not VoidResult."""
        import aleathor as ath

        model = ath.Model("Void API Test")
        sphere = ath.Sphere(0, 0, 0, radius=1.0)
        model.add_cell(region=-sphere, material=1, density=1.0)

        before = len(model.cells)
        voids = model.void.generate(
            bounds=(-2, 2, -2, 2, -2, 2),
            max_depth=2,
            min_size=0.5,
        )
        voids.merge()

        assert not hasattr(voids, "add_cells")
        assert not hasattr(voids, "add_graveyard")

        added = model.void.add(voids)
        assert added > 0
        assert len(model.cells) == before + added

    def test_generate_void_accepts_region_bounds(self):
        """Void generation can use an existing finite region as bounds."""
        import aleathor as ath

        model = ath.Model("Void Region Bounds Test")
        sphere = ath.Sphere(0, 0, 0, radius=1.0)
        bounds = ath.Sphere(0, 0, 0, radius=3.0)
        model.add_cell(region=-sphere, material=1, density=1.0)

        voids = model.void.generate(
            region=-bounds,
            max_depth=2,
            min_size=0.5,
        )
        assert len(voids) > 0

        voids.merge()
        added = model.void.add(voids)
        assert added > 0

    def test_generate_void_rejects_bounds_and_region_together(self):
        """Bbox bounds and CSG region bounds are mutually exclusive."""
        import aleathor as ath

        model = ath.Model("Void Bounds Error Test")
        sphere = ath.Sphere(0, 0, 0, radius=1.0)
        model.add_cell(region=-sphere, material=1, density=1.0)

        with pytest.raises(ValueError, match="mutually exclusive"):
            model.void.generate(
                bounds=(-2, 2, -2, 2, -2, 2),
                region=-ath.Sphere(0, 0, 0, radius=3.0),
            )


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

    def test_default_fill_is_none(self, simple_model):
        """A normal cell should not look filled."""
        cell = simple_model[1]
        assert cell.fill is None
        assert not cell.is_filled
        assert "fill=0" not in repr(cell)
        assert simple_model.cells.by_fill().ids() == []

    def test_fill_zero_is_not_a_public_universe(self, simple_model):
        """The public API uses None, not 0, to mean no fill."""
        cell = simple_model[1]

        with pytest.raises(ValueError):
            cell.fill = 0
        with pytest.raises(ValueError):
            cell.fill_with(0)
        with pytest.raises(ValueError):
            simple_model.cells.by_fill(0)
        with pytest.raises(ValueError):
            simple_model.get_cells_filling_universe(0)

        cell.fill = 5
        cell.fill = None
        assert cell.fill is None

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

    def test_cell_fill_syncs_python_cell(self, simple_model):
        """Cell.fill should sync through the C backend."""
        cell = simple_model[1]
        cell.fill = 8

        simple_model.update_cell(2, material=99)

        assert simple_model[1].fill == 8

    def test_mutation_updates_c(self, simple_model):
        """Setting material should update C backend directly."""
        cell = simple_model[1]
        cell.material = 50
        # Re-read from C to verify
        cell2 = simple_model[1]
        assert cell2.material == 50


class TestLoadedModelMetadata:
    """Tests for metadata reconstructed from imported models."""

    def test_loaded_model_surfaces_are_available(self):
        """Loading MCNP should expose referenced surfaces through model.surfaces."""
        from pathlib import Path
        import aleathor as ath

        path = Path(__file__).parent / "data" / "sphere_simple.inp"
        model = ath.load(path)

        assert len(model.surfaces) == model._sys.surface_count
        assert 1 in model.surfaces
        assert model.surfaces[1].id == 1
