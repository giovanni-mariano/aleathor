# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""Tests for Material API (composition, density, nuclides, elements, mixtures)."""

import pytest
import aleathor as ath


@pytest.fixture
def model():
    return ath.Model("Material Test")


# ---------------------------------------------------------------------------
# add_material basics
# ---------------------------------------------------------------------------

class TestAddMaterial:
    def test_returns_material(self, model):
        mat = model.add_material(1)
        assert isinstance(mat, ath.Material)

    def test_material_id(self, model):
        mat = model.add_material(42)
        assert mat.id == 42

    def test_material_name(self, model):
        mat = model.add_material(1, name="Steel")
        assert mat.name == "Steel"

    def test_material_name_none(self, model):
        mat = model.add_material(1)
        assert mat.name is None

    def test_duplicate_raises(self, model):
        model.add_material(1)
        with pytest.raises(ValueError, match="already exists"):
            model.add_material(1)

    def test_invalid_id_raises(self, model):
        with pytest.raises(ValueError, match="must be positive"):
            model.add_material(0)
        with pytest.raises(ValueError, match="must be positive"):
            model.add_material(-5)

    def test_density_kwarg(self, model):
        mat = model.add_material(1, density=7.8)
        assert mat.density == pytest.approx(7.8)


# ---------------------------------------------------------------------------
# Density
# ---------------------------------------------------------------------------

class TestMaterialDensity:
    def test_set_and_get(self, model):
        mat = model.add_material(1)
        mat.density = 10.97
        assert mat.density == pytest.approx(10.97)

    def test_no_density_returns_none(self, model):
        mat = model.add_material(1)
        assert mat.density is None


# ---------------------------------------------------------------------------
# Weight / atom fractions
# ---------------------------------------------------------------------------

class TestFractionType:
    def test_default_is_atom(self, model):
        mat = model.add_material(1)
        assert mat.weight_fractions is False

    def test_set_weight(self, model):
        mat = model.add_material(1)
        mat.weight_fractions = True
        assert mat.weight_fractions is True

    def test_set_atom(self, model):
        mat = model.add_material(1)
        mat.weight_fractions = True
        mat.weight_fractions = False
        assert mat.weight_fractions is False


# ---------------------------------------------------------------------------
# Nuclides
# ---------------------------------------------------------------------------

class TestNuclides:
    def test_add_nuclide(self, model):
        mat = model.add_material(1)
        mat.add_nuclide(92235, 0.05)
        assert len(mat.nuclides) == 1

    def test_nuclide_data(self, model):
        mat = model.add_material(1)
        mat.add_nuclide(92235, 0.05, ".80c")
        nuc = mat.nuclides[0]
        assert nuc['zaid'] == 92235
        assert nuc['fraction'] == pytest.approx(0.05)
        assert nuc['library'] == ".80c"

    def test_multiple_nuclides(self, model):
        mat = model.add_material(1)
        mat.add_nuclide(92235, 0.05)
        mat.add_nuclide(92238, 0.95)
        assert len(mat.nuclides) == 2

    def test_nuclide_no_library(self, model):
        mat = model.add_material(1)
        mat.add_nuclide(92235, 1.0)
        assert mat.nuclides[0]['library'] == ""

    def test_chaining(self, model):
        mat = model.add_material(1)
        result = mat.add_nuclide(92235, 0.05).add_nuclide(92238, 0.95)
        assert result is mat
        assert len(mat.nuclides) == 2


# ---------------------------------------------------------------------------
# Elements
# ---------------------------------------------------------------------------

class TestElements:
    def test_add_element(self, model):
        mat = model.add_material(1)
        mat.add_element(26, 1.0)
        assert len(mat.elements) == 1

    def test_element_data(self, model):
        mat = model.add_material(1)
        mat.add_element(26, 0.7)
        elem = mat.elements[0]
        assert elem['Z'] == 26
        assert elem['fraction'] == pytest.approx(0.7)

    def test_multiple_elements(self, model):
        mat = model.add_material(1)
        mat.add_element(26, 0.70)
        mat.add_element(24, 0.18)
        mat.add_element(28, 0.12)
        assert len(mat.elements) == 3

    def test_chaining(self, model):
        mat = model.add_material(1)
        result = mat.add_element(26, 0.7).add_element(24, 0.3)
        assert result is mat
        assert len(mat.elements) == 2


# ---------------------------------------------------------------------------
# Expand elements
# ---------------------------------------------------------------------------

class TestExpandElements:
    def test_expand_creates_nuclides(self, model):
        mat = model.add_material(1)
        mat.add_element(26, 1.0)  # Iron
        assert len(mat.nuclides) == 0
        mat.expand_elements()
        assert len(mat.nuclides) > 0  # Fe has 4 stable isotopes

    def test_expand_returns_self(self, model):
        mat = model.add_material(1)
        mat.add_element(8, 1.0)
        assert mat.expand_elements() is mat

    def test_expand_iron_has_fe56(self, model):
        mat = model.add_material(1)
        mat.add_element(26, 1.0)
        mat.expand_elements()
        zaids = [n['zaid'] for n in mat.nuclides]
        assert 26056 in zaids  # Fe-56 is most abundant


# ---------------------------------------------------------------------------
# Material queries on Model
# ---------------------------------------------------------------------------

class TestModelMaterialQueries:
    def test_get_material(self, model):
        model.add_material(5, name="Water")
        mat = model.get_material(5)
        assert mat.id == 5
        assert mat.name == "Water"

    def test_get_material_not_found(self, model):
        with pytest.raises(KeyError):
            model.get_material(999)

    def test_materials_list(self, model):
        model.add_material(1, name="A")
        model.add_material(2, name="B")
        mats = model.materials
        assert len(mats) >= 2
        ids = {m.id for m in mats}
        assert {1, 2}.issubset(ids)

    def test_materials_empty(self, model):
        assert len(model.materials) == 0


# ---------------------------------------------------------------------------
# Mixtures
# ---------------------------------------------------------------------------

class TestMixtures:
    def test_create_mixture_returns_id(self, model):
        model.add_material(1)
        model.add_material(2)
        mix_id = model.create_mixture([1, 2], [0.7, 0.3], new_id=10)
        assert mix_id == 10

    def test_mixture_auto_id(self, model):
        model.add_material(1)
        model.add_material(2)
        mix_id = model.create_mixture([1, 2], [0.5, 0.5])
        assert isinstance(mix_id, int)
        assert mix_id > 0

    def test_mixture_name_stored(self, model):
        model.add_material(1)
        model.add_material(2)
        model.create_mixture([1, 2], [0.5, 0.5], new_id=10, name="mix12")
        assert model._mixture_names[10] == "mix12"


# ---------------------------------------------------------------------------
# Material repr
# ---------------------------------------------------------------------------

class TestMaterialRepr:
    def test_basic_repr(self, model):
        mat = model.add_material(1)
        r = repr(mat)
        assert "Material(1" in r

    def test_repr_with_name(self, model):
        mat = model.add_material(1, name="Steel")
        r = repr(mat)
        assert "Steel" in r

    def test_repr_with_density(self, model):
        mat = model.add_material(1, density=7.8)
        r = repr(mat)
        assert "7.8" in r

    def test_repr_with_elements(self, model):
        mat = model.add_material(1)
        mat.add_element(26, 1.0)
        r = repr(mat)
        assert "1 elements" in r


# ---------------------------------------------------------------------------
# Integration: material used in cells
# ---------------------------------------------------------------------------

class TestMaterialCellIntegration:
    def test_cell_with_material(self, model):
        model.add_material(1, name="UO2", density=10.97)
        s = ath.Sphere(0, 0, 0, radius=5)
        cell = model.add_cell(region=-s, material=1, density=10.97)
        assert cell.material == 1

    def test_loaded_model_materials(self):
        """Materials from loaded MCNP models should be queryable."""
        from pathlib import Path
        data = Path(__file__).parent / "data" / "sphere_simple.inp"
        if not data.exists():
            pytest.skip("test data not found")
        m = ath.read_mcnp(data)
        mats = m.materials
        # Should have at least the materials referenced by cells
        assert isinstance(mats, list)
