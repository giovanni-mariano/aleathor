# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""Tests for the nuclear data module."""

import pytest
import aleathor as ath
from aleathor.nucdata import (
    XsDir, Nuclide, NucMaterial, Multigroup,
    parse_zaid, reaction_classify,
)


class TestParseZaid:
    """Tests for ZAID parsing."""

    def test_continuous_neutron(self):
        result = parse_zaid("92235.80c")
        assert result["Z"] == 92
        assert result["A"] == 235
        assert result["metastable"] == 0
        assert result["type"] == "continuous_neutron"

    def test_photoatomic(self):
        result = parse_zaid("1000.12p")
        assert result["Z"] == 1
        assert result["A"] == 0
        assert result["type"] == "photoatomic"

    def test_invalid_zaid(self):
        with pytest.raises(ValueError):
            parse_zaid("invalid")


class TestReactionClassify:
    """Tests for reaction classification."""

    def test_elastic(self):
        assert reaction_classify(2) == "scatter"

    def test_fission(self):
        assert reaction_classify(18) == "multiply"

    def test_capture(self):
        assert reaction_classify(102) == "absorption"

    def test_n2n(self):
        assert reaction_classify(16) == "multiply"

    def test_inelastic(self):
        # MT=51 is inelastic level 1
        assert reaction_classify(51) == "scatter"


class TestXsDir:
    """Tests for XsDir (require nuclear data files, skip if unavailable)."""

    @pytest.fixture
    def xsdir_path(self):
        """Try to find an xsdir file. Skip test if not available."""
        import os
        paths = [
            os.environ.get("DATAPATH", ""),
            os.path.expanduser("~/nuclear_data/xsdir"),
            "/opt/nuclear_data/xsdir",
        ]
        for p in paths:
            if p and os.path.exists(p):
                return p
        pytest.skip("No xsdir file found (set DATAPATH env var)")

    def test_load(self, xsdir_path):
        xsdir = XsDir(xsdir_path)
        assert xsdir.count > 0

    def test_find(self, xsdir_path):
        xsdir = XsDir(xsdir_path)
        entry = xsdir.find("92235.80c")
        if entry is None:
            pytest.skip("92235.80c not in xsdir")
        assert entry["zaid"] == "92235.80c"
        assert entry["awr"] > 0

    def test_contains(self, xsdir_path):
        xsdir = XsDir(xsdir_path)
        # At least one nuclide should exist
        assert xsdir.count > 0

    def test_invalid_path(self):
        with pytest.raises(IOError):
            XsDir("/nonexistent/xsdir")


class TestNuclide:
    """Tests for Nuclide (require nuclear data files)."""

    @pytest.fixture
    def xsdir(self):
        import os
        path = os.environ.get("DATAPATH", "")
        if not path or not os.path.exists(path):
            pytest.skip("No xsdir file found (set DATAPATH env var)")
        return XsDir(path)

    @pytest.fixture
    def u235(self, xsdir):
        if xsdir.find("92235.80c") is None:
            pytest.skip("92235.80c not available")
        return Nuclide(xsdir, "92235.80c")

    def test_properties(self, u235):
        assert u235.Z == 92
        assert u235.A == 235
        assert u235.awr > 0
        assert u235.temperature >= 0
        assert u235.n_energies > 0
        assert u235.n_reactions > 0

    def test_xs_total(self, u235):
        sigma = u235.xs_total(1.0)
        assert sigma > 0

    def test_xs_absorption(self, u235):
        sigma = u235.xs_absorption(1.0)
        assert sigma > 0

    def test_xs_elastic(self, u235):
        sigma = u235.xs_elastic(1.0)
        assert sigma > 0

    def test_fissile(self, u235):
        assert u235.is_fissile
        nu = u235.nu_bar(1.0)
        assert nu > 2.0  # U-235 nu_bar ~ 2.4

    def test_energy_grid(self, u235):
        grid = u235.energy_grid()
        assert len(grid) == u235.n_energies
        # Should be ascending
        assert grid[0] < grid[-1]

    def test_reactions(self, u235):
        rxns = u235.reactions()
        assert len(rxns) == u235.n_reactions
        for rxn in rxns:
            assert "mt" in rxn
            assert "q_value" in rxn

    def test_repr(self, u235):
        r = repr(u235)
        assert "92235" in r

    def test_invalid_zaid(self, xsdir):
        with pytest.raises(ValueError):
            Nuclide(xsdir, "99999.00c")


class TestNucMaterial:
    """Tests for NucMaterial (require nuclear data files)."""

    @pytest.fixture
    def xsdir(self):
        import os
        path = os.environ.get("DATAPATH", "")
        if not path or not os.path.exists(path):
            pytest.skip("No xsdir file found (set DATAPATH env var)")
        return XsDir(path)

    def test_create(self):
        mat = NucMaterial()
        assert mat is not None

    def test_add_and_query(self, xsdir):
        if xsdir.find("92235.80c") is None:
            pytest.skip("92235.80c not available")
        u235 = Nuclide(xsdir, "92235.80c")
        mat = NucMaterial()
        mat.add(u235, 0.048)
        sigma_t = mat.xs_total(1.0)
        assert sigma_t > 0
        mfp = mat.mean_free_path(1.0)
        assert mfp > 0

    def test_chaining(self, xsdir):
        if xsdir.find("92235.80c") is None:
            pytest.skip("92235.80c not available")
        u235 = Nuclide(xsdir, "92235.80c")
        mat = NucMaterial()
        result = mat.add(u235, 0.048)
        assert result is mat


class TestMultigroup:
    """Tests for Multigroup (require nuclear data files)."""

    @pytest.fixture
    def xsdir(self):
        import os
        path = os.environ.get("DATAPATH", "")
        if not path or not os.path.exists(path):
            pytest.skip("No xsdir file found (set DATAPATH env var)")
        return XsDir(path)

    def test_create(self):
        mg = Multigroup([20.0, 1.0, 0.1, 1e-5])
        assert mg.n_groups == 3

    def test_invalid_bounds(self):
        with pytest.raises(ValueError):
            Multigroup([1.0])  # need at least 2

    def test_collapse_and_data(self, xsdir):
        if xsdir.find("92235.80c") is None:
            pytest.skip("92235.80c not available")
        u235 = Nuclide(xsdir, "92235.80c")
        mg = Multigroup([20.0, 1.0, 0.1, 1e-5])
        mg.collapse(u235)
        data = mg.get_data()
        assert data["n_groups"] == 3
        assert len(data["sigma_t"]) == 3
        assert len(data["sigma_a"]) == 3
        assert len(data["chi"]) == 3
        assert len(data["scatter_matrix"]) == 3
        assert len(data["scatter_matrix"][0]) == 3
        # All cross sections should be positive
        for s in data["sigma_t"]:
            assert s >= 0

    def test_repr(self):
        mg = Multigroup([20.0, 1.0, 1e-5])
        assert "2" in repr(mg)
