# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""
Nuclear data module for cross-section lookup and transport calculations.

Reads ACE-format nuclear data files (continuous-energy neutron and photoatomic)
via the xsdir/xsdata directory system.

Example:
    import aleathor as ath
    from aleathor.nucdata import XsDir, Nuclide, NucMaterial, Multigroup

    # Load cross-section directory
    xsdir = XsDir("/path/to/xsdir")

    # Load a nuclide and query cross sections
    u235 = Nuclide(xsdir, "92235.80c")
    sigma_t = u235.xs_total(1.0)   # 1 MeV -> total XS in barns

    # Build a material and get macroscopic cross sections
    mat = NucMaterial()
    mat.add(u235, number_density=0.048)
    mfp = mat.mean_free_path(1.0)  # mean free path in cm

    # Multigroup collapse
    bounds = [20.0, 1.0, 0.1, 1e-5]  # 3-group, descending MeV
    mg = Multigroup(bounds)
    mg.collapse(u235)
    data = mg.get_data()  # dict with sigma_t, sigma_a, scatter_matrix, etc.
"""

from __future__ import annotations
from typing import Optional, List, Dict, Any

try:
    import _alea
except ImportError:
    _alea = None


class XsDir:
    """Cross-section directory (xsdir/xsdata) for resolving nuclide data.

    Args:
        path: Path to xsdir file or directory of .xsd files
        directory: If True, treat path as a directory of per-nuclide .xsd files
                   (FENDL-style). Default False.

    Example:
        xsdir = XsDir("/path/to/xsdir")
        xsdir = XsDir("/path/to/fendl/neutron", directory=True)
    """

    def __init__(self, path: str, *, directory: bool = False):
        if _alea is None:
            raise RuntimeError("Native extension not available")
        self._xsdir = _alea.XsDir(path, directory=directory)

    @property
    def count(self) -> int:
        """Number of entries in the directory."""
        return self._xsdir.count

    def find(self, zaid: str) -> Optional[Dict[str, Any]]:
        """Find an xsdir entry by ZAID string.

        Args:
            zaid: ZAID string, e.g. "92235.80c"

        Returns:
            Dict with zaid, awr, filename, file_type, address, temperature,
            or None if not found.
        """
        return self._xsdir.find(zaid)

    def __contains__(self, zaid: str) -> bool:
        return self._xsdir.find(zaid) is not None

    def __repr__(self) -> str:
        return f"XsDir(count={self.count})"


class Nuclide:
    """Decoded ACE nuclide with cross-section lookup.

    By default, nuclides are cached in the XsDir and shared across lookups.
    Use ``cached=False`` to get an independent copy that can be modified
    (e.g. for Doppler broadening).

    Args:
        xsdir: Cross-section directory
        zaid: ZAID string, e.g. "92235.80c"
        cached: If True (default), use cached nuclide from xsdir.
                If False, load a fresh independent copy.

    Example:
        u235 = Nuclide(xsdir, "92235.80c")
        print(u235.Z, u235.A)       # 92, 235
        print(u235.xs_total(1.0))   # total XS at 1 MeV in barns
    """

    def __init__(self, xsdir: XsDir, zaid: str, *, cached: bool = True):
        if _alea is None:
            raise RuntimeError("Native extension not available")
        self._nuc = _alea.Nuclide(xsdir._xsdir, zaid, cached=cached)

    @property
    def zaid(self) -> str:
        """ZAID string (e.g. '92235.80c')."""
        return self._nuc.zaid

    @property
    def Z(self) -> int:
        """Atomic number."""
        return self._nuc.Z

    @property
    def A(self) -> int:
        """Mass number."""
        return self._nuc.A

    @property
    def awr(self) -> float:
        """Atomic weight ratio."""
        return self._nuc.awr

    @property
    def temperature(self) -> float:
        """Temperature kT in MeV."""
        return self._nuc.temperature

    @property
    def is_fissile(self) -> bool:
        """True if nuclide has fission data."""
        return self._nuc.is_fissile

    @property
    def has_urr(self) -> bool:
        """True if nuclide has unresolved resonance region data."""
        return self._nuc.has_urr

    @property
    def has_photon(self) -> bool:
        """True if nuclide has photon interaction data."""
        return self._nuc.has_photon

    @property
    def n_reactions(self) -> int:
        """Number of non-elastic reactions."""
        return self._nuc.n_reactions

    @property
    def n_energies(self) -> int:
        """Number of energy grid points."""
        return self._nuc.n_energies

    def xs_total(self, energy: float) -> float:
        """Total cross section (barns) at energy (MeV)."""
        return self._nuc.xs_total(energy)

    def xs_absorption(self, energy: float) -> float:
        """Absorption cross section (barns) at energy (MeV)."""
        return self._nuc.xs_absorption(energy)

    def xs_elastic(self, energy: float) -> float:
        """Elastic scattering cross section (barns) at energy (MeV)."""
        return self._nuc.xs_elastic(energy)

    def xs_reaction(self, mt: int, energy: float) -> float:
        """Cross section (barns) for reaction MT at energy (MeV)."""
        return self._nuc.xs_reaction(mt, energy)

    def xs_heating(self, energy: float) -> float:
        """Heating number (MeV*barn) at energy (MeV)."""
        return self._nuc.xs_heating(energy)

    def heating_per_collision(self, energy: float) -> float:
        """Average energy deposited per collision (MeV)."""
        return self._nuc.heating_per_collision(energy)

    def nu_bar(self, energy: float) -> float:
        """Average neutrons per fission at energy (MeV).

        Raises ValueError if the nuclide is not fissile.
        """
        return self._nuc.nu_bar(energy)

    def doppler_broaden(self, kT_target: float) -> None:
        """Doppler-broaden cross sections in-place to target temperature.

        Args:
            kT_target: Target temperature in MeV (e.g. 2.53e-8 for 293.6 K)

        Only works on non-cached nuclides (``cached=False``).
        Can only broaden to a higher temperature than the current one.
        """
        self._nuc.doppler_broaden(kT_target)

    def urr_factors(self, energy: float, xi: float) -> Optional[Dict[str, float]]:
        """Get URR probability table cross-section factors.

        Args:
            energy: Incident energy (MeV)
            xi: Random number [0, 1) for band selection

        Returns:
            Dict with total, elastic, fission, capture, heating factors,
            or None if URR doesn't apply at this energy.
        """
        return self._nuc.urr_factors(energy, xi)

    def photon_xs_incoherent(self, energy: float) -> float:
        """Compton (incoherent) scattering cross section (barns)."""
        return self._nuc.photon_xs_incoherent(energy)

    def photon_xs_coherent(self, energy: float) -> float:
        """Rayleigh (coherent) scattering cross section (barns)."""
        return self._nuc.photon_xs_coherent(energy)

    def photon_xs_photoelectric(self, energy: float) -> float:
        """Photoelectric cross section (barns)."""
        return self._nuc.photon_xs_photoelectric(energy)

    def photon_xs_pair(self, energy: float) -> float:
        """Pair production cross section (barns)."""
        return self._nuc.photon_xs_pair(energy)

    def energy_grid(self) -> List[float]:
        """Get the energy grid (MeV, ascending)."""
        return self._nuc.energy_grid()

    def reactions(self) -> List[Dict[str, Any]]:
        """Get list of reactions with mt, q_value, ty."""
        return self._nuc.reactions()

    def reaction_yield(self, mt: int, energy: float) -> float:
        """Get neutron yield for reaction MT at energy (MeV)."""
        return self._nuc.reaction_yield(mt, energy)

    def __repr__(self) -> str:
        return f"Nuclide('{self.zaid}', Z={self.Z}, A={self.A})"


class NucMaterial:
    """Nuclear material composition for transport calculations.

    A material is a collection of nuclides with number densities,
    providing macroscopic cross-section lookup and collision sampling.

    Example:
        mat = NucMaterial()
        mat.add(u235, 0.048)   # 0.048 atoms/barn-cm
        mat.add(o16, 0.024)
        print(mat.mean_free_path(1.0))   # cm at 1 MeV
    """

    def __init__(self):
        if _alea is None:
            raise RuntimeError("Native extension not available")
        self._mat = _alea.NucMaterial()

    def add(self, nuclide: Nuclide, number_density: float) -> 'NucMaterial':
        """Add a nuclide component.

        Args:
            nuclide: Nuclide to add
            number_density: Number density in atoms/barn-cm

        Returns:
            self (for chaining)
        """
        self._mat.add(nuclide._nuc, number_density)
        return self

    def xs_total(self, energy: float) -> float:
        """Macroscopic total cross section (cm^-1) at energy (MeV)."""
        return self._mat.xs_total(energy)

    def xs_absorption(self, energy: float) -> float:
        """Macroscopic absorption cross section (cm^-1)."""
        return self._mat.xs_absorption(energy)

    def xs_elastic(self, energy: float) -> float:
        """Macroscopic elastic scattering cross section (cm^-1)."""
        return self._mat.xs_elastic(energy)

    def mean_free_path(self, energy: float) -> float:
        """Mean free path (cm) at energy (MeV)."""
        return self._mat.mean_free_path(energy)

    def sample_distance(self, energy: float, xi: float) -> float:
        """Sample distance to next collision (cm).

        Args:
            energy: Incident energy (MeV)
            xi: Random number [0, 1)
        """
        return self._mat.sample_distance(energy, xi)

    def sample_nuclide(self, energy: float, xi: float) -> tuple:
        """Sample which nuclide in the material is hit.

        Args:
            energy: Incident energy (MeV)
            xi: Random number [0, 1)

        Returns:
            Tuple of (component_index, zaid_string)
        """
        return self._mat.sample_nuclide(energy, xi)

    def __repr__(self) -> str:
        return "NucMaterial()"


class Multigroup:
    """Multigroup cross sections and scattering matrix.

    Collapses continuous-energy cross sections into group-averaged constants
    using a weighting spectrum (default: Maxwellian + 1/E + fission).

    Args:
        bounds: Group boundaries in descending order (MeV).
                E.g. [20.0, 1.0, 0.1, 1e-5] for 3 groups.

    Example:
        mg = Multigroup([20.0, 1.0, 0.1, 1e-5])
        mg.collapse(u235)
        data = mg.get_data()
        print(data['sigma_t'])    # total XS per group
        print(data['chi'])        # fission spectrum
    """

    def __init__(self, bounds: List[float]):
        if _alea is None:
            raise RuntimeError("Native extension not available")
        self._mg = _alea.Multigroup(bounds)

    @property
    def n_groups(self) -> int:
        """Number of energy groups."""
        return self._mg.n_groups

    def set_spectrum(self, fn) -> None:
        """Set a custom weighting spectrum for group collapse.

        The default spectrum is Maxwellian + 1/E + Watt fission. Use this
        to override it with any callable.

        Args:
            fn: A callable ``fn(energy_MeV) -> float`` returning the spectrum
                value at the given energy, or ``None`` to reset to the default.

        Example:
            # Flat (constant) weighting
            mg.set_spectrum(lambda e: 1.0)

            # 1/E everywhere
            mg.set_spectrum(lambda e: 1.0 / e)

            # Reset to default
            mg.set_spectrum(None)
        """
        self._mg.set_spectrum(fn)

    def collapse(self, nuclide: Nuclide) -> None:
        """Collapse continuous-energy cross sections into multigroup constants.

        Args:
            nuclide: Nuclide with pointwise cross sections
        """
        self._mg.collapse(nuclide._nuc)

    def scatter(self, g_from: int, g_to: int) -> float:
        """Forward scattering matrix element sigma_s(g_from -> g_to)."""
        return self._mg.scatter(g_from, g_to)

    def scatter_adjoint(self, g_from: int, g_to: int) -> float:
        """Adjoint scattering matrix element (transpose of forward)."""
        return self._mg.scatter_adjoint(g_from, g_to)

    def sample_scatter(self, g_from: int, xi: float, *, adjoint: bool = False) -> int:
        """Sample outgoing group from scattering.

        Args:
            g_from: Incoming group index
            xi: Random number [0, 1)
            adjoint: If True, use adjoint (transposed) matrix

        Returns:
            Outgoing group index
        """
        return self._mg.sample_scatter(g_from, xi, adjoint)

    def get_data(self) -> Dict[str, Any]:
        """Get all multigroup constants as a dictionary.

        Returns:
            Dict with keys: n_groups, bounds, sigma_t, sigma_a, sigma_s,
            sigma_f, nu_sigma_f, chi, scatter_matrix.
        """
        return self._mg.get_data()

    def __repr__(self) -> str:
        return f"Multigroup(n_groups={self.n_groups})"


def parse_zaid(zaid: str) -> Dict[str, Any]:
    """Parse a ZAID string into its components.

    Args:
        zaid: ZAID string, e.g. "92235.80c"

    Returns:
        Dict with Z (atomic number), A (mass number),
        metastable (0=ground), type (table type string).

    Example:
        >>> parse_zaid("92235.80c")
        {'Z': 92, 'A': 235, 'metastable': 0, 'type': 'continuous_neutron'}
    """
    if _alea is None:
        raise RuntimeError("Native extension not available")
    return _alea.parse_zaid(zaid)


def reaction_classify(mt: int) -> str:
    """Classify a reaction MT number.

    Args:
        mt: ENDF reaction MT number

    Returns:
        One of 'absorption', 'scatter', or 'multiply'.

    Example:
        >>> reaction_classify(18)   # fission
        'multiply'
        >>> reaction_classify(102)  # (n,gamma)
        'absorption'
    """
    if _alea is None:
        raise RuntimeError("Native extension not available")
    return _alea.reaction_classify(mt)
