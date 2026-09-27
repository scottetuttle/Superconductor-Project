"""
SHS Thermal Physics Module.

Defines thermal physics models used by SHS.

Material properties are provided by MaterialMap.
This module only defines thermal behavior.
"""

from dataclasses import dataclass
from math import isfinite

from shs.utils.defaults import default_section


_DEFAULTS = default_section("thermal")


@dataclass
class ThermalModel:
    """
    Thermal physics configuration.

    The physical thermal timestep supplied to the solver
    may be subdivided internally using max_substep.
    """

    bath_temperature: float

    thermal_relaxation_rate: float = _DEFAULTS["thermal_relaxation_rate"]

    max_substep: float = _DEFAULTS["max_substep"]

    stability_safety_factor: float = _DEFAULTS["stability_safety_factor"]

    model: str = _DEFAULTS["model"]
    electron_heat_capacity_fraction: float = _DEFAULTS["electron_heat_capacity_fraction"]
    electron_thermal_conductivity_fraction: float = _DEFAULTS["electron_thermal_conductivity_fraction"]
    electron_phonon_coupling_W_m3_K: float = _DEFAULTS["electron_phonon_coupling_W_m3_K"]
    phonon_escape_rate_per_s: float = _DEFAULTS["phonon_escape_rate_per_s"]

    def validate(self):
        """
        Validate the thermal model parameters.
        """

        if self.model not in {"single_temperature", "two_temperature"}:
            raise ValueError("Unknown thermal model.")
        if not 0 < self.electron_heat_capacity_fraction < 1:
            raise ValueError("Electron heat-capacity fraction must lie in (0, 1).")
        if not 0 <= self.electron_thermal_conductivity_fraction <= 1:
            raise ValueError("Electron conductivity fraction must lie in [0, 1].")
        if (not all(map(isfinite, (
                self.electron_phonon_coupling_W_m3_K, self.phonon_escape_rate_per_s)))
                or self.electron_phonon_coupling_W_m3_K < 0
                or self.phonon_escape_rate_per_s < 0):
            raise ValueError("Electron-phonon coupling and phonon escape must be finite and nonnegative.")
        if self.bath_temperature < 0.0:
            raise ValueError(
                "Bath temperature must be non-negative."
            )

        if self.thermal_relaxation_rate < 0.0:
            raise ValueError(
                "Thermal relaxation rate must be non-negative."
            )

        if self.max_substep <= 0.0:
            raise ValueError(
                "Thermal max_substep must be positive."
            )
        if self.model == "two_temperature" and self.thermal_relaxation_rate:
            raise ValueError("Use phonon_escape_rate_per_s for two-temperature bath cooling.")

        if not 0.0 < self.stability_safety_factor <= 1.0:
            raise ValueError("Thermal stability_safety_factor must lie in (0, 1].")

        return True

    @staticmethod
    def thermal_diffusivity(
        thermal_conductivity,
        heat_capacity,
    ):
        """
        Calculate thermal diffusivity.

        alpha = k / C
        """

        return (
            thermal_conductivity /
            heat_capacity
        )
