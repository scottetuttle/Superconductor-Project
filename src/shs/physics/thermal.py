"""
SHS Thermal Physics Module.

Defines thermal physics models used by SHS.

Material properties are provided by MaterialMap.
This module only defines thermal behavior.
"""

from dataclasses import dataclass


@dataclass
class ThermalModel:
    """
    Thermal physics configuration.

    The physical thermal timestep supplied to the solver
    may be subdivided internally using max_substep.
    """

    bath_temperature: float

    thermal_relaxation_rate: float = 0.0

    max_substep: float = 1e-13

    def validate(self):
        """
        Validate the thermal model parameters.
        """

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