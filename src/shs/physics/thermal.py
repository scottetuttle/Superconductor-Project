"""
SHS Thermal Physics Module

Defines thermal physics models used by SHS.

Material properties are provided by MaterialMap.
This module only defines thermal behavior and numerical
controls specific to thermal time integration.
"""

from dataclasses import dataclass


@dataclass
class ThermalModel:
    """
    Defines thermal coupling behavior.

    Material properties such as:

    - thermal conductivity
    - heat capacity

    are stored in MaterialMap.

    Parameters
    ----------
    bath_temperature:
        Thermal bath temperature (K).

    thermal_relaxation_rate:
        Coupling strength to thermal bath (1/s).

    max_substep:
        Maximum physical timestep used by the explicit thermal
        integrator. A larger requested timestep is automatically
        divided into multiple internal substeps.
    """

    bath_temperature: float

    thermal_relaxation_rate: float = 0.0

    max_substep: float = 1e-13


    @staticmethod
    def thermal_diffusivity(
        thermal_conductivity,
        heat_capacity,
    ):
        """
        Calculate thermal diffusivity.

        alpha = k / C

        Parameters
        ----------
        thermal_conductivity:
            Thermal conductivity W/(m*K)

        heat_capacity:
            Volumetric heat capacity J/(m^3*K)

        Returns
        -------
        Thermal diffusivity m^2/s
        """

        return (
            thermal_conductivity /
            heat_capacity
        )


    def validate(self):
        """
        Validate thermal model parameters.

        Raises
        ------
        ValueError
            If a thermal parameter is invalid.
        """

        if self.bath_temperature < 0.0:
            raise ValueError(
                "Thermal bath temperature must be non-negative."
            )

        if self.thermal_relaxation_rate < 0.0:
            raise ValueError(
                "Thermal relaxation rate must be non-negative."
            )

        if self.max_substep <= 0.0:
            raise ValueError(
                "Maximum thermal substep must be positive."
            )

        return True