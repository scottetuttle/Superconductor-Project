"""
SHS Thermal Physics Module

Defines thermal physics models used by SHS.

Material properties are provided by MaterialMap.
This module only defines thermal behavior.
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
    """

    bath_temperature: float

    thermal_relaxation_rate: float = 0.0


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