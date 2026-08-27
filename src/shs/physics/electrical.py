"""
SHS Electrical Physics Module

Defines electrical transport properties and coupling
between superconducting and normal current.
"""

from dataclasses import dataclass

import numpy as np


@dataclass
class ElectricalModel:
    """
    Electrical transport model.

    Separates normal and superconducting current:

        J = J_s + J_n

    with

        J_n = sigma_n E

    where the effective normal conductivity is reduced
    by the superconducting fraction.
    """

    reference_voltage: float = 1.0

    def conductivity(self, resistivity):
        """
        Convert resistivity to conductivity.
        """

        if np.any(resistivity == 0):
            raise ValueError(
                "Resistivity cannot be zero."
            )

        return 1.0 / resistivity

    def normal_conductivity(
        self,
        resistivity,
        superconducting_fraction,
    ):
        """
        Calculate the effective normal conductivity.

            sigma_eff = (1 - |psi|^2) sigma_n
        """

        conductivity = self.conductivity(
            resistivity
        )

        normal_fraction = np.maximum(
            1.0 - superconducting_fraction,
            0.0,
        )

        return (
            normal_fraction *
            conductivity
        )

    def normal_current(
        self,
        electric_field_x,
        electric_field_y,
        resistivity,
        superconducting_fraction,
    ):
        """
        Calculate the normal current contribution.

            J_n = sigma_eff E
        """

        conductivity = self.normal_conductivity(
            resistivity,
            superconducting_fraction,
        )

        current_x = (
            conductivity *
            electric_field_x
        )

        current_y = (
            conductivity *
            electric_field_y
        )

        return current_x, current_y

    def total_current(
        self,
        supercurrent_density_x,
        supercurrent_density_y,
        normal_current_x,
        normal_current_y,
    ):
        """
        Calculate total current density.

            J = J_s + J_n
        """

        total_x = (
            supercurrent_density_x
            +
            normal_current_x
        )

        total_y = (
            supercurrent_density_y
            +
            normal_current_y
        )

        return total_x, total_y

    def joule_heating(
        self,
        normal_current_x,
        normal_current_y,
        electric_field_x,
        electric_field_y,
    ):
        """
        Calculate dissipative Joule heating.

            Q_J = J_n · E
        """

        return (
            normal_current_x * electric_field_x
            +
            normal_current_y * electric_field_y
        )
