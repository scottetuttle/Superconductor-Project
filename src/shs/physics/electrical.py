"""
SHS Electrical Physics Module

Defines electrical transport properties and coupling
between superconducting and normal current.
"""

from dataclasses import dataclass

import numpy as np
from shs.utils.defaults import default_section


@dataclass
class ElectricalModel:
    """
    Electrical transport model.

    Separates normal and superconducting current:

        J = J_s + J_n

    with

        J_n = sigma_n E

    The standard generalized-TDGL model keeps the normal conductivity at its
    measured normal-state value. A condensate-depletion interpolation remains
    available as an explicitly selected legacy phenomenology.
    """

    reference_voltage: float = default_section("electromagnetic")["reference_voltage"]

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
        model="constant",
    ):
        """
        Calculate the effective normal conductivity.

            sigma_eff = (1 - |psi|^2) sigma_n
        """

        conductivity = self.conductivity(
            resistivity
        )

        if model == "constant":
            return np.broadcast_to(conductivity, np.shape(superconducting_fraction))
        if model == "condensate_depletion":
            return np.maximum(1.0 - superconducting_fraction, 0.0) * conductivity
        raise ValueError("Unknown normal conductivity model.")

    def normal_current(
        self,
        electric_field_x,
        electric_field_y,
        resistivity,
        superconducting_fraction,
        model="constant",
    ):
        """
        Calculate the normal current contribution.

            J_n = sigma_eff E
        """

        conductivity = self.normal_conductivity(
            resistivity,
            superconducting_fraction,
            model,
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
