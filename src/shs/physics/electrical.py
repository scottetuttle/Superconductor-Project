"""
SHS Electrical Physics Module

Defines electrical transport properties.

Current model:

    J = σE

and Joule heating

    Q = J · E

Future:

- superconducting conductivity
- two-fluid model
- nonlinear resistivity
- vortex dissipation
- TDGL coupling
"""

from dataclasses import dataclass

import numpy as np


@dataclass
class ElectricalModel:
    """
    Electrical transport model.
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

    def joule_heating(
        self,
        current_density_x,
        current_density_y,
        electric_field_x,
        electric_field_y,
    ):
        """
        Compute Joule heating.

        Q = J · E

        Returns
        -------
        ndarray
            Volumetric heat generation.
        """

        return (
            current_density_x * electric_field_x
            +
            current_density_y * electric_field_y
        )