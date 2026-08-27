"""
Superconducting electrical transport coupling.
"""

import numpy as np

from shs.physics.superconductivity import (
    SuperconductingTransportModel,
)


def superconducting_current_density(
    psi,
    supercurrent_density_x,
    supercurrent_density_y,
):
    """
    Return the superconducting current density.

    The TDGL model already calculates the dimensionless
    supercurrent from psi and the vector potential.

    This function provides a physical-transport interface
    for the coupled solver.
    """

    superconducting_fraction = (
        SuperconductingTransportModel
        .superconducting_fraction(psi)
    )

    current_x = (
        superconducting_fraction *
        supercurrent_density_x
    )

    current_y = (
        superconducting_fraction *
        supercurrent_density_y
    )

    return current_x, current_y