"""
Thermal diffusion solver for SHS.

Solves the heat diffusion equation:

dT/dt = alpha * laplacian(T)

"""

import numpy as np

from shs.physics import ThermalModel
from shs.geometry import Mesh


def thermal_step(
    temperature: np.ndarray,
    mesh: Mesh,
    thermal_model: ThermalModel,
    dt: float
):
    """
    Advance temperature field by one timestep.

    Parameters
    ----------
    temperature:
        Current temperature array.

    mesh:
        Simulation mesh.

    thermal_model:
        Thermal properties.

    dt:
        Time step.

    Returns
    -------
    Updated temperature field.
    """

    alpha = thermal_model.thermal_diffusivity()

    new_temperature = temperature.copy()

    laplacian = (
        (
            np.roll(temperature, 1, axis=0)
            +
            np.roll(temperature, -1, axis=0)
            +
            np.roll(temperature, 1, axis=1)
            +
            np.roll(temperature, -1, axis=1)
            -
            4 * temperature
        )
        /
        mesh.dx**2
    )


    new_temperature += (
        alpha *
        dt *
        laplacian
    )


    return new_temperature