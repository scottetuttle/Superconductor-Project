"""
Thermal solver for SHS.

Solves:

dT/dt = alpha * laplacian(T) + Q/C - cooling

where:

alpha = k/C

k and C are spatial material properties supplied
by MaterialMap.

The solver operates directly on the Simulation object.
"""

import numpy as np

from shs.config.simulation_state import Simulation
from shs.physics.thermal import ThermalModel


def thermal_step(
    simulation: Simulation,
    dt: float,
    thermal_model: ThermalModel,
):
    """
    Advance the thermal state by one timestep.

    Parameters
    ----------
    simulation:
        Complete SHS simulation state.

    dt:
        Time step.

    thermal_model:
        Thermal physics configuration.

    Returns
    -------
    Fields
        Updated simulation fields.
    """

    mesh = simulation.mesh

    fields = simulation.fields

    material_map = simulation.material_map


    # Temperature field

    T = fields.temperature.copy()


    # Spatial material properties

    k = material_map.thermal_conductivity

    C = material_map.heat_capacity


    # Thermal diffusivity

    alpha = k / C


    # Heat source

    heat_source = fields.heat_source


    # Zero-gradient boundary padding

    padded = np.pad(
        T,
        pad_width=1,
        mode="edge"
    )


    laplacian = (

        padded[2:,1:-1] +
        padded[:-2,1:-1] +
        padded[1:-1,2:] +
        padded[1:-1,:-2] -

        4 * padded[1:-1,1:-1]

    ) / mesh.dx**2


    # Diffusion term

    update = (
        alpha *
        dt *
        laplacian
    )


    # External heating

    if heat_source is not None:

        update += (
            dt *
            heat_source /
            C
        )


    # Thermal bath coupling

    if thermal_model.thermal_relaxation_rate > 0:

        cooling = (
            thermal_model.thermal_relaxation_rate
            *
            (
                T -
                thermal_model.bath_temperature
            )
        )

        update -= dt * cooling


    # Update field

    fields.temperature = T + update


    return fields