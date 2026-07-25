"""
Thermal diffusion solver for SHS.

Solves the heat diffusion equation:

dT/dt = alpha * laplacian(T)

"""

import numpy as np

from shs.physics import ThermalModel
from shs.geometry import Mesh
from shs.physics import Fields

def thermal_step(
    fields: Fields,
    mesh,
    thermal_model,
    dt
):
    """
    Advance temperature using thermal diffusion.

    Uses zero-flux boundary conditions.
    """

    alpha = thermal_model.thermal_diffusivity()

    T = fields.temperature.copy()

    heat_source = fields.heat_source

    # Pad array using edge values.
    # This creates zero-gradient boundaries.
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


    update = alpha * dt * laplacian


    if heat_source is not None:

        update += (
            dt *
            heat_source /
            thermal_model.heat_capacity
    )

    fields.temperature = T + update

    return fields