"""
Thermal solver for SHS.

Solves:

    dT/dt =
        alpha * laplacian(T)
        + Q/C
        - cooling

where:

    alpha = k/C

k and C are spatial material properties supplied
by MaterialMap.

The requested physical timestep is automatically
subdivided into internal thermal substeps when
necessary.
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
    Advance the thermal state by one physical timestep.

    The requested timestep is subdivided into internal
    thermal substeps according to thermal_model.max_substep.

    Parameters
    ----------
    simulation:
        Complete SHS simulation state.

    dt:
        Physical timestep.

    thermal_model:
        Thermal physics configuration.

    Returns
    -------
    Fields
        Updated simulation fields.
    """

    if dt <= 0.0:
        raise ValueError(
            "Thermal timestep dt must be positive."
        )

    thermal_model.validate()

    mesh = simulation.mesh

    fields = simulation.fields

    material_map = simulation.material_map

    #
    # Spatial material properties.
    #

    k = material_map.thermal_conductivity

    C = material_map.heat_capacity

    if np.any(C <= 0.0):
        raise ValueError(
            "Heat capacity must be positive."
        )

    #
    # Thermal diffusivity.
    #

    alpha = (
        thermal_model.thermal_diffusivity(
            k,
            C,
        )
    )

    #
    # Determine number of internal thermal substeps.
    #

    number_of_substeps = max(
        1,
        int(
            np.ceil(
                dt /
                thermal_model.max_substep
            )
        ),
    )

    sub_dt = (
        dt /
        number_of_substeps
    )

    #
    # External heat source.
    #
    # The heat source is evaluated from the simulation
    # state for this physical timestep.
    #

    heat_source = fields.heat_source

    #
    # Advance through internal thermal substeps.
    #

    T = fields.temperature.copy()

    for _ in range(number_of_substeps):

        #
        # Zero-gradient boundary padding.
        #

        padded = np.pad(
            T,
            pad_width=1,
            mode="edge",
        )

        #
        # Correct anisotropic finite-difference Laplacian.
        #
        # Array axis 1 corresponds to x.
        # Array axis 0 corresponds to y.
        #

        laplacian_x = (
            padded[1:-1, 2:]
            -
            2.0 * padded[1:-1, 1:-1]
            +
            padded[1:-1, :-2]
        ) / (
            mesh.dx ** 2
        )

        laplacian_y = (
            padded[2:, 1:-1]
            -
            2.0 * padded[1:-1, 1:-1]
            +
            padded[:-2, 1:-1]
        ) / (
            mesh.dy ** 2
        )

        laplacian = (
            laplacian_x +
            laplacian_y
        )

        #
        # Thermal diffusion.
        #

        update = (
            alpha *
            sub_dt *
            laplacian
        )

        #
        # External heating.
        #

        if heat_source is not None:

            update += (
                sub_dt *
                heat_source /
                C
            )

        #
        # Thermal bath coupling.
        #

        if (
            thermal_model.thermal_relaxation_rate
            > 0.0
        ):

            cooling = (
                thermal_model.thermal_relaxation_rate
                *
                (
                    T -
                    thermal_model.bath_temperature
                )
            )

            update -= (
                sub_dt *
                cooling
            )

        #
        # Explicit thermal update.
        #

        T = T + update

        #
        # Numerical sanity check.
        #

        if not np.all(np.isfinite(T)):
            raise RuntimeError(
                "Thermal solver produced "
                "non-finite temperature values."
            )

    #
    # Store completed physical timestep.
    #

    fields.temperature = T

    return fields