"""
TDGL Solver

Advances the superconducting order parameter.

Dimensionless TDGL equation:

dpsi/dt = (1/u) [
    D^2 psi
    +
    (1 - T/Tc) psi
    -
    |psi|^2 psi
]

where:

D = nabla - iA

"""

import numpy as np

from shs.config.simulation_state import Simulation

from shs.tdgl.model import TDGLModel

from shs.tdgl.operators import covariant_laplacian



def tdgl_step(
    simulation: Simulation,
    dt: float,
    tdgl_model: TDGLModel,
):
    """
    Advance the TDGL order parameter by one timestep.

    Parameters
    ----------
    simulation:
        Complete SHS simulation state.

    dt:
        Time step.

    tdgl_model:
        TDGL parameters.

    Returns
    -------
    Fields
        Updated fields.
    """


    fields = simulation.fields

    material_map = simulation.material_map

    mesh = simulation.mesh



    psi = fields.psi.copy()



    #
    # Reduced temperature
    #
    # T/Tc
    #

    reduced_temperature = (
        fields.temperature /
        material_map.Tc
    )



    #
    # Electromagnetic coupling
    #

    Ax = fields.vector_potential_x

    Ay = fields.vector_potential_y



    kinetic_term = covariant_laplacian(
        psi,
        Ax,
        Ay,
        mesh.dx,
        mesh.dy,
    )



    #
    # Ginzburg-Landau potential
    #

    potential_term = (

        (
            1.0 -
            reduced_temperature
        )
        *
        psi

        -

        (
            np.abs(psi)**2
            *
            psi
        )

    )



    #
    # TDGL evolution
    #

    dpsi_dt = (

        kinetic_term
        +
        potential_term

    ) / tdgl_model.parameters.u



    #
    # Euler time integration
    #

    fields.psi = (
        psi
        +
        dt *
        dpsi_dt
    )


    return fields