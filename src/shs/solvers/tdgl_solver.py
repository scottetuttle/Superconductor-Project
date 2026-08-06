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



    #
# Dimensionless spatial scaling
#

#
# Dimensionless TDGL scaling
#
# Current implementation assumes
# uniform material properties.
#

    xi = material_map.materials[0].coherence_length


    kinetic_term = covariant_laplacian(
        psi,
        Ax,
        Ay,
        mesh.dx / xi,
        mesh.dy / xi,
    )



    #
    # Ginzburg-Landau potential
    #
    alpha = (
        material_map.gl_alpha
        *
        (
            1 -
            reduced_temperature
        )
    )


    beta = material_map.gl_beta


    potential_term = (

        alpha *
        psi

        -

        beta *
        np.abs(psi)**2 *
        psi

    )



    #
    # TDGL evolution
    #

    dpsi_dt = (

        kinetic_term
        +
        potential_term

    ) / material_map.tdgl_u



    #
    # Euler time integration
    #

    fields.psi = (
        psi
        +
        dt *
        dpsi_dt
    )
    #
    # Numerical stability limiter
    #

    max_amplitude = 10.0


    amplitude = np.abs(
        fields.psi
    )


    mask = amplitude > max_amplitude


    if np.any(mask):

        fields.psi[mask] *= (
            max_amplitude /
            amplitude[mask]
        )

    return fields