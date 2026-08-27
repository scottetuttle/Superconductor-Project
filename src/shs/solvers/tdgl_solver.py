"""
TDGL Solver.

Advances the superconducting order parameter using
the normalized Time-Dependent Ginzburg-Landau equation:

    u dpsi/dt =
        D^2 psi
        + (1 - T/Tc) psi
        - |psi|^2 psi

where:

    D = nabla - i A

The current implementation uses explicit Euler
time integration.

The TDGL parameter u is supplied by TDGLParameters.
Material properties provide the physical scales needed
to construct the dimensionless equation.
"""

import numpy as np

from shs.config.simulation_state import Simulation

from shs.tdgl.model import TDGLModel
from shs.tdgl.operators import covariant_laplacian

from shs.tdgl.boundary import (
    apply_insulating_boundary,
)


def tdgl_step(
    simulation,
    dt,
    tdgl_model,
):
    """
    Advance the TDGL order parameter by one physical timestep.

    Parameters
    ----------
    dt:
        Physical timestep in seconds.

    The timestep is converted internally to the
    dimensionless TDGL timescale.
    """

    if dt <= 0.0:
        raise ValueError(
            "TDGL timestep dt must be positive."
        )

    fields = simulation.fields
    material_map = simulation.material_map
    mesh = simulation.mesh

    material = material_map.materials[0]

    Tc = material.Tc

    dt_dimensionless = (
        tdgl_model.dimensional_to_normalized_time(
            dt,
            Tc,
        )
    )

    if not np.all(np.isfinite(dt_dimensionless)):
        raise RuntimeError(
            "TDGL produced a non-finite normalized timestep."
        )

    max_dt = (
        tdgl_model.parameters.max_normalized_timestep
    )
    tau_GL = tdgl_model.characteristic_time(Tc)

    if max_dt <= 0.0:
        raise ValueError(
            "Maximum normalized TDGL timestep must be positive."
        )

    # Use the largest local normalized timestep so that
    # the entire spatial domain is advanced with a stable
    # explicit timestep.

    dt_dimensionless = np.max(
        dt_dimensionless
    )

    number_of_steps = max(
        1,
        int(
            np.ceil(
                dt_dimensionless /
                max_dt
            )
        )
    )

    sub_dt = (
        dt_dimensionless /
        number_of_steps
    )

    psi = fields.psi.copy()

    #
    # Validate TDGL parameters.
    #

    tdgl_model.parameters.validate()

    #
    # Reduced temperature.
    #
    # t = T / Tc
    #

    reduced_temperature = (
        fields.temperature /
        material_map.Tc
    )

    #
    # Electromagnetic vector potential.
    #

    Ax = fields.vector_potential_x
    Ay = fields.vector_potential_y

    #
    # Normalize spatial coordinates by coherence length.
    #
    # x' = x / xi
    #
    # Therefore:
    #
    # dx' = dx / xi
    #

    material = material_map.materials[0]

    xi = material.coherence_length

    if xi <= 0.0:
        raise ValueError(
            "Coherence length xi must be positive."
        )

    dx_dimensionless = (
        mesh.dx / xi
    )

    dy_dimensionless = (
        mesh.dy / xi
    )

    #
    # Gauge-covariant kinetic term.
    #

    #
    # TDGL relaxation parameter.
    #

    u = tdgl_model.parameters.u


    #
    # Explicit Euler TDGL integration.
    #
    # The physical timestep has been converted into
    # normalized TDGL time and divided into stable
    # substeps.
    #

    for _ in range(number_of_steps):

        #
        # Gauge-covariant kinetic term.
        #

        kinetic_term = covariant_laplacian(
            psi,
            Ax,
            Ay,
            dx_dimensionless,
            dy_dimensionless,
        )


        #
        # Linear GL contribution.
        #
        # (1 - T/Tc) psi
        #

        linear_term = (
            1.0 -
            reduced_temperature
        ) * psi


        #
        # Nonlinear GL contribution.
        #
        # |psi|^2 psi
        #

        nonlinear_term = (
            np.abs(psi)**2 *
            psi
        )


        #
        # TDGL evolution equation.
        #

        dpsi_dt = (
            kinetic_term
            +
            linear_term
            -
            nonlinear_term
        ) / u


        #
        # Numerical sanity check.
        #

        if not np.all(np.isfinite(dpsi_dt)):
            raise RuntimeError(
                "TDGL produced non-finite derivative values."
            )


        #
        # Explicit Euler substep.
        #

        psi = (
            psi +
            sub_dt * dpsi_dt
        )


    #
    # Store completed physical timestep.
    #

    fields.psi = psi
        #

    fields.psi = apply_insulating_boundary(
        fields.psi,
        Ax,
        Ay,
        dx_dimensionless,
        dy_dimensionless,
    )
# Supercurrent density.
#
# Calculate the superconducting current
# associated with the updated order parameter.
#

    supercurrent_x, supercurrent_y = (
        tdgl_model.supercurrent_density(
            fields.psi,
            Ax,
            Ay,
            dx_dimensionless,
            dy_dimensionless,
        )
    )

    fields.supercurrent_density_x = (
        supercurrent_x
    )

    fields.supercurrent_density_y = (
        supercurrent_y
    )
    #
    # Numerical sanity check after update.
    #

    if not np.all(np.isfinite(fields.psi)):
        raise RuntimeError(
            "TDGL produced non-finite order-parameter values."
        )

    return fields