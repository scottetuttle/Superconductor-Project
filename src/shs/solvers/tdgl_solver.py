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

The implementation uses explicit Euler time integration.

Because the requested physical timestep may be larger than
the stable normalized TDGL timestep, the solver automatically
subdivides the requested timestep into internal TDGL substeps.
"""

import numpy as np

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
    Advance the TDGL state by one physical timestep.

    Parameters
    ----------
    simulation:
        Complete SHS simulation state.

    dt:
        Physical timestep in seconds.

    tdgl_model:
        TDGL model containing the dimensionless parameters.

    Returns
    -------
    Fields
        Updated simulation fields.
    """

    if dt <= 0.0:
        raise ValueError(
            "TDGL timestep dt must be positive."
        )

    fields = simulation.fields
    material_map = simulation.material_map
    mesh = simulation.mesh

    #
    # Validate TDGL parameters before integration.
    #

    tdgl_model.parameters.validate()

    #
    # Material scales.
    #

    material = material_map.materials[0]

    Tc = material.Tc
    xi = material.coherence_length

    if Tc <= 0.0:
        raise ValueError(
            "Critical temperature Tc must be positive."
        )

    if xi <= 0.0:
        raise ValueError(
            "Coherence length xi must be positive."
        )

    #
    # Convert physical time to normalized TDGL time.
    #

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

    dt_dimensionless = float(
        np.max(dt_dimensionless)
    )

    #
    # Maximum stable/allowed normalized timestep.
    #

    max_dt = (
        tdgl_model.parameters.max_normalized_timestep
    )

    if max_dt <= 0.0:
        raise ValueError(
            "Maximum normalized TDGL timestep must be positive."
        )

    #
    # Divide the requested physical timestep into
    # stable internal normalized substeps.
    #

    number_of_steps = max(
        1,
        int(
            np.ceil(
                dt_dimensionless /
                max_dt
            )
        ),
    )

    sub_dt = (
        dt_dimensionless /
        number_of_steps
    )

    #
    # Initial order parameter.
    #

    psi = fields.psi.copy()

    #
    # Reduced temperature.
    #

    reduced_temperature = (
        fields.temperature /
        material_map.Tc
    )

    #
    # Vector potential.
    #

    Ax = fields.vector_potential_x
    Ay = fields.vector_potential_y

    #
    # Normalize spatial coordinates by coherence length.
    #

    dx_dimensionless = (
        mesh.dx /
        xi
    )

    dy_dimensionless = (
        mesh.dy /
        xi
    )

    #
    # TDGL relaxation parameter.
    #

    u = tdgl_model.parameters.u

    #
    # Explicit Euler integration.
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

        linear_term = (
            1.0 -
            reduced_temperature
        ) * psi

        #
        # Nonlinear GL contribution.
        #

        nonlinear_term = (
            np.abs(psi) ** 2
            *
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
            sub_dt *
            dpsi_dt
        )

    #
    # Store completed physical timestep.
    #

    fields.psi = psi

    #
    # Apply superconducting boundary condition.
    #

    fields.psi = apply_insulating_boundary(
        fields.psi,
        Ax,
        Ay,
        dx_dimensionless,
        dy_dimensionless,
    )

    #
    # Calculate superconducting current.
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
    # Final numerical sanity check.
    #

    if not np.all(np.isfinite(fields.psi)):
        raise RuntimeError(
            "TDGL produced non-finite order-parameter values."
        )

    return fields