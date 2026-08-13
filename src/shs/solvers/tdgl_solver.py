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
        Dimensionless TDGL timestep.

    tdgl_model:
        TDGL physical model.

    Returns
    -------
    Fields
        Updated simulation fields.

    Notes
    -----
    The normalized TDGL equation is:

        u dpsi/dt =
            D^2 psi
            + (1 - T/Tc) psi
            - |psi|^2 psi

    Spatial coordinates are normalized by the
    coherence length xi.
    """

    if dt <= 0.0:
        raise ValueError(
            "TDGL timestep dt must be positive."
        )

    fields = simulation.fields
    material_map = simulation.material_map
    mesh = simulation.mesh

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
    # TDGL relaxation parameter.
    #
    # This belongs to the dimensionless TDGL model,
    # not to the spatial material map.
    #

    u = tdgl_model.parameters.u

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
    # Numerical sanity checks.
    #

    if not np.all(np.isfinite(dpsi_dt)):
        raise RuntimeError(
            "TDGL produced non-finite derivative values."
        )

    #
    # Explicit Euler integration.
    #

    fields.psi = (
        psi +
        dt * dpsi_dt
    )
    #
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