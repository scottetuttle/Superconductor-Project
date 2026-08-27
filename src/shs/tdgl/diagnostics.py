"""
TDGL diagnostic quantities.

Provides derived quantities from the superconducting
order parameter and electromagnetic fields.

These functions do not modify simulation state.
"""

import numpy as np

from shs.tdgl.operators import covariant_gradient


def order_parameter_amplitude(
    psi,
):
    """
    Calculate |psi|.
    """

    return np.abs(psi)


def order_parameter_amplitude_squared(
    psi,
):
    """
    Calculate |psi|^2.
    """

    return np.abs(psi) ** 2


def order_parameter_phase(
    psi,
):
    """
    Calculate the phase of the order parameter.

    Returns values in the range [-pi, pi].
    """

    return np.angle(psi)


def supercurrent_magnitude(
    supercurrent_density_x,
    supercurrent_density_y,
):
    """
    Calculate the magnitude of the superconducting
    current density.
    """

    return np.sqrt(
        supercurrent_density_x**2
        +
        supercurrent_density_y**2
    )


def covariant_gradient_magnitude(
    psi,
    vector_potential_x,
    vector_potential_y,
    dx,
    dy,
):
    """
    Calculate |D psi|.

    D = nabla - i A
    """

    Dx_psi, Dy_psi = covariant_gradient(
        psi,
        vector_potential_x,
        vector_potential_y,
        dx,
        dy,
    )

    return np.sqrt(
        np.abs(Dx_psi)**2
        +
        np.abs(Dy_psi)**2
    )


def free_energy_density(
    psi,
    reduced_temperature,
    vector_potential_x=None,
    vector_potential_y=None,
    dx=None,
    dy=None,
):
    """
    Calculate the dimensionless GL free-energy density.

    The local GL contribution is

        f_local =
            -(1 - T/Tc)|psi|^2
            + 1/2 |psi|^4

    If the vector potential and grid spacing are supplied,
    the kinetic contribution is included:

        f_kinetic =
            |D psi|^2

    Therefore:

        f =
            |D psi|^2
            -(1 - T/Tc)|psi|^2
            + 1/2 |psi|^4

    Returns
    -------
    ndarray
        Dimensionless free-energy density.
    """

    amplitude_squared = (
        np.abs(psi) ** 2
    )

    local_energy = (
        -(
            1.0 -
            reduced_temperature
        )
        * amplitude_squared
        +
        0.5 *
        amplitude_squared**2
    )

    if (
        vector_potential_x is None
        or vector_potential_y is None
    ):
        return local_energy

    if dx is None or dy is None:
        raise ValueError(
            "dx and dy must be supplied when "
            "including the kinetic contribution."
        )

    gradient_magnitude_squared = (
        covariant_gradient_magnitude(
            psi,
            vector_potential_x,
            vector_potential_y,
            dx,
            dy,
        ) ** 2
    )

    return (
        gradient_magnitude_squared
        +
        local_energy
    )


def total_free_energy(
    psi,
    reduced_temperature,
    vector_potential_x=None,
    vector_potential_y=None,
    dx=None,
    dy=None,
):
    """
    Calculate the total dimensionless GL free energy.

    The energy is obtained by summing the free-energy
    density over the simulation grid.
    """

    density = free_energy_density(
        psi,
        reduced_temperature,
        vector_potential_x,
        vector_potential_y,
        dx,
        dy,
    )

    return float(
        np.sum(density)
    )


def mean_order_parameter_amplitude(
    psi,
):
    """
    Calculate the spatial mean of |psi|.
    """

    return float(
        np.mean(
            np.abs(psi)
        )
    )


def maximum_order_parameter_amplitude(
    psi,
):
    """
    Calculate the maximum value of |psi|.
    """

    return float(
        np.max(
            np.abs(psi)
        )
    )


def minimum_order_parameter_amplitude(
    psi,
):
    """
    Calculate the minimum value of |psi|.
    """

    return float(
        np.min(
            np.abs(psi)
        )
    )