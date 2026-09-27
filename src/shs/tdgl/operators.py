"""
TDGL numerical operators.

Contains differential operators acting on
the superconducting order parameter.

All operators preserve the shape of the
input field. Boundary values are currently
represented using zero derivatives for the
non-periodic outer boundary. More explicit
TDGL boundary conditions will be implemented
in the TDGL boundary module.
"""

import numpy as np


def gradient(
    psi,
    dx,
    dy,
):
    """
    Compute the spatial gradient of the complex
    superconducting order parameter.

    Parameters
    ----------
    psi : ndarray
        Complex superconducting order parameter.

    dx : float
        Grid spacing in x.

    dy : float
        Grid spacing in y.

    Returns
    -------
    dpsi_dx, dpsi_dy : ndarray
        Spatial derivatives with the same shape
        as psi.

    Notes
    -----
    Central differences are used in the interior.

    Boundary derivatives are currently set to zero.
    This is temporary infrastructure; explicit
    physical TDGL boundary conditions will later
    be applied by the TDGL boundary module.
    """

    dpsi_dx = np.zeros_like(
        psi,
        dtype=complex,
    )

    dpsi_dy = np.zeros_like(
        psi,
        dtype=complex,
    )

    if psi.shape[1] > 2:
        dpsi_dx[:, 1:-1] = (
            psi[:, 2:]
            -
            psi[:, :-2]
        ) / (
            2.0 * dx
        )

    if psi.shape[0] > 2:
        dpsi_dy[1:-1, :] = (
            psi[2:, :]
            -
            psi[:-2, :]
        ) / (
            2.0 * dy
        )

    return (
        dpsi_dx,
        dpsi_dy,
    )


def laplacian(
    psi,
    dx,
    dy,
):
    """
    Compute the Cartesian Laplacian.

    ∇²ψ =
        ∂²ψ/∂x² +
        ∂²ψ/∂y²

    A zero-normal-gradient (Neumann) boundary condition
    is used at the outer boundary.

    This corresponds to an insulating/no-flux boundary
    for the order parameter.

    The returned array always has the same shape as psi.
    """

    lap = np.zeros_like(
        psi,
        dtype=complex,
    )

    if psi.shape[0] < 2 or psi.shape[1] < 2:
        return lap

    #
    # Interior points
    #
    lap[1:-1, 1:-1] = (

        (
            psi[1:-1, 2:]
            -
            2.0 * psi[1:-1, 1:-1]
            +
            psi[1:-1, :-2]
        )
        / dx**2

        +

        (
            psi[2:, 1:-1]
            -
            2.0 * psi[1:-1, 1:-1]
            +
            psi[:-2, 1:-1]
        )
        / dy**2
    )

    #
    # Left boundary
    #
    # ∂psi/∂x = 0
    #
    # Implemented using a mirrored ghost cell:
    #
    # psi[-1] = psi[1]
    #
    lap[:, 0] = (
        2.0 * (
            psi[:, 1]
            -
            psi[:, 0]
        )
        / dx**2
    )

    #
    # Right boundary
    #
    lap[:, -1] = (
        2.0 * (
            psi[:, -2]
            -
            psi[:, -1]
        )
        / dx**2
    )

    #
    # Bottom boundary
    #
    lap[0, :] = (
        2.0 * (
            psi[1, :]
            -
            psi[0, :]
        )
        / dy**2
    )

    #
    # Top boundary
    #
    lap[-1, :] = (
        2.0 * (
            psi[-2, :]
            -
            psi[-1, :]
        )
        / dy**2
    )

    #
    # Corners receive contributions from both
    # coordinate directions.
    #
    lap[0, 0] = (
        2.0 * (
            psi[0, 1]
            -
            psi[0, 0]
        )
        / dx**2
        +
        2.0 * (
            psi[1, 0]
            -
            psi[0, 0]
        )
        / dy**2
    )

    lap[0, -1] = (
        2.0 * (
            psi[0, -2]
            -
            psi[0, -1]
        )
        / dx**2
        +
        2.0 * (
            psi[1, -1]
            -
            psi[0, -1]
        )
        / dy**2
    )

    lap[-1, 0] = (
        2.0 * (
            psi[-1, 1]
            -
            psi[-1, 0]
        )
        / dx**2
        +
        2.0 * (
            psi[-2, 0]
            -
            psi[-1, 0]
        )
        / dy**2
    )

    lap[-1, -1] = (
        2.0 * (
            psi[-1, -2]
            -
            psi[-1, -1]
        )
        / dx**2
        +
        2.0 * (
            psi[-2, -1]
            -
            psi[-1, -1]
        )
        / dy**2
    )

    return lap

def covariant_gradient(
    psi,
    vector_potential_x,
    vector_potential_y,
    dx,
    dy,
):
    """
    Compute the gauge-covariant gradient.

    D = ∇ - iA

    Therefore:

        Dx ψ = ∂ψ/∂x - i Ax ψ

        Dy ψ = ∂ψ/∂y - i Ay ψ

    All returned arrays preserve the shape of psi.
    """

    dpsi_dx, dpsi_dy = gradient(
        psi,
        dx,
        dy,
    )

    Dx = (
        dpsi_dx
        -
        1j * vector_potential_x * psi
    )

    Dy = (
        dpsi_dy
        -
        1j * vector_potential_y * psi
    )

    return (
        Dx,
        Dy,
    )


def covariant_laplacian(
    psi,
    vector_potential_x,
    vector_potential_y,
    dx,
    dy,
    active_mask=None,
):
    """
    Compute the gauge-consistent discrete covariant Laplacian.

    The continuum operator is

        D^2 psi

    with

        D = nabla - i A.

    The discrete operator uses gauge links:

        Ux = exp(-i Ax dx)
        Uy = exp(-i Ay dy)

    so that neighboring order-parameter values are
    parallel transported before forming the finite difference.

    This construction is gauge-consistent on the discrete grid
    and does not require the assumption

        div(A) = 0.

    Parameters
    ----------
    psi : ndarray
        Complex superconducting order parameter.

    vector_potential_x : ndarray
        x-component of dimensionless vector potential.

    vector_potential_y : ndarray
        y-component of dimensionless vector potential.

    dx : float
        Dimensionless grid spacing in x.

    dy : float
        Dimensionless grid spacing in y.

    Returns
    -------
    ndarray
        Gauge-consistent covariant Laplacian of psi.
    """

    if dx <= 0.0:
        raise ValueError(
            "dx must be positive."
        )

    if dy <= 0.0:
        raise ValueError(
            "dy must be positive."
        )

    if (
        psi.shape != vector_potential_x.shape
        or
        psi.shape != vector_potential_y.shape
    ):
        raise ValueError(
            "psi and vector potential fields "
            "must have the same shape."
        )

    #
    # Gauge links.
    #
    # These represent parallel transport between
    # neighboring grid points.
    #

    Ux = np.exp(
        -1j *
        vector_potential_x *
        dx
    )

    Uy = np.exp(
        -1j *
        vector_potential_y *
        dy
    )

    if active_mask is not None:
        active = np.asarray(active_mask, dtype=bool)
        if active.shape != psi.shape:
            raise ValueError("Active-domain mask must match psi.")
        lap = np.zeros_like(psi, dtype=complex)
        # Each active-active link contributes once to both endpoints. Omitting
        # links at the mask boundary is the finite-volume zero-normal-current
        # (covariant Neumann) condition for an insulating hole.
        x_links = active[:, :-1] & active[:, 1:]
        x_forward = (Ux[:, :-1] * psi[:, 1:] - psi[:, :-1]) / dx**2
        lap[:, :-1] += np.where(x_links, x_forward, 0.0)
        lap[:, 1:] += np.where(
            x_links,
            np.conjugate(Ux[:, :-1]) * psi[:, :-1] / dx**2 - psi[:, 1:] / dx**2,
            0.0,
        )
        y_links = active[:-1, :] & active[1:, :]
        y_forward = (Uy[:-1, :] * psi[1:, :] - psi[:-1, :]) / dy**2
        lap[:-1, :] += np.where(y_links, y_forward, 0.0)
        lap[1:, :] += np.where(
            y_links,
            np.conjugate(Uy[:-1, :]) * psi[:-1, :] / dy**2 - psi[1:, :] / dy**2,
            0.0,
        )
        lap[~active] = 0.0
        return lap

    lap = np.zeros_like(
        psi,
        dtype=complex,
    )

    #
    # Interior points.
    #
    # x direction
    #

    lap[:, 1:-1] += (
        Ux[:, 1:-1] * psi[:, 2:]
        +
        np.conjugate(
            Ux[:, :-2]
        ) * psi[:, :-2]
        -
        2.0 * psi[:, 1:-1]
    ) / dx**2

    #
    # y direction
    #

    lap[1:-1, :] += (
        Uy[1:-1, :] * psi[2:, :]
        +
        np.conjugate(
            Uy[:-2, :]
        ) * psi[:-2, :]
        -
        2.0 * psi[1:-1, :]
    ) / dy**2

    #
    # Boundaries.
    #
    # These currently use the same zero-normal-current
    # philosophy as the existing Laplacian infrastructure.
    #
    # Explicit superconducting boundary conditions will
    # eventually be handled by the TDGL boundary module.
    #

    # Left
    lap[:, 0] += (
        Ux[:, 0] * psi[:, 1]
        -
        psi[:, 0]
    ) / dx**2

    # Right
    lap[:, -1] += (
        np.conjugate(
            Ux[:, -2]
        ) * psi[:, -2]
        -
        psi[:, -1]
    ) / dx**2

    # Bottom
    lap[0, :] += (
        Uy[0, :] * psi[1, :]
        -
        psi[0, :]
    ) / dy**2

    # Top
    lap[-1, :] += (
        np.conjugate(
            Uy[-2, :]
        ) * psi[-2, :]
        -
        psi[-1, :]
    ) / dy**2

    return lap

def gauge_link_x(
    vector_potential_x,
    dx,
):
    """
    Construct the gauge link in the x direction.

    U_x = exp(-i A_x dx)

    The link represents the gauge phase accumulated
    across one grid spacing.
    """

    return np.exp(
        -1j *
        vector_potential_x *
        dx
    )


def gauge_link_y(
    vector_potential_y,
    dy,
):
    """
    Construct the gauge link in the y direction.

    U_y = exp(-i A_y dy)
    """

    return np.exp(
        -1j *
        vector_potential_y *
        dy
    )

def gauge_covariant_gradient(
    psi,
    vector_potential_x,
    vector_potential_y,
    dx,
    dy,
    active_mask=None,
):
    """
    Compute a gauge-covariant finite-difference gradient.

    Neighboring order parameters are compared through
    gauge links rather than by directly subtracting
    the vector potential.

    Forward-difference form:

        D_x psi =
            (U_x psi_{i+1} - psi_i) / dx

        D_y psi =
            (U_y psi_{j+1} - psi_j) / dy
    """

    Dx = np.zeros_like(
        psi,
        dtype=complex,
    )

    Dy = np.zeros_like(
        psi,
        dtype=complex,
    )

    Ux = gauge_link_x(
        vector_potential_x,
        dx,
    )

    Uy = gauge_link_y(
        vector_potential_y,
        dy,
    )

    #
    # x direction
    #

    if psi.shape[1] > 1:

        Dx[:, :-1] = (
            Ux[:, :-1] *
            psi[:, 1:]
            -
            psi[:, :-1]
        ) / dx

    #
    # y direction
    #

    if psi.shape[0] > 1:

        Dy[:-1, :] = (
            Uy[:-1, :] *
            psi[1:, :]
            -
            psi[:-1, :]
        ) / dy

    if active_mask is not None:
        active = np.asarray(active_mask, dtype=bool)
        if active.shape != psi.shape:
            raise ValueError("Active-domain mask must match psi.")
        Dx[:, :-1] *= active[:, :-1] & active[:, 1:]
        Dy[:-1, :] *= active[:-1, :] & active[1:, :]
        Dx[~active] = 0.0
        Dy[~active] = 0.0

    return Dx, Dy
