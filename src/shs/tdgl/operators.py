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
):
    """
    Compute the gauge-covariant Laplacian.

    For

        D = ∇ - iA

    and assuming

        ∇ · A = 0,

    the operator becomes:

        D²ψ =
            ∇²ψ
            - 2i A·∇ψ
            - |A|² ψ

    The returned array has the same shape as psi.

    Notes
    -----
    The current implementation assumes a
    divergence-free vector potential.

    A fully gauge-consistent implementation,
    including arbitrary gauge choices and explicit
    ∇·A contributions, will be developed as the
    electromagnetic coupling is expanded.
    """

    lap = laplacian(
        psi,
        dx,
        dy,
    )

    dpsi_dx, dpsi_dy = gradient(
        psi,
        dx,
        dy,
    )

    Ax = vector_potential_x
    Ay = vector_potential_y

    result = (

        lap

        -
        2j * (
            Ax * dpsi_dx
            +
            Ay * dpsi_dy
        )

        -
        (
            Ax**2
            +
            Ay**2
        ) * psi
    )

    return result