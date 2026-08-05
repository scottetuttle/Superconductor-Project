"""
TDGL numerical operators.

Contains differential operators acting on
the superconducting order parameter.
"""

import numpy as np



def gradient(
    psi,
    dx,
    dy,
):
    """
    Compute spatial gradient of complex order parameter.

    Returns:

    dpsi/dx
    dpsi/dy
    """

    dpsi_dx = np.zeros_like(
        psi,
        dtype=complex
    )

    dpsi_dy = np.zeros_like(
        psi,
        dtype=complex
    )


    dpsi_dx[:,1:-1] = (
        psi[:,2:]
        -
        psi[:,:-2]
    ) / (
        2 * dx
    )


    dpsi_dy[1:-1,:] = (
        psi[2:,:]
        -
        psi[:-2,:]
    ) / (
        2 * dy
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
    Standard Laplacian.
    """

    lap = np.zeros_like(
        psi,
        dtype=complex
    )


    lap[1:-1,1:-1] = (

        psi[2:,1:-1]
        +
        psi[:-2,1:-1]

        +
        psi[1:-1,2:]
        +
        psi[1:-1,:-2]

        -
        4 *
        psi[1:-1,1:-1]

    ) / (
        dx * dy
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
    Gauge covariant gradient:

    D = ∇ - iA
    """

    dpsi_dx, dpsi_dy = gradient(
        psi,
        dx,
        dy,
    )


    Dx = (
        dpsi_dx
        -
        1j *
        vector_potential_x *
        psi
    )


    Dy = (
        dpsi_dy
        -
        1j *
        vector_potential_y *
        psi
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
    Gauge covariant Laplacian.

    Computes:

        (∇ - iA)^2 ψ


    Assumes:

        div(A)=0

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
        2j *
        (
            Ax *
            dpsi_dx
            +
            Ay *
            dpsi_dy
        )

        -
        (
            Ax**2
            +
            Ay**2
        )
        *
        psi

    )


    return result