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