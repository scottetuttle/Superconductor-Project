"""
Ginzburg-Landau mathematical helpers.

These functions provide small, reusable mathematical
operations for the normalized GL formulation.

The TDGLModel remains the authoritative definition
of the physical model.
"""

import numpy as np


def reduced_temperature(
    temperature,
    Tc,
):
    """
    Calculate reduced temperature:

        t = T / Tc
    """

    if np.any(np.asarray(Tc) <= 0.0):
        raise ValueError(
            "Critical temperature must be positive."
        )

    return (
        np.asarray(temperature) /
        np.asarray(Tc)
    )


def alpha_temperature(
    temperature,
    Tc,
    alpha0=1.0,
):
    """
    Calculate the normalized temperature-dependent
    quadratic GL coefficient.

    In the normalized formulation:

        alpha(T) = alpha0 * (1 - T/Tc)
    """

    return (
        alpha0 *
        (
            1.0 -
            reduced_temperature(
                temperature,
                Tc,
            )
        )
    )


def nonlinear_term(
    psi,
    beta=1.0,
):
    """
    Calculate the cubic GL nonlinear term:

        beta |psi|^2 psi
    """

    return (
        beta *
        np.abs(psi)**2 *
        psi
    )