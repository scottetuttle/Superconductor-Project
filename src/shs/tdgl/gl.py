import numpy as np


def reduced_temperature(
    temperature,
    Tc,
):
    return temperature / Tc



def alpha_temperature(
    temperature,
    Tc,
    alpha0,
):
    """
    Temperature dependent GL alpha.
    """

    return (
        alpha0 *
        (
            1 -
            temperature/Tc
        )
    )



def nonlinear_term(
    psi,
    beta,
):
    return (
        beta *
        np.abs(psi)**2 *
        psi
    )