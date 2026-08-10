"""
TDGL initialization utilities.

Creates initial superconducting order-parameter states.
"""

import numpy as np


def uniform_superconducting_state(
    shape,
    amplitude=1.0,
    phase=0.0,
):
    """
    Create a uniform superconducting state.

    Parameters
    ----------
    shape:
        Shape of the simulation mesh.

    amplitude:
        Order-parameter amplitude.

    phase:
        Uniform phase in radians.

    Returns
    -------
    ndarray
        Complex order parameter.
    """

    return np.full(
        shape,
        amplitude * np.exp(1j * phase),
        dtype=complex,
    )


def equilibrium_superconducting_state(
    shape,
    reduced_temperature,
    phase=0.0,
):
    """
    Create a uniform equilibrium superconducting state.

    The normalized equilibrium amplitude is:

        |psi| = sqrt(1 - T/Tc)

    for T < Tc.

    Above Tc the normal-state solution is:

        psi = 0.
    """

    amplitude = np.sqrt(
        np.maximum(
            1.0 -
            np.asarray(
                reduced_temperature,
                dtype=float
            ),
            0.0,
        )
    )

    if np.ndim(amplitude) != 0:
        raise ValueError(
            "equilibrium_superconducting_state expects "
            "a scalar reduced temperature."
        )

    return uniform_superconducting_state(
        shape,
        amplitude=float(amplitude),
        phase=phase,
    )