"""
TDGL initialization utilities.

Creates initial superconducting order parameter states.
"""

import numpy as np


def uniform_superconducting_state(shape):
    """
    Create a uniform superconducting state.

    Returns:

        psi = 1 + 0j

    Corresponds to:

        |psi| = 1
        phase = 0
    """

    return np.ones(
        shape,
        dtype=complex
    )