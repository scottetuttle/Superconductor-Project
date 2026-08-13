import numpy as np

from shs.tdgl import (
    equilibrium_superconducting_state,
    uniform_superconducting_state,
)


def test_uniform_superconducting_state():

    psi = uniform_superconducting_state(
        (10, 10)
    )

    assert psi.shape == (10, 10)

    assert psi.dtype == complex

    assert np.allclose(
        psi,
        1.0 + 0.0j,
    )


def test_uniform_superconducting_state_with_amplitude():

    psi = uniform_superconducting_state(
        (10, 10),
        amplitude=0.5,
    )

    assert np.allclose(
        np.abs(psi),
        0.5,
    )


def test_equilibrium_state_at_zero_temperature():

    psi = equilibrium_superconducting_state(
        (10, 10),
        reduced_temperature=0.0,
    )

    assert np.allclose(
        np.abs(psi),
        1.0,
    )


def test_equilibrium_state_at_half_tc():

    psi = equilibrium_superconducting_state(
        (10, 10),
        reduced_temperature=0.5,
    )

    assert np.allclose(
        np.abs(psi),
        np.sqrt(0.5),
    )


def test_equilibrium_state_near_tc():

    psi = equilibrium_superconducting_state(
        (10, 10),
        reduced_temperature=0.95,
    )

    assert np.allclose(
        np.abs(psi),
        np.sqrt(0.05),
    )


def test_equilibrium_state_above_tc():

    psi = equilibrium_superconducting_state(
        (10, 10),
        reduced_temperature=1.1,
    )

    assert np.allclose(
        psi,
        0.0 + 0.0j,
    )