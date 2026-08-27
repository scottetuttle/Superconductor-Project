import numpy as np

from shs.tdgl.initialization import (
    equilibrium_superconducting_state,
)


def test_equilibrium_initialization_at_half_tc():

    psi = equilibrium_superconducting_state(
        shape=(10, 10),
        reduced_temperature=0.5,
    )

    assert psi.shape == (10, 10)

    assert np.allclose(
        np.abs(psi),
        np.sqrt(0.5),
    )


def test_equilibrium_initialization_above_tc():

    psi = equilibrium_superconducting_state(
        shape=(10, 10),
        reduced_temperature=1.1,
    )

    assert np.allclose(
        psi,
        0.0,
    )


def test_equilibrium_initialization_preserves_phase():

    psi = equilibrium_superconducting_state(
        shape=(10, 10),
        reduced_temperature=0.5,
        phase=np.pi / 2,
    )

    assert np.allclose(
        np.abs(psi),
        np.sqrt(0.5),
    )

    assert np.allclose(
        psi.real,
        0.0,
        atol=1e-12,
    )