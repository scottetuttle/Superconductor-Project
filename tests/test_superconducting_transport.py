import numpy as np

from shs.physics.superconducting_transport import (
    superconducting_current_density,
)


def test_superconducting_current_zero_when_psi_zero():

    psi = np.zeros(
        (10, 10),
        dtype=complex,
    )

    jx = np.ones(
        (10, 10)
    )

    jy = np.ones(
        (10, 10)
    )

    result_x, result_y = (
        superconducting_current_density(
            psi,
            jx,
            jy,
        )
    )

    assert np.allclose(
        result_x,
        0.0,
    )

    assert np.allclose(
        result_y,
        0.0,
    )


def test_superconducting_current_preserves_current_at_unit_psi():

    psi = np.ones(
        (10, 10),
        dtype=complex,
    )

    jx = np.full(
        (10, 10),
        2.0,
    )

    jy = np.full(
        (10, 10),
        -3.0,
    )

    result_x, result_y = (
        superconducting_current_density(
            psi,
            jx,
            jy,
        )
    )

    assert np.allclose(
        result_x,
        2.0,
    )

    assert np.allclose(
        result_y,
        -3.0,
    )


def test_superconducting_current_scales_with_order_parameter():

    psi = np.full(
        (10, 10),
        np.sqrt(0.5),
        dtype=complex,
    )

    jx = np.ones(
        (10, 10)
    )

    jy = np.ones(
        (10, 10)
    )

    result_x, result_y = (
        superconducting_current_density(
            psi,
            jx,
            jy,
        )
    )

    assert np.allclose(
        result_x,
        0.5,
    )

    assert np.allclose(
        result_y,
        0.5,
    )