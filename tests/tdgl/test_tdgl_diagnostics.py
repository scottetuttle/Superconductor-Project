import numpy as np

from shs.tdgl.diagnostics import (
    order_parameter_amplitude,
    order_parameter_amplitude_squared,
    order_parameter_phase,
    supercurrent_magnitude,
    covariant_gradient_magnitude,
    free_energy_density,
    total_free_energy,
    mean_order_parameter_amplitude,
    maximum_order_parameter_amplitude,
    minimum_order_parameter_amplitude,
)


def test_order_parameter_amplitude():

    psi = np.array(
        [
            3.0 + 4.0j,
            1.0 + 0.0j,
        ]
    )

    result = order_parameter_amplitude(psi)

    assert np.allclose(
        result,
        [5.0, 1.0],
    )


def test_order_parameter_amplitude_squared():

    psi = np.array(
        [
            3.0 + 4.0j,
            1.0 + 0.0j,
        ]
    )

    result = order_parameter_amplitude_squared(psi)

    assert np.allclose(
        result,
        [25.0, 1.0],
    )


def test_order_parameter_phase():

    psi = np.array(
        [
            1.0 + 0.0j,
            0.0 + 1.0j,
        ]
    )

    result = order_parameter_phase(psi)

    assert np.isclose(
        result[0],
        0.0,
    )

    assert np.isclose(
        result[1],
        np.pi / 2,
    )


def test_supercurrent_magnitude():

    jx = np.array(
        [
            [3.0, 0.0],
            [0.0, 1.0],
        ]
    )

    jy = np.array(
        [
            [4.0, 0.0],
            [0.0, 0.0],
        ]
    )

    result = supercurrent_magnitude(
        jx,
        jy,
    )

    assert np.allclose(
        result,
        [
            [5.0, 0.0],
            [0.0, 1.0],
        ],
    )


def test_covariant_gradient_zero_for_uniform_zero_field():

    psi = np.ones(
        (10, 10),
        dtype=complex,
    )

    Ax = np.zeros(
        (10, 10)
    )

    Ay = np.zeros(
        (10, 10)
    )

    result = covariant_gradient_magnitude(
        psi,
        Ax,
        Ay,
        1.0,
        1.0,
    )

    assert np.allclose(
        result,
        0.0,
    )


def test_covariant_gradient_detects_vector_potential():

    psi = np.ones(
        (10, 10),
        dtype=complex,
    )

    Ax = np.ones(
        (10, 10)
    )

    Ay = np.zeros(
        (10, 10)
    )

    result = covariant_gradient_magnitude(
        psi,
        Ax,
        Ay,
        1.0,
        1.0,
    )

    assert np.allclose(
        result,
        1.0,
    )


def test_free_energy_density_equilibrium():

    reduced_temperature = 0.5

    psi = np.full(
        (10, 10),
        np.sqrt(0.5),
        dtype=complex,
    )

    result = free_energy_density(
        psi,
        reduced_temperature,
    )

    expected = -0.5 * (0.5)**2

    assert np.allclose(
        result,
        expected,
    )


def test_free_energy_density_normal_state():

    psi = np.zeros(
        (10, 10),
        dtype=complex,
    )

    result = free_energy_density(
        psi,
        1.5,
    )

    assert np.allclose(
        result,
        0.0,
    )


def test_total_free_energy():

    psi = np.ones(
        (4, 4),
        dtype=complex,
    )

    result = total_free_energy(
        psi,
        0.0,
    )

    expected = (
        4 * 4 *
        (
            -1.0
            +
            0.5
        )
    )

    assert np.isclose(
        result,
        expected,
    )


def test_amplitude_statistics():

    psi = np.array(
        [
            [1.0 + 0.0j, 2.0 + 0.0j],
            [3.0 + 0.0j, 4.0 + 0.0j],
        ]
    )

    assert np.isclose(
        mean_order_parameter_amplitude(psi),
        2.5,
    )

    assert np.isclose(
        minimum_order_parameter_amplitude(psi),
        1.0,
    )

    assert np.isclose(
        maximum_order_parameter_amplitude(psi),
        4.0,
    )