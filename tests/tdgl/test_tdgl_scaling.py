import numpy as np
import pytest

from shs.tdgl.scaling import TDGLScales


def create_scales():
    """
    Create representative NbN TDGL scales for testing.
    """

    return TDGLScales(
        xi=5.0e-9,
        lambda_=2.0e-7,
        critical_temperature=15.5,
        time_scale=1.0e-12,
    )


def test_scaling_requires_positive_xi():

    with pytest.raises(ValueError):

        TDGLScales(
            xi=0.0,
            lambda_=2.0e-7,
            critical_temperature=15.5,
            time_scale=1.0e-12,
        )


def test_scaling_requires_positive_lambda():

    with pytest.raises(ValueError):

        TDGLScales(
            xi=5.0e-9,
            lambda_=0.0,
            critical_temperature=15.5,
            time_scale=1.0e-12,
        )


def test_scaling_requires_positive_tc():

    with pytest.raises(ValueError):

        TDGLScales(
            xi=5.0e-9,
            lambda_=2.0e-7,
            critical_temperature=0.0,
            time_scale=1.0e-12,
        )


def test_scaling_requires_positive_time():

    with pytest.raises(ValueError):

        TDGLScales(
            xi=5.0e-9,
            lambda_=2.0e-7,
            critical_temperature=15.5,
            time_scale=0.0,
        )


def test_length_conversion():

    scales = create_scales()

    physical_length = 5.0e-9

    dimensionless = (
        scales.length_to_dimensionless(
            physical_length
        )
    )

    assert np.isclose(
        dimensionless,
        1.0,
    )


def test_length_conversion_is_reversible():

    scales = create_scales()

    lengths = np.array([
        1.0e-9,
        5.0e-9,
        1.0e-7,
        5.0e-6,
    ])

    dimensionless = (
        scales.length_to_dimensionless(
            lengths
        )
    )

    recovered = (
        scales.length_to_physical(
            dimensionless
        )
    )

    assert np.allclose(
        recovered,
        lengths,
    )


def test_temperature_conversion():

    scales = create_scales()

    temperature = 15.5

    reduced = (
        scales.temperature_to_dimensionless(
            temperature
        )
    )

    assert np.isclose(
        reduced,
        1.0,
    )


def test_temperature_conversion_is_reversible():

    scales = create_scales()

    temperatures = np.array([
        0.0,
        3.0,
        7.75,
        15.5,
        20.0,
    ])

    reduced = (
        scales.temperature_to_dimensionless(
            temperatures
        )
    )

    recovered = (
        scales.temperature_to_physical(
            reduced
        )
    )

    assert np.allclose(
        recovered,
        temperatures,
    )


def test_time_conversion_is_reversible():

    scales = create_scales()

    times = np.array([
        1.0e-15,
        1.0e-12,
        1.0e-9,
        1.0e-6,
    ])

    dimensionless = (
        scales.time_to_dimensionless(
            times
        )
    )

    recovered = (
        scales.time_to_physical(
            dimensionless
        )
    )

    assert np.allclose(
        recovered,
        times,
    )


def test_vector_potential_scale():

    scales = create_scales()

    expected = (
        scales.flux_quantum /
        (
            2.0 *
            np.pi *
            scales.xi
        )
    )

    assert np.isclose(
        scales.vector_potential_scale,
        expected,
    )


def test_magnetic_field_scale():

    scales = create_scales()

    expected = (
        scales.flux_quantum /
        (
            2.0 *
            np.pi *
            scales.xi**2
        )
    )

    assert np.isclose(
        scales.magnetic_field_scale,
        expected,
    )


def test_current_density_scale():

    scales = create_scales()

    mu0 = 4.0e-7 * np.pi

    expected = (
        scales.flux_quantum /
        (
            2.0 *
            np.pi *
            mu0 *
            scales.lambda_**2 *
            scales.xi
        )
    )

    assert np.isclose(
        scales.current_density_scale,
        expected,
    )


def test_vector_potential_conversion_is_reversible():

    scales = create_scales()

    values = np.array([
        0.0,
        1.0e-9,
        1.0e-6,
        1.0e-3,
    ])

    dimensionless = (
        scales.vector_potential_to_dimensionless(
            values
        )
    )

    recovered = (
        scales.vector_potential_to_physical(
            dimensionless
        )
    )

    assert np.allclose(
        recovered,
        values,
    )


def test_magnetic_field_conversion_is_reversible():

    scales = create_scales()

    values = np.array([
        0.0,
        1.0e-6,
        1.0e-3,
        1.0,
    ])

    dimensionless = (
        scales.magnetic_field_to_dimensionless(
            values
        )
    )

    recovered = (
        scales.magnetic_field_to_physical(
            dimensionless
        )
    )

    assert np.allclose(
        recovered,
        values,
    )


def test_current_density_conversion_is_reversible():

    scales = create_scales()

    values = np.array([
        0.0,
        1.0e6,
        1.0e8,
        1.0e10,
    ])

    dimensionless = (
        scales.current_density_to_dimensionless(
            values
        )
    )

    recovered = (
        scales.current_density_to_physical(
            dimensionless
        )
    )

    assert np.allclose(
        recovered,
        values,
    )


def test_array_conversion_preserves_shape():

    scales = create_scales()

    values = np.ones(
        (20, 20)
    ) * 5.0e-9

    result = (
        scales.length_to_dimensionless(
            values
        )
    )

    assert result.shape == values.shape