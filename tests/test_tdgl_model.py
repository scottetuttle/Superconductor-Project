import numpy as np

from shs.tdgl import (
    TDGLParameters,
    TDGLModel,
)


def test_tdgl_parameters_are_valid():

    parameters = TDGLParameters()

    assert parameters.validate()


def test_tdgl_equilibrium_amplitude_zero_temperature():

    model = TDGLModel(
        TDGLParameters()
    )

    assert np.isclose(
        model.equilibrium_amplitude(0.0),
        1.0,
    )


def test_tdgl_equilibrium_amplitude_half_tc():

    model = TDGLModel(
        TDGLParameters()
    )

    expected = np.sqrt(0.5)

    assert np.isclose(
        model.equilibrium_amplitude(0.5),
        expected,
    )


def test_tdgl_equilibrium_amplitude_near_tc():

    model = TDGLModel(
        TDGLParameters()
    )

    expected = np.sqrt(0.05)

    assert np.isclose(
        model.equilibrium_amplitude(0.95),
        expected,
    )


def test_tdgl_equilibrium_amplitude_at_tc():

    model = TDGLModel(
        TDGLParameters()
    )

    assert np.isclose(
        model.equilibrium_amplitude(1.0),
        0.0,
    )


def test_tdgl_equilibrium_amplitude_above_tc():

    model = TDGLModel(
        TDGLParameters()
    )

    assert np.isclose(
        model.equilibrium_amplitude(1.2),
        0.0,
    )


def test_tdgl_equilibrium_amplitude_squared():

    model = TDGLModel(
        TDGLParameters()
    )

    assert np.isclose(
        model.equilibrium_amplitude_squared(0.5),
        0.5,
    )


def test_tdgl_equilibrium_order_parameter():

    model = TDGLModel(
        TDGLParameters()
    )

    psi = model.equilibrium_order_parameter(
        reduced_temperature=0.75,
        phase=0.0,
    )

    assert np.isclose(
        np.abs(psi),
        np.sqrt(0.25),
    )

    assert np.isclose(
        psi.imag,
        0.0,
    )


def test_tdgl_equilibrium_order_parameter_phase():

    model = TDGLModel(
        TDGLParameters()
    )

    psi = model.equilibrium_order_parameter(
        reduced_temperature=0.0,
        phase=np.pi / 2,
    )

    assert np.isclose(
        psi.real,
        0.0,
        atol=1e-12,
    )

    assert np.isclose(
        psi.imag,
        1.0,
        atol=1e-12,
    )