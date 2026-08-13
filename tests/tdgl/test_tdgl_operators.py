import numpy as np

from shs.tdgl import (
    gradient, 
    laplacian,
    covariant_gradient,
    covariant_laplacian,
)

def test_tdgl_gradient_constant():

    psi = np.ones(
        (10,10),
        dtype=complex
    )


    dx, dy = gradient(
        psi,
        1.0,
        1.0
    )


    assert np.allclose(
        dx,
        0.0
    )


    assert np.allclose(
        dy,
        0.0
    )



def test_tdgl_gradient_linear():

    x = np.arange(10)

    psi = np.tile(
        x,
        (10,1)
    ).astype(complex)


    dx, dy = gradient(
        psi,
        1.0,
        1.0
    )


    assert np.isclose(
        dx[5,5].real,
        1.0
    )


    assert np.isclose(
        dy[5,5].real,
        0.0
    )

import numpy as np

from shs.tdgl import covariant_gradient



def test_covariant_gradient_zero_field():

    psi = np.ones(
        (10,10),
        dtype=complex
    )


    Ax = np.zeros(
        (10,10)
    )

    Ay = np.zeros(
        (10,10)
    )


    dx, dy = covariant_gradient(
        psi,
        Ax,
        Ay,
        1.0,
        1.0,
    )


    assert np.allclose(
        dx,
        0.0
    )


    assert np.allclose(
        dy,
        0.0
    )



def test_covariant_gradient_phase_field():

    psi = np.ones(
        (10,10),
        dtype=complex
    )


    Ax = np.ones(
        (10,10)
    )

    Ay = np.zeros(
        (10,10)
    )


    dx, dy = covariant_gradient(
        psi,
        Ax,
        Ay,
        1.0,
        1.0,
    )


    # For constant psi:
    #
    # ∂psi/∂x = 0
    #
    # Dx psi = -iA psi


    assert np.isclose(
        dx[5,5].imag,
        -1.0
    )


    assert np.isclose(
        dy[5,5].imag,
        0.0
    )

def test_tdgl_laplacian_preserves_shape():

    psi = np.ones(
        (100, 100),
        dtype=complex
    )

    result = laplacian(
        psi,
        1.0,
        1.0,
    )

    assert result.shape == psi.shape

def test_covariant_laplacian_uniform_zero_field():

    psi = np.ones(
        (20, 20),
        dtype=complex,
    )

    Ax = np.zeros(
        (20, 20),
        dtype=float,
    )

    Ay = np.zeros(
        (20, 20),
        dtype=float,
    )

    result = covariant_laplacian(
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

def test_covariant_laplacian_uniform_vector_potential():

    psi = np.ones(
        (20, 20),
        dtype=complex,
    )

    Ax = np.ones(
        (20, 20),
        dtype=float,
    )

    Ay = np.zeros(
        (20, 20),
        dtype=float,
    )

    dx = 1.0
    dy = 1.0

    result = covariant_laplacian(
        psi,
        Ax,
        Ay,
        dx,
        dy,
    )

    expected = (
        2.0 *
        np.cos(1.0)
        -
        2.0
    )

    assert np.isclose(
        result[10, 10].real,
        expected,
    )

    assert np.isclose(
        result[10, 10].imag,
        0.0,
    )

def test_covariant_laplacian_gauge_transformation():

    nx = 20
    ny = 20

    dx = 0.2
    dy = 0.2

    x = np.arange(nx) * dx

    k = 0.3

    psi = np.ones(
        (ny, nx),
        dtype=complex,
    )

    Ax = np.zeros(
        (ny, nx),
        dtype=float,
    )

    Ay = np.zeros(
        (ny, nx),
        dtype=float,
    )

    original = covariant_laplacian(
        psi,
        Ax,
        Ay,
        dx,
        dy,
    )

    #
    # Gauge transformation:
    #
    # psi' = exp(i chi) psi
    #
    # A' = A + grad(chi)
    #

    chi = k * x

    psi_transformed = np.tile(
        np.exp(1j * chi),
        (ny, 1),
    )

    Ax_transformed = np.full(
        (ny, nx),
        k,
        dtype=float,
    )

    Ay_transformed = np.zeros(
        (ny, nx),
        dtype=float,
    )

    transformed = covariant_laplacian(
        psi_transformed,
        Ax_transformed,
        Ay_transformed,
        dx,
        dy,
    )

    expected = (
        np.exp(1j * chi)[None, :] *
        original
    )

    assert np.allclose(
        transformed[1:-1, 1:-1],
        expected[1:-1, 1:-1],
        atol=1e-10,
    )