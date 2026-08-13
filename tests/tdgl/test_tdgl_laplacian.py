import numpy as np

from shs.tdgl import (
    covariant_laplacian, 
    laplacian
)


def test_covariant_laplacian_constant_state():

    psi = np.ones(
        (20,20),
        dtype=complex
    )


    Ax = np.zeros(
        (20,20)
    )

    Ay = np.zeros(
        (20,20)
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
        0.0
    )



def test_covariant_laplacian_uniform_vector_potential():

    psi = np.ones(
        (20,20),
        dtype=complex
    )


    Ax = np.ones(
        (20,20)
    )

    Ay = np.zeros(
        (20,20)
    )


    result = covariant_laplacian(
        psi,
        Ax,
        Ay,
        1.0,
        1.0,
    )


    # Constant psi:
    #
    # D²psi = -A²psi

    expected = (
        2.0 * np.cos(1.0)
        - 2.0
    )

    assert np.isclose(
        result[10, 10].real,
        expected
    )

def test_laplacian_constant_state_including_boundaries():

    psi = np.ones(
        (20, 20),
        dtype=complex
    )

    result = laplacian(
        psi,
        1.0,
        1.0,
    )

    assert np.allclose(
        result,
        0.0
    )