import numpy as np

from shs.tdgl import gradient


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