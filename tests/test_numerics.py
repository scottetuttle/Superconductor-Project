import numpy as np

from shs.numerics import (
    gradient,
    laplacian,
    gauss_seidel,
)



def test_gradient():

    field = np.zeros((5,5))

    for y in range(5):
        for x in range(5):
            field[y,x] = x


    gx, gy = gradient(
        field,
        dx=1,
        dy=1,
    )


    assert np.isclose(
        gx[2,2],
        1.0
    )


    assert np.isclose(
        gy[2,2],
        0.0
    )



def test_laplacian():

    field = np.ones((5,5))


    result = laplacian(
        field,
        1,
        1,
    )


    assert np.allclose(
        result,
        0
    )



def test_gauss_seidel():

    shape = (10,10)


    V = np.zeros(shape)

    sigma = np.ones(shape)

    source = np.zeros(shape)


    boundaries = np.zeros(
        shape,
        dtype=bool
    )


    values = np.zeros(shape)


    boundaries[:,0] = True
    values[:,0] = 1


    result = gauss_seidel(
        V,
        sigma,
        source,
        boundaries,
        values,
        1,
        1,
    )


    assert result.field[:,0].mean() == 1

    assert result.converged

    assert result.iterations > 0