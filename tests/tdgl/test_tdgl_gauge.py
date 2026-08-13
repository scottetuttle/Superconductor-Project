import numpy as np

from shs.tdgl.operators import (
    gauge_link_x,
    gauge_link_y,
    gauge_covariant_gradient,
)


def test_gauge_links_have_unit_magnitude():

    Ax = np.random.random((20, 20))
    Ay = np.random.random((20, 20))

    Ux = gauge_link_x(
        Ax,
        0.5,
    )

    Uy = gauge_link_y(
        Ay,
        0.5,
    )

    assert np.allclose(
        np.abs(Ux),
        1.0,
    )

    assert np.allclose(
        np.abs(Uy),
        1.0,
    )


def test_zero_vector_potential_reduces_to_gradient():

    x = np.arange(20)
    y = np.arange(20)

    X, Y = np.meshgrid(
        x,
        y,
    )

    psi = (
        X.astype(float)
        +
        1j * Y.astype(float)
    )

    Ax = np.zeros_like(X, dtype=float)
    Ay = np.zeros_like(Y, dtype=float)

    dx = 1.0
    dy = 1.0

    Dx, Dy = gauge_covariant_gradient(
        psi,
        Ax,
        Ay,
        dx,
        dy,
    )

    assert np.allclose(
        Dx[:, :-1],
        1.0,
    )

    assert np.allclose(
        Dy[:-1, :],
        1j,
    )
    
def test_gauge_transformation_preserves_covariant_gradient():

    nx = 30
    ny = 30

    x = np.arange(nx)
    y = np.arange(ny)

    X, Y = np.meshgrid(
        x,
        y,
    )

    dx = 0.5
    dy = 0.5

    #
    # Original fields
    #

    psi = np.exp(
        1j * (
            0.2 * X
            +
            0.15 * Y
        )
    )

    Ax = np.full(
        (ny, nx),
        0.3,
    )

    Ay = np.full(
        (ny, nx),
        -0.2,
    )

    #
    # Gauge transformation
    #
    # chi(x,y)
    #

    chi = (
        0.1 * X**2
        +
        0.05 * Y**2
    )

    psi_transformed = (
        psi *
        np.exp(1j * chi)
    )

    #
    # Approximate gauge transformation of A
    #

    dchi_dx = np.zeros_like(chi)

    dchi_dy = np.zeros_like(chi)

    dchi_dx[:, :-1] = (
        chi[:, 1:] -
        chi[:, :-1]
    ) / dx

    dchi_dy[:-1, :] = (
        chi[1:, :] -
        chi[:-1, :]
    ) / dy

    Ax_transformed = (
        Ax +
        dchi_dx
    )

    Ay_transformed = (
        Ay +
        dchi_dy
    )

    #
    # Covariant derivatives
    #

    Dx, Dy = gauge_covariant_gradient(
        psi,
        Ax,
        Ay,
        dx,
        dy,
    )

    Dx_transformed, Dy_transformed = (
        gauge_covariant_gradient(
            psi_transformed,
            Ax_transformed,
            Ay_transformed,
            dx,
            dy,
        )
    )

    #
    # Gauge covariance:
    #
    # D'psi' = exp(i chi) D psi
    #

    expected_Dx = (
        np.exp(1j * chi) *
        Dx
    )

    expected_Dy = (
        np.exp(1j * chi) *
        Dy
    )

    assert np.allclose(
        Dx_transformed[:, :-1],
        expected_Dx[:, :-1],
        atol=1e-10,
    )

    assert np.allclose(
        Dy_transformed[:-1, :],
        expected_Dy[:-1, :],
        atol=1e-10,
    )