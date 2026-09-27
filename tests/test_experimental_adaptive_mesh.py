import numpy as np

from shs.experimental.adaptive_mesh import (
    composite_coarse_exterior,
    gauge_covariant_reconstruction,
    gauge_covariant_restriction,
    history_predict,
    history_predict_order_parameter,
)


def test_history_predicts_smooth_trend_and_backs_off_on_reversal():
    smooth = history_predict(np.array([0.0]), np.array([1.0]), np.array([2.0]),
                             maximum_weight=.8)
    reversal = history_predict(np.array([0.0]), np.array([1.0]), np.array([0.0]),
                               maximum_weight=.8)
    assert smooth.value[0] > 2.7
    assert smooth.weight[0] > .7
    assert reversal.weight[0] < .02
    assert abs(reversal.value[0]) < .02


def test_order_parameter_prediction_respects_wrapped_phase_increment():
    phase = np.array([3.0, 3.1, -3.083185307179586])
    result = history_predict_order_parameter(
        *[np.array([np.exp(1j * value)]) for value in phase])
    predicted_increment = np.angle(result.value[0] * np.exp(-1j * phase[-1]))
    assert predicted_increment > 0
    assert predicted_increment < .1


def test_composite_keeps_fine_patch_exact_and_coarsens_exterior():
    y, x = np.mgrid[:9, :9]
    field = x**2 + y**2
    composite, weight, coarse = composite_coarse_exterior(
        field, factor=2, fine_center_cell=[4, 4], fine_radius_cells=1,
        transition_cells=1)
    assert np.allclose(composite[weight == 1], field[weight == 1])
    assert np.allclose(composite[weight == 0], coarse[weight == 0])
    assert np.any(np.abs(composite - field) > 0)


def test_gauge_covariant_reconstruction_commutes_with_gauge_change():
    y, x = np.mgrid[:5, :5]
    psi = (1 + .05*x) * np.exp(1j * (.2*x - .1*y))
    ax = np.full((5, 5), .15)
    ay = np.full((5, 5), -.08)
    chi = .03*x*x + .04*x*y - .02*y*y
    dx = dy = .25
    transformed_psi = psi * np.exp(1j * chi)
    transformed_ax = ax + np.diff(chi, axis=1, append=chi[:, -1:]) / dx
    transformed_ay = ay + np.diff(chi, axis=0, append=chi[-1:, :]) / dy
    first = gauge_covariant_reconstruction(
        psi, ax, ay, factor=2, dx_dimensionless=dx, dy_dimensionless=dy)
    second = gauge_covariant_reconstruction(
        transformed_psi, transformed_ax, transformed_ay, factor=2,
        dx_dimensionless=dx, dy_dimensionless=dy)
    assert np.allclose(second, first * np.exp(1j * chi), atol=1e-12)


def test_gauge_restriction_composes_link_holonomy():
    psi = np.ones((5, 5), dtype=complex)
    ax = np.arange(25, dtype=float).reshape(5, 5) * .01
    ay = -ax
    restricted = gauge_covariant_restriction(
        psi, ax, ay, factor=2, dx_dimensionless=.2, dy_dimensionless=.2)
    fine_product = np.prod(np.exp(-1j * ax[0, :2] * .2))
    coarse_link = np.exp(-1j * restricted.ax_dimensionless[0, 0] * .4)
    assert np.allclose(coarse_link, fine_product)
    assert restricted.psi.shape == (3, 3)
