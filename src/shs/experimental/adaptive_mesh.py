"""Conservative prototypes for nested meshes and temporal warm starts.

Nothing in this module is imported by the production solver path.
"""

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import RegularGridInterpolator


@dataclass(frozen=True)
class HistoryPrediction:
    value: np.ndarray
    weight: np.ndarray
    normalized_curvature: np.ndarray


@dataclass(frozen=True)
class RestrictedGaugeFields:
    psi: np.ndarray
    ax_dimensionless: np.ndarray
    ay_dimensionless: np.ndarray


def history_predict(previous2, previous, current, *, maximum_weight=0.8,
                    scale_floor=1e-12):
    """Damped linear extrapolation with local confidence from second differences.

    Smooth, consistent changes approach ``maximum_weight``. Abrupt changes reduce
    the extrapolation toward the last accepted value. A stationary location has
    high confidence but a zero extrapolation increment.
    """
    a, b, c = (np.asarray(value) for value in (previous2, previous, current))
    if a.shape != b.shape or b.shape != c.shape:
        raise ValueError("History arrays must have identical shapes.")
    if not 0 <= maximum_weight <= 1 or scale_floor <= 0:
        raise ValueError("maximum_weight must lie in [0,1] and scale_floor be positive.")
    d0, d1 = b - a, c - b
    scale = np.abs(d0) + np.abs(d1) + scale_floor
    curvature = np.abs(d1 - d0) / scale
    weight = maximum_weight * np.exp(-4.0 * curvature)
    return HistoryPrediction(c + weight * d1, weight, curvature)


def history_predict_order_parameter(previous2, previous, current, *, maximum_weight=0.8,
                                    amplitude_floor=1e-8):
    """Predict complex psi using amplitude and wrapped temporal phase increments.

    This avoids interpolating a wrapped phase directly. It is still only
    temporally gauge-consistent when the scalar-potential gauge is unchanged.
    Production use requires a gauge-covariant acceptance check.
    """
    a, b, c = (np.asarray(value, dtype=complex) for value in
               (previous2, previous, current))
    amplitude_prediction = history_predict(
        np.abs(a), np.abs(b), np.abs(c), maximum_weight=maximum_weight,
        scale_floor=amplitude_floor)
    phase0 = np.angle(b * np.conjugate(a))
    phase1 = np.angle(c * np.conjugate(b))
    phase_curvature = np.abs(np.angle(np.exp(1j * (phase1 - phase0)))) / np.pi
    phase_weight = maximum_weight * np.exp(-4.0 * phase_curvature)
    weight = np.minimum(amplitude_prediction.weight, phase_weight)
    amplitude = np.maximum(0.0, np.abs(c) + weight * (np.abs(c) - np.abs(b)))
    phase = np.angle(c) + weight * phase1
    inactive = (np.abs(a) < amplitude_floor) | (np.abs(b) < amplitude_floor) | (np.abs(c) < amplitude_floor)
    weight = np.where(inactive, 0.0, weight)
    value = np.where(inactive, c, amplitude * np.exp(1j * phase))
    return HistoryPrediction(value, weight,
                             np.maximum(amplitude_prediction.normalized_curvature,
                                        phase_curvature))


def _coarse_reconstruction(field, factor):
    ny, nx = field.shape
    if factor < 2 or (nx - 1) % factor or (ny - 1) % factor:
        raise ValueError("factor must divide both node-interval counts.")
    y = np.arange(ny, dtype=float)
    x = np.arange(nx, dtype=float)
    yc, xc = y[::factor], x[::factor]

    def interpolate(values):
        interpolator = RegularGridInterpolator((yc, xc), values[::factor, ::factor],
                                               method="linear", bounds_error=True)
        yy, xx = np.meshgrid(y, x, indexing="ij")
        return interpolator(np.column_stack((yy.ravel(), xx.ravel()))).reshape(field.shape)

    if np.iscomplexobj(field):
        return interpolate(np.real(field)) + 1j * interpolate(np.imag(field))
    return interpolate(np.asarray(field, dtype=float))


def gauge_covariant_reconstruction(psi, ax_dimensionless, ay_dimensionless, *,
                                   factor, dx_dimensionless, dy_dimensionless):
    """Prolong injected coarse psi with fine-link parallel transport.

    ``ax_dimensionless`` and ``ay_dimensionless`` define the target-grid gauge
    connection. Coarse nodes are the every-``factor`` fine nodes. Four coarse
    corners are transported along deterministic x-then-y paths to each target
    node before bilinear interpolation, making the result transform covariantly
    under a gauge change of the target connection.
    """
    psi = np.asarray(psi, dtype=complex)
    ax = np.asarray(ax_dimensionless, dtype=float)
    ay = np.asarray(ay_dimensionless, dtype=float)
    if psi.ndim != 2 or ax.shape != psi.shape or ay.shape != psi.shape:
        raise ValueError("psi and vector potential must share a two-dimensional shape.")
    ny, nx = psi.shape
    factor = int(factor)
    if factor < 2 or (nx - 1) % factor or (ny - 1) % factor:
        raise ValueError("factor must divide both node-interval counts.")
    if min(dx_dimensionless, dy_dimensionless) <= 0:
        raise ValueError("Dimensionless spacings must be positive.")
    ux = np.exp(-1j * ax[:, :-1] * dx_dimensionless)
    uy = np.exp(-1j * ay[:-1, :] * dy_dimensionless)

    def transport(value, ys, xs, yt, xt):
        # Move in x at the source row, then in y at the target column.
        if xt > xs:
            value *= np.prod(np.conjugate(ux[ys, xs:xt]))
        elif xt < xs:
            value *= np.prod(ux[ys, xt:xs])
        if yt > ys:
            value *= np.prod(np.conjugate(uy[ys:yt, xt]))
        elif yt < ys:
            value *= np.prod(uy[yt:ys, xt])
        return value

    result = np.empty_like(psi)
    for y in range(ny):
        y0 = (y // factor) * factor
        y1 = min(y0 + factor, ny - 1)
        fy = 0.0 if y1 == y0 else (y - y0) / (y1 - y0)
        for x in range(nx):
            x0 = (x // factor) * factor
            x1 = min(x0 + factor, nx - 1)
            fx = 0.0 if x1 == x0 else (x - x0) / (x1 - x0)
            corners = ((y0, x0, (1-fx)*(1-fy)),
                       (y0, x1, fx*(1-fy)),
                       (y1, x0, (1-fx)*fy),
                       (y1, x1, fx*fy))
            result[y, x] = sum(weight * transport(psi[ys, xs], ys, xs, y, x)
                               for ys, xs, weight in corners if weight)
    return result


def gauge_covariant_restriction(psi, ax_dimensionless, ay_dimensionless, *,
                                factor, dx_dimensionless, dy_dimensionless):
    """Inject coarse psi and compose fine link variables into coarse links."""
    psi = np.asarray(psi, dtype=complex)
    ax = np.asarray(ax_dimensionless, dtype=float)
    ay = np.asarray(ay_dimensionless, dtype=float)
    if psi.ndim != 2 or ax.shape != psi.shape or ay.shape != psi.shape:
        raise ValueError("psi and vector potential must share a two-dimensional shape.")
    ny, nx = psi.shape
    factor = int(factor)
    if factor < 2 or (nx-1) % factor or (ny-1) % factor:
        raise ValueError("factor must divide both node-interval counts.")
    if min(dx_dimensionless, dy_dimensionless) <= 0:
        raise ValueError("Dimensionless spacings must be positive.")
    coarse_psi = psi[::factor, ::factor].copy()
    cy, cx = coarse_psi.shape
    coarse_ax = np.zeros((cy, cx))
    coarse_ay = np.zeros((cy, cx))
    for j, y in enumerate(range(0, ny, factor)):
        for i, x in enumerate(range(0, nx-factor, factor)):
            holonomy = np.prod(np.exp(-1j * ax[y, x:x+factor] * dx_dimensionless))
            coarse_ax[j, i] = -np.angle(holonomy) / (factor * dx_dimensionless)
        coarse_ax[j, -1] = coarse_ax[j, -2]
    for j, y in enumerate(range(0, ny-factor, factor)):
        for i, x in enumerate(range(0, nx, factor)):
            holonomy = np.prod(np.exp(-1j * ay[y:y+factor, x] * dy_dimensionless))
            coarse_ay[j, i] = -np.angle(holonomy) / (factor * dy_dimensionless)
    coarse_ay[-1] = coarse_ay[-2]
    return RestrictedGaugeFields(coarse_psi, coarse_ax, coarse_ay)


def composite_coarse_exterior(field, *, factor, fine_center_cell, fine_radius_cells,
                              transition_cells=2):
    """Reconstruct a coarse exterior on the fine grid and retain a fine patch.

    The returned uniform array is a validation representation, not an efficient
    AMR data structure. A cosine transition avoids a hard coarse/fine seam.
    """
    original = np.asarray(field)
    if original.ndim != 2:
        raise ValueError("field must be two-dimensional.")
    center = np.asarray(fine_center_cell, dtype=float)
    if center.shape != (2,) or fine_radius_cells <= 0 or transition_cells < 0:
        raise ValueError("Invalid fine-patch geometry.")
    coarse = _coarse_reconstruction(original, int(factor))
    yy, xx = np.meshgrid(np.arange(original.shape[0]), np.arange(original.shape[1]),
                         indexing="ij")
    radius = np.hypot(xx - center[0], yy - center[1])
    if transition_cells == 0:
        fine_weight = (radius <= fine_radius_cells).astype(float)
    else:
        u = np.clip((radius - fine_radius_cells) / transition_cells, 0.0, 1.0)
        fine_weight = 0.5 * (1.0 + np.cos(np.pi * u))
        fine_weight[radius <= fine_radius_cells] = 1.0
        fine_weight[radius >= fine_radius_cells + transition_cells] = 0.0
    composite = fine_weight * original + (1.0 - fine_weight) * coarse
    return composite, fine_weight, coarse
