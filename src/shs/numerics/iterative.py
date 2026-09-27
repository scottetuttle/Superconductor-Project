"""Variable-coefficient Poisson solvers on a rectangular node grid.

Unspecified outer faces are insulating. Prescribed cells are Dirichlet.
Residual is the maximum diagonal-scaled equation defect (solution units),
not the change produced by a relaxation sweep.
"""
from dataclasses import dataclass
import numpy as np
from scipy.sparse import bmat, coo_matrix, csr_matrix
from scipy.sparse.linalg import factorized, spsolve
from shs.utils.defaults import default_section


_ELECTRICAL_SOLVER_DEFAULTS = default_section("electrical")["solver"]


@dataclass
class SolverResult:
    field: np.ndarray
    iterations: int
    residual: float
    converged: bool


def face_coefficients(coefficient, dx, dy):
    a = np.asarray(coefficient, dtype=float)
    x = np.zeros((a.shape[0], a.shape[1] - 1))
    y = np.zeros((a.shape[0] - 1, a.shape[1]))
    np.divide(2*a[:, :-1]*a[:, 1:], a[:, :-1]+a[:, 1:], out=x,
              where=(a[:, :-1]+a[:, 1:]) > 0)
    np.divide(2*a[:-1]*a[1:], a[:-1]+a[1:], out=y,
              where=(a[:-1]+a[1:]) > 0)
    return x / dx**2, y / dy**2


def _system(coefficient, dx, dy):
    x, y = face_coefficients(coefficient, dx, dy)
    diagonal = np.zeros_like(coefficient, dtype=float)
    diagonal[:, :-1] += x
    diagonal[:, 1:] += x
    diagonal[:-1] += y
    diagonal[1:] += y
    return x, y, diagonal


def _neighbors(v, x, y):
    total = np.zeros_like(v)
    total[:, :-1] += x*v[:, 1:]
    total[:, 1:] += x*v[:, :-1]
    total[:-1] += y*v[1:]
    total[1:] += y*v[:-1]
    return total


def sparse_direct(solution, coefficient, source, boundary_mask, boundary_values,
                  dx, dy, tolerance=_ELECTRICAL_SOLVER_DEFAULTS["tolerance"]):
    """Solve the same variable-coefficient Poisson system by sparse LU."""
    v = np.asarray(solution, dtype=float).copy()
    a, b = np.asarray(coefficient), np.asarray(source)
    fixed = np.asarray(boundary_mask, dtype=bool)
    values = np.asarray(boundary_values)
    if v.ndim != 2 or min(v.shape) < 2:
        raise ValueError('Poisson grid must be at least 2 by 2.')
    if any(arr.shape != v.shape for arr in (a, b, fixed, values)):
        raise ValueError('Poisson arrays must have identical shapes.')
    if not all(np.all(np.isfinite(arr)) for arr in (v, a, b, values)):
        raise ValueError('Poisson data must be finite.')
    if np.any(a < 0) or not np.isfinite(dx + dy) or min(dx, dy) <= 0:
        raise ValueError('Conductivity must be nonnegative and spacing positive.')
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError('Invalid Poisson convergence controls.')

    x, y, diagonal = _system(a, dx, dy)
    active = ~fixed & (diagonal > 0)
    if np.any((~fixed & (diagonal == 0)) & (b != 0)):
        raise ValueError('Nonzero source in an electrically isolated cell.')
    v[fixed] = values[fixed]
    count = int(np.count_nonzero(active))
    if count == 0:
        return SolverResult(v, 0, 0.0, True)

    indices = np.full(v.shape, -1, dtype=int)
    indices[active] = np.arange(count)
    rows = [np.arange(count)]
    cols = [np.arange(count)]
    data = [diagonal[active]]
    rhs = -b[active].copy()
    rhs_grid = np.zeros_like(v)
    rhs_grid[active] = rhs

    def add_faces(first, second, weights):
        first_active = active[first]
        second_active = active[second]
        both = first_active & second_active
        p = indices[first][both]
        q = indices[second][both]
        c = weights[both]
        rows.extend((p, q))
        cols.extend((q, p))
        data.extend((-c, -c))
        first_fixed = first_active & fixed[second]
        rhs_grid[first][first_fixed] += weights[first_fixed] * values[second][first_fixed]
        second_fixed = second_active & fixed[first]
        rhs_grid[second][second_fixed] += weights[second_fixed] * values[first][second_fixed]

    add_faces((slice(None), slice(None, -1)), (slice(None), slice(1, None)), x)
    add_faces((slice(None, -1), slice(None)), (slice(1, None), slice(None)), y)
    matrix = coo_matrix(
        (np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
        shape=(count, count),
    ).tocsr()
    v[active] = spsolve(matrix, rhs_grid[active])
    if not np.all(np.isfinite(v)):
        raise RuntimeError('Sparse Poisson solve produced nonfinite values.')
    defect = _neighbors(v, x, y) - diagonal*v - b
    residual = float(np.max(np.abs(defect[active]/diagonal[active]), initial=0))
    return SolverResult(v, 1, residual, residual <= tolerance)


def sparse_neumann(solution, coefficient, source, dx, dy, reference_value=0.0,
                   tolerance=_ELECTRICAL_SOLVER_DEFAULTS["tolerance"], thickness=None,
                   active_mask=None, factorization_cache=None):
    """Solve a finite-volume pure-Neumann problem with fixed mean voltage."""
    v = np.asarray(solution, dtype=float).copy()
    a, b = np.asarray(coefficient), np.asarray(source)
    if v.ndim != 2 or min(v.shape) < 2 or a.shape != v.shape or b.shape != v.shape:
        raise ValueError("Neumann Poisson arrays must be matching two-dimensional grids.")
    if not all(np.all(np.isfinite(array)) for array in (v, a, b)):
        raise ValueError("Neumann Poisson data must be finite.")
    active = np.ones_like(v, dtype=bool) if active_mask is None else np.asarray(active_mask, dtype=bool)
    if active.shape != v.shape or not np.any(active):
        raise ValueError("Neumann active mask must match the grid and contain active nodes.")
    if np.any(a[active] <= 0) or np.any(a[~active] < 0) or min(dx, dy) <= 0 or not np.isfinite(dx + dy):
        raise ValueError("Pure-Neumann conductivity must be positive in the active domain.")
    if not np.isfinite(reference_value) or tolerance <= 0:
        raise ValueError("Neumann reference value and tolerance must be finite and valid.")
    t = np.ones_like(a) if thickness is None else np.asarray(thickness, dtype=float)
    if t.shape != v.shape or np.any(t[active] <= 0) or not np.all(np.isfinite(t)):
        raise ValueError("Neumann thickness must be positive, finite, and match the grid.")
    wx, wy = np.ones(v.shape[1]), np.ones(v.shape[0])
    wx[[0, -1]], wy[[0, -1]] = 0.5, 0.5
    volumes = t * wy[:, None] * wx[None, :] * dx * dy * active
    integrated_source = b * volumes
    compatibility_scale = max(float(np.max(np.abs(integrated_source))), 1.0) * b.size
    if abs(float(np.sum(integrated_source))) > 256 * np.finfo(float).eps * compatibility_scale:
        raise ValueError("Pure-Neumann Poisson source does not satisfy compatibility.")

    cache = {} if factorization_cache is None else factorization_cache
    key = (a.shape, float(dx), float(dy), a.tobytes(), t.tobytes(), active.tobytes())
    cached = cache.get("neumann")
    sigma_x, sigma_y = face_coefficients(a, 1.0, 1.0)
    sigma_x *= active[:, :-1] & active[:, 1:]
    sigma_y *= active[:-1, :] & active[1:, :]
    thickness_x = 0.5 * (t[:, :-1] + t[:, 1:])
    thickness_y = 0.5 * (t[:-1] + t[1:])
    x = sigma_x * thickness_x * wy[:, None] * dy / dx
    y = sigma_y * thickness_y * wx[None, :] * dx / dy
    diagonal = np.zeros_like(a)
    diagonal[:, :-1] += x
    diagonal[:, 1:] += x
    diagonal[:-1] += y
    diagonal[1:] += y
    if np.any(active & (diagonal == 0)):
        raise ValueError("Neumann active domain contains an isolated node.")
    count = int(np.count_nonzero(active))
    indices = np.full(v.shape, -1, dtype=int)
    indices[active] = np.arange(count)
    rows = [np.arange(count)]
    cols = [np.arange(count)]
    data = [diagonal[active]]

    def add_faces(first, second, weights):
        connected = weights > 0
        p = indices[first][connected]
        q = indices[second][connected]
        c = weights[connected]
        rows.extend((p, q))
        cols.extend((q, p))
        data.extend((-c, -c))

    add_faces((slice(None), slice(None, -1)), (slice(None), slice(1, None)), x)
    add_faces((slice(None, -1), slice(None)), (slice(1, None), slice(None)), y)
    if cached is not None and cached[0] == key:
        solve = cached[1]
    else:
        matrix = coo_matrix(
            (np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
            shape=(count, count),
        ).tocsr()
        ones = csr_matrix(np.ones((count, 1)))
        augmented = bmat([[matrix, ones], [ones.T, None]], format="csc")
        solve = factorized(augmented)
        cache["neumann"] = (key, solve)
    rhs = np.concatenate((-integrated_source[active], [reference_value * count]))
    solved = solve(rhs)
    v[active] = solved[:-1]
    v[~active] = reference_value
    if not np.all(np.isfinite(v)):
        raise RuntimeError("Sparse Neumann solve produced nonfinite values.")
    defect = _neighbors(v, x, y) - diagonal * v - integrated_source
    residual = float(np.max(np.abs(defect[active] / diagonal[active]), initial=0.0))
    return SolverResult(v, 1, residual, residual <= tolerance)


def _solve(solution, coefficient, source, boundary_mask, boundary_values,
           dx, dy, tolerance, max_iterations, omega, red_black,
           residual_check_interval):
    v = np.asarray(solution, dtype=float).copy()
    a, b = np.asarray(coefficient), np.asarray(source)
    fixed = np.asarray(boundary_mask, dtype=bool)
    values = np.asarray(boundary_values)
    if v.ndim != 2 or min(v.shape) < 2:
        raise ValueError('Poisson grid must be at least 2 by 2.')
    if any(arr.shape != v.shape for arr in (a, b, fixed, values)):
        raise ValueError('Poisson arrays must have identical shapes.')
    if not all(np.all(np.isfinite(arr)) for arr in (v, a, b, values)):
        raise ValueError('Poisson data must be finite.')
    if np.any(a < 0) or not np.isfinite(dx + dy) or min(dx, dy) <= 0:
        raise ValueError('Conductivity must be nonnegative and spacing positive.')
    if not np.isfinite(tolerance) or tolerance <= 0 or max_iterations < 1:
        raise ValueError('Invalid Poisson convergence controls.')
    if residual_check_interval < 1:
        raise ValueError('Residual check interval must be positive.')
    if not np.isfinite(omega) or not 0 < omega < 2:
        raise ValueError('SOR omega must lie between 0 and 2.')
    x, y, diagonal = _system(a, dx, dy)
    active = ~fixed & (diagonal > 0)
    if np.any((~fixed & (diagonal == 0)) & (b != 0)):
        raise ValueError('Nonzero source in an electrically isolated cell.')
    v[fixed] = values[fixed]
    colors = np.indices(v.shape).sum(axis=0) % 2
    masks = [active & (colors == color) for color in (0, 1)]
    points = np.argwhere(active) if not red_black else None
    residual = 0.0
    for iteration in range(1, max_iterations + 1):
        if red_black:
            for mask in masks:
                rhs = _neighbors(v, x, y) - b
                v[mask] += omega*(rhs[mask]/diagonal[mask] - v[mask])
        else:
            for i, j in points:
                rhs = -b[i, j]
                if j: rhs += x[i, j-1]*v[i, j-1]
                if j+1 < v.shape[1]: rhs += x[i, j]*v[i, j+1]
                if i: rhs += y[i-1, j]*v[i-1, j]
                if i+1 < v.shape[0]: rhs += y[i, j]*v[i+1, j]
                v[i, j] += omega*(rhs/diagonal[i, j] - v[i, j])
        if not np.all(np.isfinite(v)):
            raise RuntimeError('Poisson iteration produced nonfinite values.')
        if (
            iteration == 1
            or iteration % residual_check_interval == 0
            or iteration == max_iterations
        ):
            defect = _neighbors(v, x, y) - diagonal*v - b
            residual = float(np.max(np.abs(defect[active]/diagonal[active]), initial=0))
            if residual <= tolerance:
                return SolverResult(v, iteration, residual, True)
    return SolverResult(v, max_iterations, residual, False)


def red_black_sor(solution, coefficient, source, boundary_mask, boundary_values,
                  dx, dy,
                  tolerance=_ELECTRICAL_SOLVER_DEFAULTS["tolerance"],
                  max_iterations=_ELECTRICAL_SOLVER_DEFAULTS["max_iterations"],
                  omega=_ELECTRICAL_SOLVER_DEFAULTS["omega"],
                  residual_check_interval=_ELECTRICAL_SOLVER_DEFAULTS["residual_check_interval"]):
    return _solve(solution, coefficient, source, boundary_mask, boundary_values,
                  dx, dy, tolerance, max_iterations, omega, True,
                  residual_check_interval)


def gauss_seidel(solution, coefficient, source, boundary_mask, boundary_values,
                 dx, dy,
                 tolerance=_ELECTRICAL_SOLVER_DEFAULTS["tolerance"],
                 max_iterations=_ELECTRICAL_SOLVER_DEFAULTS["max_iterations"],
                 omega=1.0,
                 residual_check_interval=_ELECTRICAL_SOLVER_DEFAULTS["residual_check_interval"]):
    return _solve(solution, coefficient, source, boundary_mask, boundary_values,
                  dx, dy, tolerance, max_iterations, omega, False,
                  residual_check_interval)
