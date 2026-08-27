"""
Iterative numerical solvers.

Contains reusable PDE solvers.

Current:

- Gauss-Seidel
- Red-Black SOR

Future:

- Conjugate gradient
- Multigrid
- FEM/FVM methods
"""


from dataclasses import dataclass

import numpy as np



@dataclass
class SolverResult:
    """
    Result from an iterative solver.
    """

    field: np.ndarray

    iterations: int

    residual: float

    converged: bool



def red_black_sor(
    solution,
    coefficient,
    source,
    boundary_mask,
    boundary_values,
    dx,
    dy,
    tolerance=1e-8,
    max_iterations=10000,
    omega=1.7,
):
    """
    Solve:

        ∇ · (σ∇V)=source


    using Red-Black Successive
    Over Relaxation.


    Red-black ordering allows
    vectorized updates while keeping
    Gauss-Seidel convergence behavior.
    """


    V = solution.copy()


    V[boundary_mask] = (
        boundary_values[boundary_mask]
    )


    dx2 = dx**2
    dy2 = dy**2


    ny, nx = V.shape


    #
    # Create checkerboard masks
    #

    y_grid, x_grid = np.indices(
        V.shape
    )


    red = (
        (x_grid + y_grid) % 2 == 0
    )


    black = ~red


    red &= ~boundary_mask
    black &= ~boundary_mask



    residual = np.inf



    for iteration in range(
        1,
        max_iterations + 1
    ):

        old = V.copy()


        #
        # Update red cells
        #

        _sor_update(
            V,
            red,
            coefficient,
            source,
            dx2,
            dy2,
            omega,
        )


        #
        # Update black cells
        #

        _sor_update(
            V,
            black,
            coefficient,
            source,
            dx2,
            dy2,
            omega,
        )


        #
        # Reinforce boundaries
        #

        V[boundary_mask] = (
            boundary_values[boundary_mask]
        )



        residual = np.max(
            np.abs(
                V - old
            )
        )


        if residual < tolerance:

            return SolverResult(
                field=V,
                iterations=iteration,
                residual=residual,
                converged=True,
            )



    return SolverResult(
        field=V,
        iterations=max_iterations,
        residual=residual,
        converged=False,
    )




def _sor_update(
    V,
    mask,
    coefficient,
    source,
    dx2,
    dy2,
    omega,
):
    """
    Perform one vectorized SOR color update.

    The coefficient represents normal conductivity.

    Cells with zero conductivity are treated as insulating
    for the normal-current problem and are not updated.

    Face conductivities use the harmonic mean. Therefore a
    face between a conducting cell and a zero-conductivity
    cell has zero conductivity and cannot carry normal current.
    """

    #
    # Interior cells belonging to this SOR color
    #

    interior = mask[1:-1, 1:-1]

    if not np.any(interior):
        return

    #
    # Center conductivity
    #

    sigma_c = coefficient[1:-1, 1:-1]

    #
    # Neighbor conductivities
    #

    sigma_e = coefficient[1:-1, 2:]
    sigma_w = coefficient[1:-1, :-2]
    sigma_n = coefficient[2:, 1:-1]
    sigma_s = coefficient[:-2, 1:-1]

    #
    # Harmonic face conductivity
    #
    # sigma_face = 2*sigma_a*sigma_b / (sigma_a + sigma_b)
    #
    # When either side is zero, the face conductivity is zero.
    #

    denominator_e = sigma_c + sigma_e
    denominator_w = sigma_c + sigma_w
    denominator_n = sigma_c + sigma_n
    denominator_s = sigma_c + sigma_s

    ce = np.zeros_like(sigma_c)
    cw = np.zeros_like(sigma_c)
    cn = np.zeros_like(sigma_c)
    cs = np.zeros_like(sigma_c)

    np.divide(
        2.0 * sigma_c * sigma_e,
        denominator_e,
        out=ce,
        where=denominator_e > 0.0,
    )

    np.divide(
        2.0 * sigma_c * sigma_w,
        denominator_w,
        out=cw,
        where=denominator_w > 0.0,
    )

    np.divide(
        2.0 * sigma_c * sigma_n,
        denominator_n,
        out=cn,
        where=denominator_n > 0.0,
    )

    np.divide(
        2.0 * sigma_c * sigma_s,
        denominator_s,
        out=cs,
        where=denominator_s > 0.0,
    )

    #
    # Construct the discretized equation:
    #
    # ∇ · (σ ∇V) = source
    #

    numerator = (
        ce * V[1:-1, 2:] / dx2
        +
        cw * V[1:-1, :-2] / dx2
        +
        cn * V[2:, 1:-1] / dy2
        +
        cs * V[:-2, 1:-1] / dy2
        -
        source[1:-1, 1:-1]
    )

    denominator = (
        (ce + cw) / dx2
        +
        (cn + cs) / dy2
    )

    #
    # Only cells that actually participate in the
    # normal-conductivity problem should be updated.
    #
    # A zero denominator means there is no conducting
    # connection to any neighboring cell.
    #

    valid = (
        interior
        &
        (sigma_c > 0.0)
        &
        (denominator > 0.0)
    )

    if not np.any(valid):
        return

    #
    # Calculate the local SOR solution only where
    # the discretized equation is well-defined.
    #

    new_values = np.zeros_like(
        denominator
    )

    np.divide(
        numerator,
        denominator,
        out=new_values,
        where=denominator > 0.0,
    )

    #
    # SOR relaxation
    #

    current = V[1:-1, 1:-1]

    current[valid] += omega * (
        new_values[valid]
        -
        current[valid]
    )
def gauss_seidel(
    solution,
    coefficient,
    source,
    boundary_mask,
    boundary_values,
    dx,
    dy,
    tolerance=1e-8,
    max_iterations=10000,
    omega=1.0,
):
    """
    Reference Gauss-Seidel/SOR solver.

    Kept as a baseline solver for validation.
    """

    V = solution.copy()

    V[boundary_mask] = boundary_values[boundary_mask]

    dx2 = dx**2
    dy2 = dy**2

    residual = np.inf


    for iteration in range(1, max_iterations + 1):

        old = V.copy()

        for y in range(1, V.shape[0]-1):

            for x in range(1, V.shape[1]-1):

                if boundary_mask[y,x]:
                    continue


                sigma_e = (
                    coefficient[y,x]
                    +
                    coefficient[y,x+1]
                ) / 2

                sigma_w = (
                    coefficient[y,x]
                    +
                    coefficient[y,x-1]
                ) / 2

                sigma_n = (
                    coefficient[y,x]
                    +
                    coefficient[y+1,x]
                ) / 2

                sigma_s = (
                    coefficient[y,x]
                    +
                    coefficient[y-1,x]
                ) / 2


                denominator = (
                    (sigma_e + sigma_w)/dx2
                    +
                    (sigma_n + sigma_s)/dy2
                )


                numerator = (
                    sigma_e*V[y,x+1]/dx2
                    +
                    sigma_w*V[y,x-1]/dx2
                    +
                    sigma_n*V[y+1,x]/dy2
                    +
                    sigma_s*V[y-1,x]/dy2
                    -
                    source[y,x]
                )


                new_value = numerator / denominator


                V[y,x] += omega * (
                    new_value - V[y,x]
                )


        V[boundary_mask] = (
            boundary_values[boundary_mask]
        )


        residual = np.max(
            np.abs(V-old)
        )


        if residual < tolerance:

            return SolverResult(
                field=V,
                iterations=iteration,
                residual=residual,
                converged=True,
            )


    return SolverResult(
        field=V,
        iterations=max_iterations,
        residual=residual,
        converged=False,
    )