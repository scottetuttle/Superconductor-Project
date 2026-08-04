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
    """


    sigma_e = (
        coefficient[:,1:]
        +
        coefficient[:,:-1]
    ) / 2


    sigma_w = sigma_e.copy()


    sigma_n = (
        coefficient[1:,:]
        +
        coefficient[:-1,:]
    ) / 2


    sigma_s = sigma_n.copy()


    #
    # Interior slices
    #

    interior = mask[1:-1,1:-1]


    if not np.any(interior):
        return



    ce = (
        coefficient[1:-1,2:]
        +
        coefficient[1:-1,1:-1]
    ) / 2


    cw = (
        coefficient[1:-1,:-2]
        +
        coefficient[1:-1,1:-1]
    ) / 2


    cn = (
        coefficient[2:,1:-1]
        +
        coefficient[1:-1,1:-1]
    ) / 2


    cs = (
        coefficient[:-2,1:-1]
        +
        coefficient[1:-1,1:-1]
    ) / 2



    numerator = (

        ce *
        V[1:-1,2:]
        / dx2

        +

        cw *
        V[1:-1,:-2]
        / dx2

        +

        cn *
        V[2:,1:-1]
        / dy2

        +

        cs *
        V[:-2,1:-1]
        / dy2

        -

        source[1:-1,1:-1]

    )


    denominator = (

        (ce + cw)
        / dx2

        +

        (cn + cs)
        / dy2

    )


    new_values = (
        numerator /
        denominator
    )


    current = V[1:-1,1:-1]


    current[interior] += omega * (
        new_values[interior]
        -
        current[interior]
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