"""
Iterative numerical solvers.

Current:

- Gauss-Seidel
- Successive Over Relaxation (SOR)

Future:

- Conjugate gradient
- Multigrid
- Newton methods
"""


from dataclasses import dataclass

import numpy as np



@dataclass
class SolverResult:
    """
    Result returned by iterative solvers.
    """

    field: np.ndarray

    iterations: int

    residual: float

    converged: bool



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
    Solve variable coefficient elliptic equation.

    Solves:

        ∇ · (σ ∇V) = source


    Parameters
    ----------

    omega:
        Relaxation factor.

        omega = 1.0:
            Gauss-Seidel

        1 < omega < 2:
            SOR

    """

    V = solution.copy()


    # Apply initial boundaries

    V[boundary_mask] = boundary_values[boundary_mask]


    dx2 = dx**2
    dy2 = dy**2


    residual = np.inf


    for iteration in range(1, max_iterations + 1):

        old = V.copy()


        for y in range(1, V.shape[0]-1):

            for x in range(1, V.shape[1]-1):


                if boundary_mask[y, x]:

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

                    (sigma_e + sigma_w)
                    / dx2

                    +

                    (sigma_n + sigma_s)
                    / dy2

                )


                numerator = (

                    sigma_e *
                    V[y,x+1]
                    / dx2

                    +

                    sigma_w *
                    V[y,x-1]
                    / dx2

                    +

                    sigma_n *
                    V[y+1,x]
                    / dy2

                    +

                    sigma_s *
                    V[y-1,x]
                    / dy2

                    -

                    source[y,x]

                )


                new_value = (
                    numerator /
                    denominator
                )


                #
                # SOR correction
                #

                V[y,x] = (

                    V[y,x]
                    +
                    omega *
                    (
                        new_value
                        -
                        V[y,x]
                    )

                )


        #
        # Enforce boundaries
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