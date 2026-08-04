"""
Iterative numerical solvers.

Current:

- Gauss-Seidel

Future:

- Successive over relaxation
- Conjugate gradient
- Multigrid
"""


import numpy as np



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
):
    """
    Solve variable coefficient elliptic equation.

    Currently solves:

        ∇ · (σ ∇V)=source

    using Gauss-Seidel iteration.

    Parameters
    ----------

    solution:
        Initial guess array.

    coefficient:
        Spatial coefficient σ.

    source:
        Right hand side.

    boundary_mask:
        Boolean array marking fixed cells.

    boundary_values:
        Values applied at boundaries.

    """

    V = solution.copy()

# Apply initial boundary conditions
    V[boundary_mask] = boundary_values[boundary_mask]


    dx2 = dx**2
    dy2 = dy**2



    for iteration in range(max_iterations):

        old = V.copy()


        for y in range(1, V.shape[0]-1):

            for x in range(1, V.shape[1]-1):


                if boundary_mask[y,x]:
                    V[y,x] = boundary_values[y,x]
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


                V[y,x] = (
                    numerator /
                    denominator
                )



# Reapply boundary conditions
        V[boundary_mask] = boundary_values[boundary_mask]


        error = np.max(
            np.abs(
                V-old
    )
)


        if error < tolerance:
            break



    return V