import time
import numpy as np

from shs.numerics import (
    gauss_seidel,
    red_black_sor,
)


def create_test_problem():

    shape = (20, 20)

    solution = np.zeros(shape)

    coefficient = np.ones(shape)

    source = np.zeros(shape)

    boundary_mask = np.zeros(
        shape,
        dtype=bool
    )

    boundary_values = np.zeros(shape)


    # Left voltage boundary

    boundary_mask[:, 0] = True
    boundary_values[:, 0] = 1.0


    # Right voltage boundary

    boundary_mask[:, -1] = True
    boundary_values[:, -1] = 0.0


    return (
        solution,
        coefficient,
        source,
        boundary_mask,
        boundary_values,
    )



def test_solver_comparison():

    problem = create_test_problem()


    start = time.perf_counter()

    gs = gauss_seidel(
        *problem,
        dx=1.0,
        dy=1.0,
        max_iterations=50000,
    )

    gs_time = time.perf_counter() - start



    problem = create_test_problem()


    start = time.perf_counter()

    rb = red_black_sor(
        *problem,
        dx=1.0,
        dy=1.0,
    )

    rb_time = time.perf_counter() - start



    print()

    print(
        "Gauss-Seidel:"
    )

    print(
        f"Iterations: {gs.iterations}"
    )

    print(
        f"Runtime: {gs_time:.3f}s"
    )


    print()

    print(
        "Red-Black SOR:"
    )

    print(
        f"Iterations: {rb.iterations}"
    )

    print(
        f"Runtime: {rb_time:.3f}s"
    )


    assert gs.residual < 1e-7

    assert rb.residual < 1e-7


    difference = np.max(
        np.abs(
            gs.field -
            rb.field
        )
    )


    print(
        f"Difference: {difference:.3e}"
)

    print(
        f"GS residual: {gs.residual:.3e}"
)

    print(
        f"RB residual: {rb.residual:.3e}"
)


    assert difference < 1e-4