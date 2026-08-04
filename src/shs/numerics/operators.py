"""
Finite difference operators.

These operators are independent of SHS physics.

They operate on numerical fields only.

Future extensions:

- nonuniform grids
- higher order schemes
- spectral methods
- finite volume methods
"""


import numpy as np



def gradient(
    field,
    dx,
    dy,
):
    """
    Calculate spatial gradient.

    Returns:

    d(field)/dx
    d(field)/dy
    """

    gradient_x = np.zeros_like(field)

    gradient_y = np.zeros_like(field)


    gradient_x[:, 1:-1] = (
        field[:, 2:]
        -
        field[:, :-2]
    ) / (
        2 * dx
    )


    gradient_y[1:-1, :] = (
        field[2:, :]
        -
        field[:-2, :]
    ) / (
        2 * dy
    )


    return gradient_x, gradient_y



def divergence(
    field_x,
    field_y,
    dx,
    dy,
):
    """
    Calculate divergence:

        ∇ · F

    """

    result = np.zeros_like(field_x)


    result[:,1:-1] += (
        field_x[:,2:]
        -
        field_x[:,:-2]
    ) / (
        2 * dx
    )


    result[1:-1,:] += (
        field_y[2:,:]
        -
        field_y[:-2,:]
    ) / (
        2 * dy
    )


    return result



def laplacian(
    field,
    dx,
    dy,
):
    """
    Calculate scalar Laplacian:

        ∇²f
    """


    result = np.zeros_like(field)


    result[1:-1,1:-1] = (

        field[2:,1:-1]
        +
        field[:-2,1:-1]

        +

        field[1:-1,2:]
        +
        field[1:-1,:-2]

        -

        4 * field[1:-1,1:-1]

    ) / dx**2


    return result