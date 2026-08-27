"""
Electrical transport solver.

Solves:

    ∇ · (σ ∇V) = 0

then calculates:

    E = -∇V

    J = σE


Uses iterative numerical solvers
through the SHS numerics layer.
"""


import numpy as np


from shs.physics import Fields

from shs.mapping.contact_map import ContactMap

from shs.mapping.material_map import MaterialMap

from shs.numerics import red_black_sor



def electrical_step(
    fields: Fields,
    mesh,
    material_map: MaterialMap,
    contact_map: ContactMap,
    voltage_left: float = 1e-3,
    voltage_right: float = 0.0,
    superconducting_fraction=None,
    superconducting_current_x=None,
    superconducting_current_y=None,
):
    """
    Perform one electrical transport solve.

    Steps:

    1. Build voltage boundary conditions
    2. Solve conductivity PDE
    3. Calculate electric field
    4. Calculate current density

    """


    V = fields.voltage.copy()


    #
    # Contact masks
    #

    left = contact_map.contact_masks.get(
        "left_current"
    )

    right = contact_map.contact_masks.get(
        "right_current"
    )


    if left is None or right is None:

        raise ValueError(
            "Current contacts missing."
        )



    #
    # Boundary conditions
    #

    boundary_mask = (
        left |
        right
    )


    boundary_values = np.zeros_like(
        V
    )


    boundary_values[left] = (
        voltage_left
    )

    boundary_values[right] = (
        voltage_right
    )



    #
    # Electrical conductivity
    #

    if superconducting_current_x is None:
        superconducting_current_x = np.zeros_like(V)

    if superconducting_current_y is None:
        superconducting_current_y = np.zeros_like(V)

    if superconducting_fraction is None:
        superconducting_fraction = np.zeros_like(V)

    normal_fraction = np.maximum(
        1.0 - superconducting_fraction,
        0.0,
    )

    sigma = (
        material_map.electrical_conductivity
        *
        normal_fraction
    )
    #
# Perfectly superconducting limit
#
# If there is no normal conductivity anywhere,
# there is no normal-current potential equation
# to solve.
#

    if np.all(sigma <= 0.0):

        fields.electric_field_x = np.zeros_like(V)
        fields.electric_field_y = np.zeros_like(V)

        fields.normal_current_density_x = (
            np.zeros_like(V)
        )

        fields.normal_current_density_y = (
            np.zeros_like(V)
        )

        fields.current_density_x = (
            superconducting_current_x.copy()
        )

        fields.current_density_y = (
            superconducting_current_y.copy()
        )

        fields.heat_source = np.zeros_like(V)

        fields.electrical_solver_iterations = 0
        fields.electrical_solver_residual = 0.0

        return fields





    #
    # Solve:
    #
    # ∇ · σ∇V = 0
    #
    # J=Js​+Jn
    #
    # ∇⋅(Js​+σn​E)=0.
    #

    source = np.zeros_like(V)

    source[1:-1, 1:-1] = -(
        (
            superconducting_current_x[1:-1, 2:]
            -
            superconducting_current_x[1:-1, :-2]
        )
        /
        (2.0 * mesh.dx)
        +
        (
            superconducting_current_y[2:, 1:-1]
            -
            superconducting_current_y[:-2, 1:-1]
        )
        /
        (2.0 * mesh.dy)
    )

    result = red_black_sor(
        solution=V,
        coefficient=sigma,
        source=source,
        boundary_mask=boundary_mask,
        boundary_values=boundary_values,
        dx=mesh.dx,
        dy=mesh.dy,
        tolerance=1e-8,
        max_iterations=10000,
        omega=1.7,
    )



    fields.voltage = (
        result.field
    )



    #
    # Electric field
    #

    V = fields.voltage


    Ex = np.zeros_like(V)

    Ey = np.zeros_like(V)



    Ex[:,1:-1] = -(
        V[:,2:]
        -
        V[:,:-2]
    ) / (
        2 * mesh.dx
    )


    Ey[1:-1,:] = -(
        V[2:,:]
        -
        V[:-2,:]
    ) / (
        2 * mesh.dy
    )



    fields.electric_field_x = Ex

    fields.electric_field_y = Ey



    #
    # Current density
    #

    normal_current_x = (
        sigma *
        Ex
    )

    normal_current_y = (
        sigma *
        Ey
    )
    fields.normal_current_density_x = (
        normal_current_x
    )

    fields.normal_current_density_y = (
        normal_current_y
    )    


    fields.current_density_x = (
        normal_current_x
        +
        superconducting_current_x
    )

    fields.current_density_y = (
        normal_current_y
        +
        superconducting_current_y
    )
        #
    # Joule heating
    #
    # Q = J^2 rho
    #

    fields.heat_source = (
        normal_current_x * Ex
        +
        normal_current_y * Ey
    )


    #
    # Store solver information
    #
    # Future:
    # add to SimulationState
    #

    fields.electrical_solver_iterations = (
        result.iterations
    )


    fields.electrical_solver_residual = (
        result.residual
    )



    return fields