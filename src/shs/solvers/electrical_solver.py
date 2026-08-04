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
    voltage_left: float = 1.0,
    voltage_right: float = 0.0,
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

    sigma = (
        material_map.electrical_conductivity
    )



    #
    # Solve:
    #
    # ∇ · σ∇V = 0
    #

    result = red_black_sor(
        solution=V,
        coefficient=sigma,
        source=np.zeros_like(V),
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

    fields.current_density_x = (
        sigma *
        Ex
    )


    fields.current_density_y = (
        sigma *
        Ey
    )
        #
    # Joule heating
    #
    # Q = J^2 rho
    #

    J_squared = (
        fields.current_density_x**2
        +
        fields.current_density_y**2
    )


    rho = (
        material_map.normal_resistivity
    )


    fields.heat_source = (
        J_squared *
        rho
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