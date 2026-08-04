"""
Electrical transport solver.

Solves:

    ∇ · (σ ∇V) = 0

where:

    σ = electrical conductivity
    V = electric potential

Then calculates:

    E = -∇V

    J = σE

    Q = J · E
"""


import numpy as np


from shs.physics import Fields
from shs.physics.electrical import ElectricalModel

from shs.numerics import (
    gradient,
    gauss_seidel,
)

from shs.mapping.contact_map import ContactMap
from shs.mapping.material_map import MaterialMap



def electrical_step(
    fields: Fields,
    mesh,
    material_map: MaterialMap,
    contact_map: ContactMap,
    voltage_left: float = 1.0,
    voltage_right: float = 0.0,
):
    """
    Solve electrical transport.

    Uses:

        ∇ · (σ∇V)=0

    """

    #
    # Initial voltage guess
    #

    V = fields.voltage.copy()


    #
    # Build voltage boundary arrays
    #

    boundary_mask = np.zeros_like(
        V,
        dtype=bool
    )


    boundary_values = np.zeros_like(
        V,
        dtype=float
    )


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


    boundary_mask[left] = True
    boundary_values[left] = voltage_left


    boundary_mask[right] = True
    boundary_values[right] = voltage_right



    #
    # Conductivity map
    #

    sigma = material_map.electrical_conductivity



    #
    # Solve potential equation
    #

    result = gauss_seidel(
        solution=V,
        coefficient=sigma,
        source=np.zeros_like(V),
        boundary_mask=boundary_mask,
        boundary_values=boundary_values,
        dx=mesh.dx,
        dy=mesh.dy,
        omega=1.7,
)


    V = result.field

    print(
    f"Electrical solver: "
    f"{result.iterations} iterations, "
    f"residual={result.residual:.3e}, "
    f"converged={result.converged}"
)


    fields.voltage = V



    #
    # Electric field
    #

    Ex, Ey = gradient(
        V,
        mesh.dx,
        mesh.dy,
    )


    Ex = -Ex
    Ey = -Ey


    fields.electric_field_x = Ex
    fields.electric_field_y = Ey



    #
    # Current density
    #

    Jx = sigma * Ex
    Jy = sigma * Ey


    fields.current_density_x = Jx
    fields.current_density_y = Jy



    #
    # Joule heating
    #

    electrical_model = ElectricalModel()


    fields.heat_source = electrical_model.joule_heating(
        Jx,
        Jy,
        Ex,
        Ey,
    )


    return fields