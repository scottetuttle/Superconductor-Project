"""
Electrical transport solver.

Solves:

    J = sigma E

using voltage gradients.
"""


import numpy as np


from shs.physics import Fields
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
    Perform one electrical solve step.

    Current implementation:

    - linear voltage gradient
    - ohmic conductivity
    - contact driven transport

    """

    V = fields.voltage.copy()


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


    # Apply voltage boundary conditions

    V[left] = voltage_left

    V[right] = voltage_right


    # Simple initial potential interpolation

    for i in range(
        mesh.nx
    ):

        fraction = i / (
            mesh.nx - 1
        )

        V[:, i] = (
            voltage_left
            +
            fraction *
            (
                voltage_right
                -
                voltage_left
            )
        )


    fields.voltage = V


    # Electric field

    Ex = np.zeros_like(V)

    Ey = np.zeros_like(V)


    Ex[:,1:-1] = -(
        V[:,2:]
        -
        V[:,:-2]
    ) / (
        2 *
        mesh.dx
    )


    Ey[1:-1,:] = -(
        V[2:,:]
        -
        V[:-2,:]
    ) / (
        2 *
        mesh.dy
    )


    fields.electric_field_x = Ex

    fields.electric_field_y = Ey


    # Conductivity

    sigma = (
        1 /
        material_map.normal_resistivity
    )


    # Current density

    fields.current_density_x = (
        sigma *
        Ex
    )


    fields.current_density_y = (
        sigma *
        Ey
    )


    return fields