"""
Electrical transport solver.

Solves:

    J = σE

using voltage gradients.
"""

import numpy as np

from shs.physics import Fields
from shs.physics.electrical import ElectricalModel

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

    Current implementation

    - linear voltage interpolation
    - Ohmic transport
    - Joule heating calculation
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

    #
    # Apply voltage contacts
    #

    V[left] = voltage_left
    V[right] = voltage_right

    #
    # Initial linear voltage profile
    #

    for i in range(mesh.nx):

        fraction = i / (mesh.nx - 1)

        V[:, i] = (
            voltage_left
            +
            fraction
            * (voltage_right - voltage_left)
        )

    fields.voltage = V

    #
    # Electric field
    #

    Ex = np.zeros_like(V)
    Ey = np.zeros_like(V)

    Ex[:, 1:-1] = -(
        V[:, 2:]
        -
        V[:, :-2]
    ) / (
        2 * mesh.dx
    )

    Ey[1:-1, :] = -(
        V[2:, :]
        -
        V[:-2, :]
    ) / (
        2 * mesh.dy
    )

    fields.electric_field_x = Ex
    fields.electric_field_y = Ey

    #
    # Conductivity
    #

    sigma = material_map.electrical_conductivity

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