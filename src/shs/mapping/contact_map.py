"""
Contact mapping.

Converts geometry contacts into mesh cell masks.

Each contact is represented as a boolean mask over the mesh.

Future versions will support:

- Partial-edge contacts
- Interior contacts
- Arbitrary contact shapes
- Contact materials
"""

from dataclasses import dataclass

import numpy as np

from shs.geometry.geometry import Geometry
from shs.geometry.mesh import Mesh


@dataclass
class ContactMap:
    """
    Maps electrical contacts onto the numerical mesh.
    """

    contact_masks: dict[str, np.ndarray]


def build_contact_map(
    geometry: Geometry,
    mesh: Mesh,
) -> ContactMap:
    """
    Build mesh masks for every contact.

    Parameters
    ----------
    geometry
        Device geometry.

    mesh
        Numerical mesh.

    Returns
    -------
    ContactMap
    """

    shape = (mesh.ny, mesh.nx)

    contact_masks = {}

    for contact in geometry.contacts:

        mask = np.zeros(shape, dtype=bool)

        x0 = contact.x
        x1 = contact.x + contact.x_size
        y0 = contact.y
        y1 = contact.y + contact.y_size

        X, Y = np.meshgrid(
            mesh.x,
            mesh.y
        )

        mask = (
            (X >= x0)
            &
            (X <= x1)
            &
            (Y >= y0)
            &
            (Y <= y1)
        )

        contact_masks[contact.name] = mask

    return ContactMap(
        contact_masks=contact_masks
    )