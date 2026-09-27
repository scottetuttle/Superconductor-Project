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
from shs.utils.defaults import default_section

from shs.geometry.geometry import Geometry
from shs.geometry.mesh import Mesh


@dataclass
class ContactMap:
    """
    Maps electrical contacts onto the numerical mesh.
    """

    contact_masks: dict[str, np.ndarray]

    contact_types: dict[str, str]



def build_contact_map(
    geometry: Geometry,
    mesh: Mesh,
    roundoff_ulps=None,
) -> ContactMap:
    """
    Build mesh masks for every contact.
    """

    shape = (
        mesh.ny,
        mesh.nx
    )


    contact_masks = {}

    contact_types = {}


    X, Y = np.meshgrid(
        mesh.x,
        mesh.y
    )


    if roundoff_ulps is None:
        roundoff_ulps = default_section("numerics")["contact_roundoff_ulps"]
    tolerance = roundoff_ulps * np.spacing(max(abs(mesh.x[-1]), abs(mesh.y[-1])))

    for contact in geometry.contacts:


        mask = (
            (X >= contact.x - tolerance)
            &
            (X <= contact.x + contact.x_size + tolerance)

            &

            (Y >= contact.y - tolerance)
            &
            (Y <= contact.y + contact.y_size + tolerance)
        )


        if not np.any(mask):
            raise ValueError(
                f"Contact {contact.name} contains no mesh cells."
            )


        contact_masks[contact.name] = mask


        contact_types[contact.name] = (
            contact.contact_type
        )


    return ContactMap(

        contact_masks=contact_masks,

        contact_types=contact_types

    )
