"""
Region mapping.

Converts geometric regions into a numerical map over the mesh.

Currently every cell belongs to the superconducting film.

Future versions will support:

- Contacts
- Substrates
- Defects
- Patterned devices
- Oxides
"""

from dataclasses import dataclass

import numpy as np

from shs.geometry.mesh import Mesh
from shs.geometry.geometry import Geometry


@dataclass
class RegionMap:
    """
    Numerical mapping of mesh cells to geometry regions.
    """

    region_ids: np.ndarray

    region_names: dict[int, str]


def build_region_map(
    geometry: Geometry,
    mesh: Mesh,
) -> RegionMap:
    """
    Build a region map from the device geometry.

    Parameters
    ----------
    geometry
        Device geometry.

    mesh
        Numerical mesh.

    Returns
    -------
    RegionMap
    """

    ids = np.zeros(
        (mesh.ny, mesh.nx),
        dtype=np.int32
    )

    region_names = {
        0: "film"
    }

    return RegionMap(
        region_ids=ids,
        region_names=region_names
    )