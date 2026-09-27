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

    region_types: dict[int, str]

    active_mask: np.ndarray | None = None




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

    active_mask = np.ones_like(ids, dtype=bool)
    xx, yy = np.meshgrid(mesh.x, mesh.y)
    hole_names = [hole.name for hole in geometry.holes or []]
    if len(hole_names) != len(set(hole_names)):
        raise ValueError("Geometry hole names must be unique.")
    for hole in geometry.holes or []:
        if hole.shape == "circle":
            if hole.radius is None or hole.radius <= 0:
                raise ValueError(f"Circular hole {hole.name!r} requires a positive radius.")
            half_width = half_height = hole.radius
            # Work in mesh units so roundoff at a grid-aligned circle does not
            # include one side of a symmetric boundary and exclude the other.
            radius_x = hole.radius / mesh.dx
            radius_y = hole.radius / mesh.dy
            distance_squared = (
                ((xx - hole.x) / mesh.dx) ** 2 / radius_x**2
                + ((yy - hole.y) / mesh.dy) ** 2 / radius_y**2
            )
            inside = distance_squared <= 1.0 + 64.0 * np.finfo(float).eps
        elif hole.shape == "rectangle":
            if hole.width is None or hole.height is None or min(hole.width, hole.height) <= 0:
                raise ValueError(f"Rectangular hole {hole.name!r} requires positive width and height.")
            half_width, half_height = 0.5 * hole.width, 0.5 * hole.height
            inside = ((np.abs(xx - hole.x) <= 0.5 * hole.width)
                      & (np.abs(yy - hole.y) <= 0.5 * hole.height))
        else:
            raise ValueError(f"Unsupported hole shape {hole.shape!r}.")
        if not (mesh.x[0] < hole.x - half_width
                and hole.x + half_width < mesh.x[-1]
                and mesh.y[0] < hole.y - half_height
                and hole.y + half_height < mesh.y[-1]):
            raise ValueError(f"Hole {hole.name!r} must lie strictly inside the film.")
        if np.any(inside[[0, -1], :]) or np.any(inside[:, [0, -1]]):
            raise ValueError(f"Hole {hole.name!r} must not touch the outer film boundary.")
        active_mask &= ~inside
    ids[~active_mask] = -1

    region_names = {
        -1: "vacuum", 0: "film"
}

    region_types = {
        -1: "inactive", 0: "superconductor"
}
    return RegionMap(
        region_ids=ids,
        region_names=region_names,
        region_types=region_types,
        active_mask=active_mask,
    )
