"""
SHS Mesh Module

Creates numerical grids from physical geometries.
"""

from dataclasses import dataclass
import numpy as np

from .geometry import RectangularFilm


@dataclass
class Mesh:
    """
    Numerical representation of geometry.
    """

    x: np.ndarray
    y: np.ndarray

    nx: int
    ny: int


def create_mesh(
    geometry: RectangularFilm
) -> Mesh:
    """
    Generate a Cartesian mesh.

    Parameters
    ----------
    geometry:
        Physical geometry object.

    Returns
    -------
    Mesh
        Numerical grid.
    """

    x = np.linspace(
        0,
        geometry.width,
        geometry.nx
    )

    y = np.linspace(
        0,
        geometry.height,
        geometry.ny
    )


    return Mesh(
        x=x,
        y=y,
        nx=geometry.nx,
        ny=geometry.ny
    )