"""
SHS Mesh Module

Converts physical geometry into numerical grids.
"""

from dataclasses import dataclass
import numpy as np


from .geometry import RectangularFilm



@dataclass
class Mesh:

    x: np.ndarray
    y: np.ndarray

    nx: int
    ny: int

    dx: float
    dy: float



def create_mesh(geometry):

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
        ny=geometry.ny,

        dx=geometry.dx,
        dy=geometry.dy
    )