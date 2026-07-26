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
        geometry.film.width,
        geometry.film.nx
    )

    y = np.linspace(
        0,
        geometry.film.height,
        geometry.film.ny
    )


    return Mesh(
        x=x,
        y=y,

        nx=geometry.film.nx,
        ny=geometry.film.ny,

        dx=geometry.film.nx,
        dy=geometry.film.ny
    )