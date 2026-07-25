"""
SHS Physics Fields

Stores physical quantities defined on the simulation mesh.
"""

from dataclasses import dataclass

import numpy as np

from shs.geometry.mesh import Mesh


@dataclass
class Fields:
    """
    Container for simulation fields.

    All arrays are defined on the simulation mesh.
    """

    temperature: np.ndarray

    voltage: np.ndarray

    current_density_x: np.ndarray

    current_density_y: np.ndarray

    heat_source: np.ndarray


    @classmethod
    def create(cls, mesh: Mesh, initial_temperature: float):
        """
        Create empty fields on a mesh.

        Parameters
        ----------
        mesh:
            Simulation mesh.

        initial_temperature:
            Starting temperature in Kelvin.
        """

        shape = (
            mesh.nx,
            mesh.ny
        )

        return cls(
            temperature=np.full(
                shape,
                initial_temperature
            ),

            voltage=np.zeros(shape),

            current_density_x=np.zeros(shape),

            current_density_y=np.zeros(shape),

            heat_source=np.zeros(shape)
        )