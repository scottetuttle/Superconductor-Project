"""
SHS Physics Fields

Stores all evolving physical quantities defined on the simulation mesh.

Fields represent simulation state.

Material properties do not belong here.
They are stored in MaterialMap.
"""

from dataclasses import dataclass

import numpy as np

from shs.geometry.mesh import Mesh

from shs.tdgl.initialization import equilibrium_superconducting_state


@dataclass
class Fields:
    """
    Container for evolving simulation fields.

    All arrays are defined on the simulation mesh.
    """

    # Thermal

    temperature: np.ndarray

    heat_source: np.ndarray


    # Electrical

    voltage: np.ndarray

    electric_field_x: np.ndarray
    electric_field_y: np.ndarray

    current_density_x: np.ndarray
    current_density_y: np.ndarray


    # Magnetic

    magnetic_field_x: np.ndarray
    magnetic_field_y: np.ndarray

    vector_potential_x: np.ndarray
    vector_potential_y: np.ndarray

    #tdgl

    psi: np.ndarray
    
    supercurrent_density_x: np.ndarray
    supercurrent_density_y: np.ndarray


    @classmethod
    def create(
        cls,
        mesh: Mesh,
        initial_temperature: float,
        initial_psi=None,
    ):
        """
        Create empty simulation fields.
        """

        shape = (
            mesh.ny,
            mesh.nx
        )

        if initial_psi is None:

            initial_psi = equilibrium_superconducting_state(
                shape,
                reduced_temperature=0.0,
            )
        return cls(

            # Thermal

            temperature=np.full(
                shape,
                initial_temperature,
                dtype=float
            ),

            heat_source=np.zeros(
                shape,
                dtype=float
            ),


            # Electrical

            voltage=np.zeros(
                shape,
                dtype=float
            ),

            electric_field_x=np.zeros(
                shape,
                dtype=float
            ),

            electric_field_y=np.zeros(
                shape,
                dtype=float
            ),

            current_density_x=np.zeros(
                shape,
                dtype=float
            ),

            current_density_y=np.zeros(
                shape,
                dtype=float
            ),


            # Magnetic

            magnetic_field_x=np.zeros(
                shape,
                dtype=float
            ),

            magnetic_field_y=np.zeros(
                shape,
                dtype=float
            ),

            vector_potential_x=np.zeros(
                shape,
                dtype=float
            ),

            vector_potential_y=np.zeros(
                shape,
                dtype=float
            ),

            #tdgl

            psi=np.asarray(
                initial_psi,
                dtype=complex,
            ).copy(),
            supercurrent_density_x=np.zeros(
                shape,
                dtype=float
            ),
            supercurrent_density_y=np.zeros(
                shape,
                dtype=float
            ),         
        )