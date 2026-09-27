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

    All arrays are shaped (ny, nx). Temperature is K, voltage V,
    heat sources W/m^3, currents A/m^2, electric field V/m,
    magnetic field T, and vector potential T m. psi is dimensionless.
    Directional electrical/current arrays use outgoing links; the last
    x column and last y row have no outgoing link.
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
    magnetic_field_z: np.ndarray

    vector_potential_x: np.ndarray
    vector_potential_y: np.ndarray

    applied_vector_potential_x: np.ndarray
    applied_vector_potential_y: np.ndarray

    induced_vector_potential_x: np.ndarray
    induced_vector_potential_y: np.ndarray

    #tdgl

    psi: np.ndarray
    
    supercurrent_density_x: np.ndarray
    supercurrent_density_y: np.ndarray

    normal_current_density_x: np.ndarray
    normal_current_density_y: np.ndarray

    phonon_temperature: np.ndarray | None = None

    # Volumetric power densities [W/m^3]. heat_source is their total.
    external_heat_source: np.ndarray | None = None
    joule_heat_source: np.ndarray | None = None
    laser_heat_source: np.ndarray | None = None
    laser_position_x_m: float | None = None
    laser_position_y_m: float | None = None

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

            phonon_temperature=np.full(shape, initial_temperature, dtype=float),

            magnetic_field_z=np.zeros(
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

            applied_vector_potential_x=np.zeros(shape, dtype=float),
            applied_vector_potential_y=np.zeros(shape, dtype=float),
            induced_vector_potential_x=np.zeros(shape, dtype=float),
            induced_vector_potential_y=np.zeros(shape, dtype=float),

            #tdgl

            # TDGL

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

            normal_current_density_x=np.zeros(
                shape,
                dtype=float
            ),

            normal_current_density_y=np.zeros(
                shape,
                dtype=float
            ),      
        )
