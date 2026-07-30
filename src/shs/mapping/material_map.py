"""
Material mapping.

Converts RegionMap into spatially varying material property arrays.

Current implementation assumes the entire superconducting film is one
material.

Future versions will support:

- Contacts
- Multiple superconductors
- Patterned devices
- Defects
- Oxides
- Substrates
"""

from dataclasses import dataclass

import numpy as np

from shs.materials.materials import Material
from shs.mapping.region_map import RegionMap


@dataclass
class MaterialMap:
    """
    Maps every mesh cell to material properties.
    """

    material_ids: np.ndarray

    materials: dict[int, Material]

    thermal_conductivity: np.ndarray
    heat_capacity: np.ndarray
    normal_resistivity: np.ndarray
    electrical_conductivity: np.ndarray

    thickness: np.ndarray

    Tc: np.ndarray
    coherence_length: np.ndarray
    penetration_depth: np.ndarray


def build_material_map(
    region_map: RegionMap,
    material: Material,
) -> MaterialMap:
    """
    Build a MaterialMap from a RegionMap.

    Parameters
    ----------
    region_map
        Region map of the device.

    material
        Material assigned to the superconducting film.

    Returns
    -------
    MaterialMap
    """

    shape = region_map.region_ids.shape

    material_ids = np.zeros(shape, dtype=np.int32)

    thermal_conductivity = np.full(
        shape,
        material.thermal_conductivity,
        dtype=float
    )

    heat_capacity = np.full(
        shape,
        material.heat_capacity,
        dtype=float
    )

    normal_resistivity = np.full(
        shape,
        material.normal_resistivity,
        dtype=float
    )
    electrical_conductivity = np.full(
        shape,
        1.0 / material.normal_resistivity,
        dtype=float,
    )

    thickness = np.full(
        shape,
        material.thickness,
        dtype=float
    )

    Tc = np.full(
        shape,
        material.Tc,
        dtype=float
    )

    coherence_length = np.full(
        shape,
        material.coherence_length,
        dtype=float
    )

    penetration_depth = np.full(
        shape,
        material.penetration_depth,
        dtype=float
    )

    return MaterialMap(
        material_ids=material_ids,
        materials={
            0: material
        },

        thermal_conductivity=thermal_conductivity,
        heat_capacity=heat_capacity,
        
        normal_resistivity=normal_resistivity,
        electrical_conductivity=electrical_conductivity,

        thickness=thickness,
        Tc=Tc,
        coherence_length=coherence_length,
        penetration_depth=penetration_depth,
    )