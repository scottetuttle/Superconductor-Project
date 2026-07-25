from dataclasses import dataclass


@dataclass
class Material:
    """
    Defines physical properties of a simulated material.
    """

    name: str

    # Superconducting properties
    Tc: float
    coherence_length: float
    penetration_depth: float

    # Thermal properties
    thermal_conductivity: float
    heat_capacity: float

    # Electrical properties
    normal_resistivity: float

    # Geometry
    thickness: float