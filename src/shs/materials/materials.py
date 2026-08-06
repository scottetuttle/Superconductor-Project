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


    # GL parameters

    gl_alpha: float = 1.0

    gl_beta: float = 1.0

    tdgl_u: float = 5.79


    # Thermal properties

    thermal_conductivity: float = 0.0

    heat_capacity: float = 0.0


    # Electrical properties

    normal_resistivity: float = 0.0


    # Geometry

    thickness: float = 0.0