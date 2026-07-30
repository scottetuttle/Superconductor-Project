"""
SHS Electrical Physics Module

Defines electrical transport properties.

Current model:

Ohmic transport:

    J = sigma E

where:

    sigma = electrical conductivity
    J     = current density
    E     = electric field


Future:

- superconducting conductivity
- two-fluid model
- nonlinear resistivity
- vortex dissipation
- TDGL coupling
"""


from dataclasses import dataclass



@dataclass
class ElectricalModel:
    """
    Electrical transport parameters.
    """


    reference_voltage: float = 1.0


    def conductivity(
        self,
        resistivity
    ):
        """
        Convert resistivity to conductivity.

        sigma = 1/rho
        """

        if resistivity == 0:
            raise ValueError(
                "Resistivity cannot be zero."
            )


        return 1.0 / resistivity