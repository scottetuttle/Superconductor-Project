"""
TDGL physical model.

Defines the dimensionless time-dependent
Ginzburg-Landau parameters and equations.

The numerical solver is responsible for
advancing the solution in time.
"""

from dataclasses import dataclass

from .parameters import TDGLParameters


@dataclass
class TDGLModel:
    """
    Dimensionless TDGL physics model.
    """

    parameters: TDGLParameters


    def reduced_temperature(
        self,
        temperature,
        critical_temperature,
    ):
        """
        Calculate reduced temperature.

        T_red = T / Tc

        Values:

        T_red < 1:
            superconducting regime

        T_red >= 1:
            normal regime
        """

        return temperature / critical_temperature



    def equilibrium_amplitude(
        self,
        reduced_temperature,
    ):
        """
        Equilibrium order parameter amplitude.

        Simplified GL result:

        |psi|^2 = 1 - T/Tc

        """

        if reduced_temperature >= 1:
            return 0.0

        return (
            1.0 -
            reduced_temperature
        )