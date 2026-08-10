"""
Dimensionless TDGL parameter definitions.

This module defines the parameters used by the normalized
Time-Dependent Ginzburg-Landau model.

The current TDGL formulation is dimensionless and uses:

    u dpsi/dt =
        D^2 psi
        + (1 - T/Tc) psi
        - |psi|^2 psi

where:

    D = nabla - i A

The order parameter is normalized such that:

    |psi| = 1

at zero temperature in the absence of fields.
"""

from dataclasses import dataclass


@dataclass
class TDGLParameters:
    """
    Parameters for the normalized TDGL model.

    Parameters
    ----------
    u:
        TDGL relaxation parameter.

    gamma:
        Optional amplitude/phase coupling parameter.
        The current normalized implementation does not
        yet use this parameter dynamically.

    kappa:
        Ginzburg-Landau parameter:

            kappa = lambda / xi

        This parameter will become important when magnetic
        self-consistency is implemented.

    reduced_temperature:
        Default reduced temperature T/Tc.
        This is primarily provided for standalone model
        calculations and tests. The full solver obtains
        temperature from the simulation fields.
    """

    u: float = 5.79

    gamma: float = 0.0

    kappa: float = 1.0

    reduced_temperature: float = 0.0

    def validate(self):
        """
        Validate TDGL parameters.

        Raises
        ------
        ValueError
            If a parameter is physically invalid.
        """

        if self.u <= 0.0:
            raise ValueError(
                "TDGL relaxation parameter u must be positive."
            )

        if self.kappa <= 0.0:
            raise ValueError(
                "Ginzburg-Landau parameter kappa must be positive."
            )

        return True