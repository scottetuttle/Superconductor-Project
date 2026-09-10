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
        The full solver obtains temperature from the
        simulation fields.

    max_normalized_timestep:
        Maximum normalized timestep used by the explicit
        TDGL integrator. Larger physical timesteps are
        automatically divided into internal substeps.
    """

    u: float = 5.79

    gamma: float = 0.0

    kappa: float = 1.0

    reduced_temperature: float = 0.0

    max_normalized_timestep: float = 0.1

    def validate(self):
        """
        Validate TDGL parameters.

        Raises
        ------
        ValueError
            If a parameter is physically or numerically invalid.
        """

        if self.u <= 0.0:
            raise ValueError(
                "TDGL relaxation parameter u must be positive."
            )

        if self.kappa <= 0.0:
            raise ValueError(
                "Ginzburg-Landau parameter kappa must be positive."
            )

        if self.max_normalized_timestep <= 0.0:
            raise ValueError(
                "Maximum normalized TDGL timestep must be positive."
            )

        return True