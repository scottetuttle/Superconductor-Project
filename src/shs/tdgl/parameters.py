"""
Dimensionless TDGL parameter definitions.

This module converts material information into
parameters used by the normalized TDGL equations.
"""


from dataclasses import dataclass


@dataclass
class TDGLParameters:
    """
    Dimensionless TDGL parameters.

    These parameters define the superconducting
    dynamics independent of SI units.
    """


    # Relaxation parameter
    u: float = 5.79


    # Coupling between amplitude and phase dynamics
    gamma: float = 0.0


    # GL parameter:
    #
    # kappa = lambda / xi
    #
    kappa: float = 1.0


    # Temperature normalization
    #
    # Reduced temperature:
    #
    # t = T / Tc
    #
    reduced_temperature: float = 0.0