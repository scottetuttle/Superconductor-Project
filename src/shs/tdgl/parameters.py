"""
Dimensionless TDGL parameter definitions.

This module defines the parameters used by the normalized
Time-Dependent Ginzburg-Landau model.

The current TDGL formulation is dimensionless and uses:

    u (d/dt + i phi) psi =
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

from shs.utils.defaults import default_section


_DEFAULTS = default_section("tdgl")


@dataclass
class TDGLParameters:
    """
    Parameters for the normalized TDGL model.

    Parameters
    ----------
    u:
        TDGL relaxation parameter.

    gamma:
        Kramer-Watts-Tobin amplitude/phase coupling parameter. Zero recovers
        the simpler relaxational TDGL equation.

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

    u: float = _DEFAULTS["u"]
    time_integrator: str = _DEFAULTS["time_integrator"]

    gamma: float = _DEFAULTS["gamma"]

    kappa: float = _DEFAULTS["kappa"]

    normalization: str = _DEFAULTS["normalization"]

    temperature_model: str = _DEFAULTS["temperature_model"]

    include_scalar_potential: bool = _DEFAULTS["include_scalar_potential"]

    reduced_temperature: float = 0.0

    max_normalized_timestep: float = _DEFAULTS["max_normalized_timestep"]

    stability_safety_factor: float = _DEFAULTS["stability_safety_factor"]

    def validate(self):
        """
        Validate TDGL parameters.

        Raises
        ------
        ValueError
            If a parameter is physically or numerically invalid.
        """

        if self.time_integrator not in {"euler", "heun"}:
            raise ValueError("TDGL time_integrator must be euler or heun.")
        if self.u <= 0.0:
            raise ValueError(
                "TDGL relaxation parameter u must be positive."
            )

        if self.kappa <= 0.0:
            raise ValueError(
                "Ginzburg-Landau parameter kappa must be positive."
            )

        if self.normalization not in {"pytdgl", "legacy_gl"}:
            raise ValueError("TDGL normalization must be pytdgl or legacy_gl.")

        if self.temperature_model not in {
            "tc_over_t_minus_one", "one_minus_t_over_tc"
        }:
            raise ValueError("Unknown TDGL temperature coefficient model.")

        if self.max_normalized_timestep <= 0.0:
            raise ValueError(
                "Maximum normalized TDGL timestep must be positive."
            )

        if not isinstance(self.include_scalar_potential, bool):
            raise ValueError("include_scalar_potential must be boolean.")

        if not 0.0 < self.stability_safety_factor <= 1.0:
            raise ValueError("TDGL stability_safety_factor must lie in (0, 1].")

        return True
