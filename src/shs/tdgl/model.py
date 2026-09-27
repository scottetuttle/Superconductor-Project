"""
TDGL physical model.

Defines the normalized Time-Dependent Ginzburg-Landau
physics used by SHS.

The supported dimensionless TDGL equation is:

    u (d/dt + i phi) psi =
        D^2 psi
        + epsilon(T) psi
        - |psi|^2 psi

where:

    D = nabla - i A

The default pyTDGL convention uses epsilon=clip(Tc/T-1, -1, 1),
tau0=mu0*sigma*lambda^2, and its matched electromagnetic scales. The former
1-T/Tc coefficient and microscopic relaxation time are available only through
the explicit legacy_gl configuration. The scalar-potential term is optional. Generalized-gamma dynamics are handled
by the time integrator; self-consistent vector-potential evolution is not yet
implemented.

The numerical solver is responsible for advancing
the solution in time.
"""

from dataclasses import dataclass
import numpy as np
from shs.utils.constants import (
    REDUCED_PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT,
    VACUUM_PERMEABILITY,
)

from .parameters import TDGLParameters

from .operators import gauge_covariant_gradient

@dataclass
class TDGLModel:
    """
    Normalized TDGL physics model.
    """

    parameters: TDGLParameters

    def __post_init__(self):
        self.parameters.validate()

    @staticmethod
    def reduced_temperature(
        temperature,
        critical_temperature,
    ):
        """
        Calculate reduced temperature.

        T_red = T / Tc

        Parameters
        ----------
        temperature:
            Physical temperature.

        critical_temperature:
            Critical temperature Tc.

        Returns
        -------
        float or ndarray
            Reduced temperature.
        """

        if np.any(critical_temperature <= 0.0):
            raise ValueError(
                "Critical temperature must be positive."
            )

        return (
            np.asarray(temperature) /
            np.asarray(critical_temperature)
        )

    def temperature_coefficient(self, reduced_temperature):
        """Return the configured dimensionless linear TDGL coefficient."""
        reduced = np.asarray(reduced_temperature, dtype=float)
        if self.parameters.temperature_model == "tc_over_t_minus_one":
            with np.errstate(divide="ignore"):
                coefficient = 1.0 / reduced - 1.0
            coefficient = np.clip(coefficient, -1.0, 1.0)
        else:
            coefficient = 1.0 - reduced
        return float(coefficient) if coefficient.ndim == 0 else coefficient

    def equilibrium_amplitude_squared(
        self, reduced_temperature,
    ):
        """
        Calculate equilibrium |psi|^2.

        For the normalized GL potential:

            V(psi) =
                -(1-t)|psi|^2
                + 1/2 |psi|^4

        the equilibrium solution is:

            |psi|^2 = 1 - t

        for:

            t < 1

        and:

            |psi|^2 = 0

        for:

            t >= 1.

        Parameters
        ----------
        reduced_temperature:
            T/Tc.

        Returns
        -------
        float or ndarray
            Equilibrium order-parameter magnitude squared.
        """

        t = np.asarray(
            reduced_temperature,
            dtype=float
        )

        result = np.maximum(self.temperature_coefficient(t), 0.0)

        if result.ndim == 0:
            return float(result)

        return result

    def equilibrium_amplitude(
        self, reduced_temperature,
    ):
        """
        Calculate equilibrium order-parameter amplitude.

        The normalized equilibrium solution is:

            |psi| = sqrt(1 - T/Tc)

        for:

            T < Tc

        and:

            |psi| = 0

        for:

            T >= Tc.
        """

        amplitude_squared = (
            self.equilibrium_amplitude_squared(
                reduced_temperature
            )
        )

        return np.sqrt(
            amplitude_squared
        )

    def equilibrium_order_parameter(
        self,
        reduced_temperature=None,
        phase=0.0,
    ):
        """
        Construct a uniform equilibrium order parameter.

        Parameters
        ----------
        reduced_temperature:
            T/Tc.

            If omitted, the value stored in the model
            parameters is used.

        phase:
            Uniform superconducting phase in radians.

        Returns
        -------
        complex
            Equilibrium order parameter.
        """

        if reduced_temperature is None:
            reduced_temperature = (
                self.parameters.reduced_temperature
            )

        amplitude = self.equilibrium_amplitude(
            reduced_temperature
        )

        return (
            amplitude *
            np.exp(1j * phase)
        )
    @staticmethod
    def supercurrent_density(
        psi,
        vector_potential_x,
        vector_potential_y,
        dx,
        dy,
        active_mask=None,
    ):
        """
        Calculate the dimensionless superconducting current density.

        Using the TDGL convention

            D = ∇ - iA

        the dimensionless supercurrent is

            j_s = Im(psi* D psi)

        with components

            j_s,x = Im(psi* Dx psi)

            j_s,y = Im(psi* Dy psi)

        Parameters
        ----------
        psi : ndarray
            Complex superconducting order parameter.

        vector_potential_x : ndarray
            x-component of the dimensionless vector potential.

        vector_potential_y : ndarray
            y-component of the dimensionless vector potential.

        dx : float
            Dimensionless grid spacing in x.

        dy : float
            Dimensionless grid spacing in y.

        Returns
        -------
        jx, jy : ndarray
            Dimensionless superconducting current-density
            components.
        """

        Dx_psi, Dy_psi = gauge_covariant_gradient(
            psi,
            vector_potential_x,
            vector_potential_y,
            dx,
            dy,
            active_mask,
        )

        jx = np.imag(
            np.conjugate(psi) * Dx_psi
        )

        jy = np.imag(
            np.conjugate(psi) * Dy_psi
        )

        return jx, jy
    def characteristic_time(
        self,
        critical_temperature,
        normal_conductivity=None,
        penetration_depth=None,
    ):
        """
        Return the reference physical TDGL characteristic time.

        This is the fixed time scale used to convert physical
        time into the normalized TDGL time used by the solver.
        """


        critical_temperature = float(
            critical_temperature
        )

        if critical_temperature <= 0.0:
            raise ValueError(
                "Critical temperature must be positive."
            )

        if self.parameters.normalization == "pytdgl":
            if normal_conductivity is None or penetration_depth is None:
                raise ValueError(
                    "pytdgl normalization requires normal conductivity and "
                    "penetration depth."
                )
            if normal_conductivity <= 0 or penetration_depth <= 0:
                raise ValueError("Conductivity and penetration depth must be positive.")
            return VACUUM_PERMEABILITY * normal_conductivity * penetration_depth**2
        return np.pi * REDUCED_PLANCK_CONSTANT / (
            8.0 * BOLTZMANN_CONSTANT * critical_temperature
        )
    def dimensional_to_normalized_time(
        self,
        dt,
        critical_temperature,
        normal_conductivity=None,
        penetration_depth=None,
    ):
        """
        Convert a physical timestep into normalized TDGL time.
        """

        if dt <= 0.0:
            raise ValueError(
                "Physical timestep must be positive."
            )

        tau_GL = self.characteristic_time(
            critical_temperature,
            normal_conductivity,
            penetration_depth,
        )

        return dt / tau_GL
