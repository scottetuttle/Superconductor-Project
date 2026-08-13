"""
TDGL physical model.

Defines the normalized Time-Dependent Ginzburg-Landau
physics used by SHS.

The current dimensionless TDGL equation is:

    u dpsi/dt =
        D^2 psi
        + (1 - T/Tc) psi
        - |psi|^2 psi

where:

    D = nabla - i A

The numerical solver is responsible for advancing
the solution in time.
"""

from dataclasses import dataclass
import numpy as np

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

    @staticmethod
    def equilibrium_amplitude_squared(
        reduced_temperature,
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

        result = np.maximum(
            1.0 - t,
            0.0
        )

        if result.ndim == 0:
            return float(result)

        return result

    @staticmethod
    def equilibrium_amplitude(
        reduced_temperature,
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
            TDGLModel.equilibrium_amplitude_squared(
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
        )

        jx = np.imag(
            np.conjugate(psi) * Dx_psi
        )

        jy = np.imag(
            np.conjugate(psi) * Dy_psi
        )

        return jx, jy

