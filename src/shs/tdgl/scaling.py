"""
TDGL physical and dimensionless scaling.

This module defines the conversion between physical SI quantities
used by SHS and the dimensionless quantities used internally by
the normalized TDGL equations.

The current normalized TDGL formulation uses the coherence length
as the characteristic spatial scale:

    x' = x / xi

and the critical temperature as the characteristic temperature scale:

    t = T / Tc

The electromagnetic normalization follows the standard GL scales:

    A0 = Phi0 / (2 pi xi)

    B0 = Phi0 / (2 pi xi^2)

    J0 = 4 xi Bc2 / (mu0 lambda^2)

For the default pyTDGL convention the caller supplies
tau0=mu0*sigma*lambda^2. The legacy convention remains selectable for old
comparison cases.

This separation is intentional:

    physical SHS quantities
            |
            v
        TDGLScales
            |
            v
    dimensionless TDGL quantities
"""


from dataclasses import dataclass

import numpy as np
from shs.utils.constants import SUPERCONDUCTING_FLUX_QUANTUM, VACUUM_PERMEABILITY


@dataclass(frozen=True)
class TDGLScales:
    """
    Characteristic physical scales used by the normalized TDGL model.

    Parameters
    ----------
    xi:
        Characteristic Ginzburg-Landau coherence length [m].

    lambda_:
        London penetration depth [m].

    critical_temperature:
        Critical temperature Tc [K].

    time_scale:
        Characteristic TDGL time [s].

        This must be supplied explicitly because the precise TDGL
        time normalization depends on the chosen physical TDGL
        formulation.

    flux_quantum:
        Magnetic flux quantum Phi0 [Wb].

        Defaults to the superconducting flux quantum.
    """

    xi: float
    lambda_: float
    critical_temperature: float
    time_scale: float

    normalization: str = "legacy_gl"

    flux_quantum: float = SUPERCONDUCTING_FLUX_QUANTUM

    def __post_init__(self):
        """
        Validate the physical scales.
        """

        if self.xi <= 0.0:
            raise ValueError(
                "Coherence length xi must be positive."
            )

        if self.lambda_ <= 0.0:
            raise ValueError(
                "Penetration depth lambda must be positive."
            )

        if self.critical_temperature <= 0.0:
            raise ValueError(
                "Critical temperature must be positive."
            )

        if self.time_scale <= 0.0:
            raise ValueError(
                "TDGL time scale must be positive."
            )

        if self.flux_quantum <= 0.0:
            raise ValueError(
                "Flux quantum must be positive."
            )
        if self.normalization not in {"pytdgl", "legacy_gl"}:
            raise ValueError("Unknown TDGL normalization.")

    # ------------------------------------------------------------------
    # Characteristic electromagnetic scales
    # ------------------------------------------------------------------

    @property
    def vector_potential_scale(self):
        """
        Characteristic vector-potential scale.

        A0 = Phi0 / (2 pi xi)

        Units:
            T m
        """

        return (
            self.flux_quantum /
            (2.0 * np.pi * self.xi)
        )

    @property
    def magnetic_field_scale(self):
        """
        Characteristic magnetic-field scale.

        B0 = Phi0 / (2 pi xi^2)

        Units:
            T
        """

        return (
            self.flux_quantum /
            (
                2.0 *
                np.pi *
                self.xi**2
            )
        )
    @property
    def current_density_scale(self):
        """
        Characteristic GL current-density scale.

        pyTDGL: J0 = 4 xi Bc2 / (mu0 lambda^2)
        legacy: J0 = Phi0 / (2 pi mu0 lambda^2 xi)

        Units:
            A/m^2
        """

        base = (
            self.flux_quantum /
            (
                2.0 *
                np.pi *
                VACUUM_PERMEABILITY *
                self.lambda_**2 *
                self.xi
            )
        )
        return 4.0 * base if self.normalization == "pytdgl" else base

    @property
    def electric_field_scale(self):
        """
        Characteristic electric-field scale.

        E0 = A0 / t0

        Units:
            V/m
        """

        return self.scalar_potential_scale / self.xi

    @property
    def scalar_potential_scale(self):
        """Characteristic electric scalar-potential scale [V].

        phi0 = xi E0 = Phi0 / (2 pi t0)
        """
        base = self.flux_quantum / (2.0 * np.pi * self.time_scale)
        return 4.0 * base if self.normalization == "pytdgl" else base

    # ------------------------------------------------------------------
    # Length
    # ------------------------------------------------------------------

    def length_to_dimensionless(
        self,
        length,
    ):
        """
        Convert physical length [m] to dimensionless length.

        x' = x / xi
        """

        return np.asarray(length) / self.xi

    def length_to_physical(
        self,
        dimensionless_length,
    ):
        """
        Convert dimensionless length to physical length [m].
        """

        return (
            np.asarray(dimensionless_length) *
            self.xi
        )

    # ------------------------------------------------------------------
    # Temperature
    # ------------------------------------------------------------------

    def temperature_to_dimensionless(
        self,
        temperature,
    ):
        """
        Convert physical temperature [K] to reduced temperature.

        t = T / Tc
        """

        return (
            np.asarray(temperature) /
            self.critical_temperature
        )

    def temperature_to_physical(
        self,
        reduced_temperature,
    ):
        """
        Convert reduced temperature to physical temperature [K].
        """

        return (
            np.asarray(reduced_temperature) *
            self.critical_temperature
        )

    # ------------------------------------------------------------------
    # Time
    # ------------------------------------------------------------------

    def time_to_dimensionless(
        self,
        time,
    ):
        """
        Convert physical time [s] to dimensionless TDGL time.
        """

        return (
            np.asarray(time) /
            self.time_scale
        )

    def time_to_physical(
        self,
        dimensionless_time,
    ):
        """
        Convert dimensionless TDGL time to physical time [s].
        """

        return (
            np.asarray(dimensionless_time) *
            self.time_scale
        )

    # ------------------------------------------------------------------
    # Vector potential
    # ------------------------------------------------------------------

    def vector_potential_to_dimensionless(
        self,
        vector_potential,
    ):
        """
        Convert physical vector potential [T m] to
        dimensionless vector potential.
        """

        return (
            np.asarray(vector_potential) /
            self.vector_potential_scale
        )

    def vector_potential_to_physical(
        self,
        dimensionless_vector_potential,
    ):
        """
        Convert dimensionless vector potential to [T m].
        """

        return (
            np.asarray(dimensionless_vector_potential) *
            self.vector_potential_scale
        )

    # ------------------------------------------------------------------
    # Magnetic field
    # ------------------------------------------------------------------

    def magnetic_field_to_dimensionless(
        self,
        magnetic_field,
    ):
        """
        Convert physical magnetic field [T] to dimensionless
        magnetic field.
        """

        return (
            np.asarray(magnetic_field) /
            self.magnetic_field_scale
        )

    def magnetic_field_to_physical(
        self,
        dimensionless_magnetic_field,
    ):
        """
        Convert dimensionless magnetic field to [T].
        """

        return (
            np.asarray(dimensionless_magnetic_field) *
            self.magnetic_field_scale
        )

    # ------------------------------------------------------------------
    # Current density
    # ------------------------------------------------------------------

    def current_density_to_dimensionless(
        self,
        current_density,
    ):
        """
        Convert physical current density [A/m^2] to
        dimensionless current density.
        """

        return (
            np.asarray(current_density) /
            self.current_density_scale
        )

    def current_density_to_physical(
        self,
        dimensionless_current_density,
    ):
        """
        Convert dimensionless current density to [A/m^2].
        """

        return (
            np.asarray(dimensionless_current_density) *
            self.current_density_scale
        )

    def supercurrent_density_to_physical(self, covariant_current):
        """Convert ``Im(psi* D psi)`` to physical supercurrent density.

        In the pyTDGL units, ``J0`` contains a factor of four relative to the
        GL covariant-current scale. The legacy convention uses the base scale
        directly.
        """
        factor = 0.25 if self.normalization == "pytdgl" else 1.0
        return np.asarray(covariant_current) * self.current_density_scale * factor

    # ------------------------------------------------------------------
    # Electric field
    # ------------------------------------------------------------------

    def electric_field_to_dimensionless(
        self,
        electric_field,
    ):
        """
        Convert physical electric field [V/m] to
        dimensionless electric field.
        """

        return (
            np.asarray(electric_field) /
            self.electric_field_scale
        )

    def electric_field_to_physical(
        self,
        dimensionless_electric_field,
    ):
        """
        Convert dimensionless electric field to [V/m].
        """

        return (
            np.asarray(dimensionless_electric_field) *
            self.electric_field_scale
        )

    def scalar_potential_to_dimensionless(self, scalar_potential):
        """Convert physical electric potential [V] to TDGL units."""
        return np.asarray(scalar_potential) / self.scalar_potential_scale

    def scalar_potential_to_physical(self, dimensionless_potential):
        """Convert dimensionless TDGL scalar potential to volts."""
        return np.asarray(dimensionless_potential) * self.scalar_potential_scale
