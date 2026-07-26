"""
SHS Thermal Physics Module
==========================

Purpose:
--------
Defines the thermal physics models required to simulate temperature
evolution in superconducting devices.

The thermal module is responsible for modeling the interaction between
localized heating sources and thermal transport throughout the device.

Primary physical processes:
--------------------------

1. Heat diffusion

    Temperature evolves according to thermal conduction:

        C(T) dT/dt = ∇ · (k(T) ∇T) + Q

    where:

        C(T)  = heat capacity
        k(T)  = thermal conductivity
        Q     = volumetric heat sources


2. Joule heating

Electrical current can locally heat regions that become resistive:

        Q_J = J²ρ(T,J)

where:

        J     = current density
        ρ     = electrical resistivity


3. Optical heating

Laser excitation is modeled as a spatially and temporally varying
heat source:

        Q_opt(x,y,t)

Supported profiles will include:

- Gaussian beams
- Pulsed excitation
- Continuous-wave lasers
- Moving hotspots
- Arbitrary user-defined profiles


4. Thermal relaxation

The module may include coupling to a substrate or bath:

        Q_loss = -G(T - T_bath)

where:

        G       = thermal boundary conductance
        T_bath  = substrate temperature


Future capabilities:
--------------------

- Electron-phonon coupling
- Two-temperature models
- Thermal boundary resistance
- Nonlinear material properties
- Temperature-dependent material parameters


Dependencies:
-------------

Requires:

- Material properties from shs.materials
- Geometry and mesh from shs.geometry


Current implementation status:
------------------------------

Architecture only.

No numerical solver implemented yet.
"""




from dataclasses import dataclass


@dataclass
class ThermalModel:
    """
    Defines thermal properties of a superconducting device.

    Units:
    - thermal_conductivity: W/(m*K)
    - heat_capacity: J/(m^3*K)
    - bath_temperature: K
    - thermal_boundary_conductance: W/(m^2*K)
    """

    thermal_conductivity: float

    heat_capacity: float

    bath_temperature: float

    thermal_boundary_conductance: float | None = None


"""
Thermal physics model for SHS.
"""

from dataclasses import dataclass


@dataclass
class ThermalModel:

    thermal_conductivity: float

    heat_capacity: float

    bath_temperature: float

    thermal_relaxation_rate: float | None = None


    """
    thermal_relaxation_rate

    Coupling strength between the superconducting film
    and the thermal bath (substrate).

    Larger values remove heat more quickly.

    Units:
        1/s
    """
    thermal_relaxation_rate: float


    def thermal_diffusivity(self):
        """
        Calculate thermal diffusivity.

        alpha = k / C

        Units:
        m^2/s
        """

        return (
            self.thermal_conductivity /
            self.heat_capacity
        )