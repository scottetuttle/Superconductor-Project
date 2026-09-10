"""
Simulation configuration loader.

Defines the configuration objects used to construct an SHS
simulation from a JSON configuration file.

The configuration is divided into:

    - physical simulation parameters
    - thermal model parameters
    - electrical parameters
    - TDGL parameters
    - boundary conditions

Coupling-controller parameters are intentionally kept
separate from this physical simulation configuration.
"""

import json

from dataclasses import dataclass

from pathlib import Path


@dataclass
class ThermalConfig:
    """
    Thermal model configuration.
    """

    bath_temperature: float

    thermal_relaxation_rate: float = 0.0

    max_substep: float = 1e-13

    def validate(self):
        """
        Validate thermal configuration.
        """

        if self.bath_temperature < 0.0:
            raise ValueError(
                "Bath temperature must be non-negative."
            )

        if self.thermal_relaxation_rate < 0.0:
            raise ValueError(
                "Thermal relaxation rate must be non-negative."
            )

        if self.max_substep <= 0.0:
            raise ValueError(
                "Thermal max_substep must be positive."
            )

        return True


@dataclass
class ElectricalSolverConfig:
    """
    Numerical configuration for the electrical solver.
    """

    tolerance: float = 1e-12

    max_iterations: int = 10000

    omega: float = 1.7

    def validate(self):
        """
        Validate electrical solver controls.
        """

        if self.tolerance <= 0.0:
            raise ValueError(
                "Electrical solver tolerance must be positive."
            )

        if self.max_iterations <= 0:
            raise ValueError(
                "Electrical solver maximum iterations must be positive."
            )

        if not 0.0 < self.omega < 2.0:
            raise ValueError(
                "Electrical solver omega must satisfy 0 < omega < 2."
            )

        return True


@dataclass
class ElectricalConfig:
    """
    Electrical transport configuration.

    voltage_left and voltage_right are physical simulation
    parameters.

    solver contains numerical parameters used to solve
    the electrical potential equation.
    """

    voltage_left: float = 1e-3

    voltage_right: float = 0.0

    solver: ElectricalSolverConfig = None

    def __post_init__(self):
        if self.solver is None:
            self.solver = ElectricalSolverConfig()

    def validate(self):
        """
        Validate electrical configuration.
        """

        self.solver.validate()

        return True


@dataclass
class TDGLConfig:
    """
    TDGL model configuration.

    Parameters
    ----------
    u:
        TDGL relaxation parameter.

    gamma:
        Amplitude/phase coupling parameter.

    kappa:
        Ginzburg-Landau parameter.

    max_normalized_timestep:
        Maximum normalized timestep used by the explicit
        TDGL integrator.
    """

    u: float = 5.79

    gamma: float = 0.0

    kappa: float = 1.0

    max_normalized_timestep: float = 0.01

    def validate(self):
        """
        Validate TDGL configuration.
        """

        if self.u <= 0.0:
            raise ValueError(
                "TDGL parameter u must be positive."
            )

        if self.kappa <= 0.0:
            raise ValueError(
                "TDGL parameter kappa must be positive."
            )

        if self.max_normalized_timestep <= 0.0:
            raise ValueError(
                "TDGL max_normalized_timestep must be positive."
            )

        return True


@dataclass
class SimulationConfig:
    """
    Complete physical simulation configuration.
    """

    name: str

    geometry: str

    material: str

    temperature: float

    current: float

    duration: float

    dt: float

    boundaries: dict

    tdgl_boundaries: dict

    thermal: ThermalConfig

    electrical: ElectricalConfig

    tdgl: TDGLConfig

    def validate(self):
        """
        Validate the complete simulation configuration.
        """

        if self.temperature < 0.0:
            raise ValueError(
                "Simulation temperature must be non-negative."
            )

        if self.duration <= 0.0:
            raise ValueError(
                "Simulation duration must be positive."
            )

        if self.dt <= 0.0:
            raise ValueError(
                "Simulation timestep must be positive."
            )

        self.thermal.validate()

        self.electrical.validate()

        self.tdgl.validate()

        return True


def load_simulation(filepath):
    """
    Load a simulation configuration from JSON.

    Parameters
    ----------
    filepath:
        Path to the simulation configuration.

    Returns
    -------
    SimulationConfig
        Parsed and validated simulation configuration.
    """

    filepath = Path(filepath)

    with open(
        filepath,
        "r",
        encoding="utf-8",
    ) as file:
        data = json.load(file)

    #
    # Thermal configuration.
    #

    thermal_data = data.get(
        "thermal",
        {},
    )

    thermal = ThermalConfig(
        bath_temperature=float(
            thermal_data.get(
                "bath_temperature",
                data["temperature"],
            )
        ),

        thermal_relaxation_rate=float(
            thermal_data.get(
                "thermal_relaxation_rate",
                0.0,
            )
        ),

        max_substep=float(
            thermal_data.get(
                "max_substep",
                1e-13,
            )
        ),
    )

    #
    # Electrical configuration.
    #

    electrical_data = data.get(
        "electrical",
        {},
    )

    electrical_solver_data = electrical_data.get(
        "solver",
        {},
    )

    electrical_solver = ElectricalSolverConfig(
        tolerance=float(
            electrical_solver_data.get(
                "tolerance",
                1e-10,
            )
        ),

        max_iterations=int(
            electrical_solver_data.get(
                "max_iterations",
                10000,
            )
        ),

        omega=float(
            electrical_solver_data.get(
                "omega",
                1.7,
            )
        ),
    )

    electrical = ElectricalConfig(
        voltage_left=float(
            electrical_data.get(
                "voltage_left",
                1e-3,
            )
        ),

        voltage_right=float(
            electrical_data.get(
                "voltage_right",
                0.0,
            )
        ),

        solver=electrical_solver,
    )

    #
    # TDGL configuration.
    #

    tdgl_data = data.get(
        "tdgl",
        {},
    )

    tdgl = TDGLConfig(
        u=float(
            tdgl_data.get(
                "u",
                5.79,
            )
        ),

        gamma=float(
            tdgl_data.get(
                "gamma",
                0.0,
            )
        ),

        kappa=float(
            tdgl_data.get(
                "kappa",
                1.0,
            )
        ),

        max_normalized_timestep=float(
            tdgl_data.get(
                "max_normalized_timestep",
                0.01,
            )
        ),
    )

    #
    # Boundary configurations.
    #

    boundaries = data.get(
        "boundaries",
        {}
    )

    tdgl_boundaries = data.get(
        "tdgl_boundaries",
        {}
    )

    #
    # Complete configuration.
    #

    config = SimulationConfig(
        name=data["name"],

        geometry=data["geometry"],

        material=data["material"],

        temperature=float(
            data["temperature"]
        ),

        current=float(
            data["current"]["value"]
        ),

        duration=float(
            data["time"]["duration"]
        ),

        dt=float(
            data["time"]["dt"]
        ),

        boundaries=boundaries,

        tdgl_boundaries=tdgl_boundaries,

        thermal=thermal,

        electrical=electrical,

        tdgl=tdgl,
    )

    config.validate()

    return config