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

Coupling-controller parameters occupy a separate numerical subsection.
"""

import json

from dataclasses import dataclass, field

import numpy as np

from pathlib import Path

from .defaults import default_section


_THERMAL_DEFAULTS = default_section("thermal")
_LASER_DEFAULTS = default_section("laser")
_PINNING_DEFAULTS = default_section("pinning")
_JOSEPHSON_DEFAULTS = default_section("josephson")
_ELECTRICAL_DEFAULTS = default_section("electrical")
_TDGL_DEFAULTS = default_section("tdgl")
_ELECTROMAGNETIC_DEFAULTS = default_section("electromagnetic")
_COUPLING_DEFAULTS = default_section("coupling")
_NUMERICS_DEFAULTS = default_section("numerics")


@dataclass
class ThermalConfig:
    """
    Thermal model configuration.
    """

    bath_temperature: float

    model: str = _THERMAL_DEFAULTS["model"]
    electron_heat_capacity_fraction: float = _THERMAL_DEFAULTS["electron_heat_capacity_fraction"]
    electron_thermal_conductivity_fraction: float = _THERMAL_DEFAULTS["electron_thermal_conductivity_fraction"]
    electron_phonon_coupling_W_m3_K: float = _THERMAL_DEFAULTS["electron_phonon_coupling_W_m3_K"]
    phonon_escape_rate_per_s: float = _THERMAL_DEFAULTS["phonon_escape_rate_per_s"]

    thermal_relaxation_rate: float = _THERMAL_DEFAULTS["thermal_relaxation_rate"]

    max_substep: float = _THERMAL_DEFAULTS["max_substep"]

    stability_safety_factor: float = _THERMAL_DEFAULTS["stability_safety_factor"]

    def validate(self):
        """
        Validate thermal configuration.
        """

        if self.model not in {"single_temperature", "two_temperature"}:
            raise ValueError("Unknown thermal model.")
        if not 0 < self.electron_heat_capacity_fraction < 1:
            raise ValueError("Electron heat-capacity fraction must lie in (0, 1).")
        if not 0 <= self.electron_thermal_conductivity_fraction <= 1:
            raise ValueError("Electron thermal-conductivity fraction must lie in [0, 1].")
        if (not np.isfinite(self.electron_phonon_coupling_W_m3_K)
                or self.electron_phonon_coupling_W_m3_K < 0
                or not np.isfinite(self.phonon_escape_rate_per_s)
                or self.phonon_escape_rate_per_s < 0):
            raise ValueError("Electron-phonon coupling and phonon escape must be finite and nonnegative.")
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
        if self.model == "two_temperature" and self.thermal_relaxation_rate:
            raise ValueError("Use phonon_escape_rate_per_s for two-temperature bath cooling.")

        if not 0.0 < self.stability_safety_factor <= 1.0:
            raise ValueError("Thermal stability_safety_factor must lie in (0, 1].")

        return True


@dataclass
class ElectricalSolverConfig:
    """
    Numerical configuration for the electrical solver.
    """

    backend: str = _ELECTRICAL_DEFAULTS["solver"]["backend"]

    tolerance: float = _ELECTRICAL_DEFAULTS["solver"]["tolerance"]

    max_iterations: int = _ELECTRICAL_DEFAULTS["solver"]["max_iterations"]

    omega: float = _ELECTRICAL_DEFAULTS["solver"]["omega"]

    residual_check_interval: int = _ELECTRICAL_DEFAULTS["solver"]["residual_check_interval"]

    def validate(self):
        """
        Validate electrical solver controls.
        """

        if self.backend not in {"red_black_sor", "sparse_direct"}:
            raise ValueError("Unknown electrical solver backend.")

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

        if self.residual_check_interval < 1:
            raise ValueError("Electrical residual_check_interval must be positive.")

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

    drive_mode: str = _ELECTRICAL_DEFAULTS["drive_mode"]

    source_contact: str = _ELECTRICAL_DEFAULTS["source_contact"]

    sink_contact: str = _ELECTRICAL_DEFAULTS["sink_contact"]

    reference_voltage: float = _ELECTRICAL_DEFAULTS["reference_voltage"]

    voltage_left: float = _ELECTRICAL_DEFAULTS["voltage_left"]

    voltage_right: float = _ELECTRICAL_DEFAULTS["voltage_right"]

    normal_conductivity_model: str = _ELECTRICAL_DEFAULTS["normal_conductivity_model"]

    solver: ElectricalSolverConfig = None

    def __post_init__(self):
        if self.solver is None:
            self.solver = ElectricalSolverConfig()

    def validate(self):
        """
        Validate electrical configuration.
        """

        self.solver.validate()

        if self.drive_mode not in {"voltage", "current"}:
            raise ValueError("Electrical drive_mode must be voltage or current.")
        if not self.source_contact or not self.sink_contact:
            raise ValueError("Electrical source and sink contact names are required.")
        if self.source_contact == self.sink_contact:
            raise ValueError("Electrical source and sink contacts must differ.")
        if not np.isfinite(self.reference_voltage):
            raise ValueError("Electrical reference_voltage must be finite.")

        if self.normal_conductivity_model not in {"constant", "condensate_depletion"}:
            raise ValueError("Unknown normal conductivity model.")

        return True


@dataclass(frozen=True)
class LaserWaypointConfig:
    time_s: float
    x_m: float
    y_m: float
    power_fraction: float = 1.0


@dataclass
class LaserConfig:
    """Absorbed Gaussian laser power and its piecewise-linear path."""

    enabled: bool = _LASER_DEFAULTS["enabled"]
    absorbed_power_W: float = _LASER_DEFAULTS["absorbed_power_W"]
    sigma_m: float = _LASER_DEFAULTS["sigma_m"]
    repeat_path: bool = _LASER_DEFAULTS["repeat_path"]
    waypoints: list[LaserWaypointConfig] = field(default_factory=list)

    def validate(self):
        if not isinstance(self.enabled, bool) or not isinstance(self.repeat_path, bool):
            raise ValueError("Laser enabled and repeat_path must be boolean.")
        if not np.isfinite(self.absorbed_power_W) or self.absorbed_power_W < 0:
            raise ValueError("Laser absorbed power must be finite and nonnegative.")
        if not np.isfinite(self.sigma_m) or self.sigma_m <= 0:
            raise ValueError("Laser sigma must be positive and finite.")
        if self.enabled and not self.waypoints:
            raise ValueError("An enabled laser requires at least one waypoint.")
        times = [point.time_s for point in self.waypoints]
        if any(not np.isfinite(value) for point in self.waypoints
               for value in (point.time_s, point.x_m, point.y_m, point.power_fraction)):
            raise ValueError("Laser waypoint values must be finite.")
        if any(not 0.0 <= point.power_fraction <= 1.0 for point in self.waypoints):
            raise ValueError("Laser waypoint power_fraction must lie in [0, 1].")
        if any(time < 0 for time in times) or any(b <= a for a, b in zip(times, times[1:])):
            raise ValueError("Laser waypoint times must be nonnegative and strictly increasing.")


@dataclass(frozen=True)
class PinningSiteConfig:
    """Gaussian material defect that locally reduces the linear GL coefficient."""

    x_m: float
    y_m: float
    sigma_m: float
    strength: float


@dataclass
class PinningConfig:
    """Static phenomenological material-pinning landscape."""

    enabled: bool = _PINNING_DEFAULTS["enabled"]
    maximum_suppression: float = _PINNING_DEFAULTS["maximum_suppression"]
    sites: list[PinningSiteConfig] = field(default_factory=list)

    def validate(self):
        if not isinstance(self.enabled, bool):
            raise ValueError("Pinning enabled must be boolean.")
        if not np.isfinite(self.maximum_suppression) or self.maximum_suppression <= 0:
            raise ValueError("Pinning maximum_suppression must be positive and finite.")
        for site in self.sites:
            if not all(np.isfinite(value) for value in (
                    site.x_m, site.y_m, site.sigma_m, site.strength)):
                raise ValueError("Pinning-site values must be finite.")
            if site.sigma_m <= 0:
                raise ValueError("Pinning-site sigma_m must be positive.")
            if site.strength < 0:
                raise ValueError("Pinning-site strength must be nonnegative.")


@dataclass
class JosephsonConfig:
    """Resolved vertical weak link represented by local GL suppression."""
    enabled: bool = _JOSEPHSON_DEFAULTS["enabled"]
    center_x_m: float = _JOSEPHSON_DEFAULTS["center_x_m"]
    width_m: float = _JOSEPHSON_DEFAULTS["width_m"]
    coefficient_suppression: float = _JOSEPHSON_DEFAULTS["coefficient_suppression"]

    def validate(self):
        if not isinstance(self.enabled, bool):
            raise ValueError("Josephson enabled must be boolean.")
        if not all(np.isfinite(v) for v in (
                self.center_x_m, self.width_m, self.coefficient_suppression)):
            raise ValueError("Josephson weak-link values must be finite.")
        if self.width_m <= 0:
            raise ValueError("Josephson weak-link width_m must be positive.")
        if self.coefficient_suppression < 0:
            raise ValueError("Josephson coefficient_suppression must be nonnegative.")


@dataclass
class ElectromagneticConfig:
    """Magnetic-field coupling and self-consistent screening controls."""

    include_self_field: bool = _ELECTROMAGNETIC_DEFAULTS["include_self_field"]
    include_displacement_current: bool = _ELECTROMAGNETIC_DEFAULTS["include_displacement_current"]
    screening_tolerance: float = _ELECTROMAGNETIC_DEFAULTS["screening_tolerance"]
    screening_max_iterations: int = _ELECTROMAGNETIC_DEFAULTS["screening_max_iterations"]
    screening_step_size: float = _ELECTROMAGNETIC_DEFAULTS["screening_step_size"]
    screening_source_stride: int = _ELECTROMAGNETIC_DEFAULTS["screening_source_stride"]

    def validate(self):
        if not isinstance(self.include_self_field, bool):
            raise ValueError("include_self_field must be boolean.")
        if not isinstance(self.include_displacement_current, bool):
            raise ValueError("include_displacement_current must be boolean.")
        if not np.isfinite(self.screening_tolerance) or self.screening_tolerance <= 0:
            raise ValueError("screening_tolerance must be positive and finite.")
        if not isinstance(self.screening_max_iterations, int) or self.screening_max_iterations < 1:
            raise ValueError("screening_max_iterations must be a positive integer.")
        if not 0 < self.screening_step_size <= 1:
            raise ValueError("screening_step_size must lie in (0, 1].")
        if not isinstance(self.screening_source_stride, int) or self.screening_source_stride < 1:
            raise ValueError("screening_source_stride must be a positive integer.")


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

    u: float = _TDGL_DEFAULTS["u"]
    time_integrator: str = _TDGL_DEFAULTS["time_integrator"]

    gamma: float = _TDGL_DEFAULTS["gamma"]

    kappa: float = _TDGL_DEFAULTS["kappa"]

    normalization: str = _TDGL_DEFAULTS["normalization"]

    temperature_model: str = _TDGL_DEFAULTS["temperature_model"]

    include_scalar_potential: bool = _TDGL_DEFAULTS["include_scalar_potential"]

    max_normalized_timestep: float = _TDGL_DEFAULTS["max_normalized_timestep"]

    stability_safety_factor: float = _TDGL_DEFAULTS["stability_safety_factor"]

    def validate(self):
        """
        Validate TDGL configuration.
        """

        if self.time_integrator not in {"euler", "heun"}:
            raise ValueError("TDGL time_integrator must be euler or heun.")
        if self.u <= 0.0:
            raise ValueError(
                "TDGL parameter u must be positive."
            )

        if self.kappa <= 0.0:
            raise ValueError(
                "TDGL parameter kappa must be positive."
            )

        if self.normalization not in {"pytdgl", "legacy_gl"}:
            raise ValueError("TDGL normalization must be pytdgl or legacy_gl.")

        if self.temperature_model not in {
            "tc_over_t_minus_one", "one_minus_t_over_tc"
        }:
            raise ValueError("Unknown TDGL temperature coefficient model.")

        if self.max_normalized_timestep <= 0.0:
            raise ValueError(
                "TDGL max_normalized_timestep must be positive."
            )

        if not isinstance(self.include_scalar_potential, bool):
            raise ValueError("TDGL include_scalar_potential must be boolean.")

        if not 0.0 < self.stability_safety_factor <= 1.0:
            raise ValueError("TDGL stability_safety_factor must lie in (0, 1].")

        return True


@dataclass
class CouplingConfig:
    """Numerical coupling controls; field scales use the shared fields' SI units."""

    tolerance: float = _COUPLING_DEFAULTS["tolerance"]
    max_iterations: int = _COUPLING_DEFAULTS["max_iterations"]
    prediction_window: int = _COUPLING_DEFAULTS["prediction_window"]
    reasonable_iterations: int = _COUPLING_DEFAULTS["reasonable_iterations"]
    check_interval: int = _COUPLING_DEFAULTS["check_interval"]
    minimum_iterations: int = _COUPLING_DEFAULTS["minimum_iterations"]
    fast_ratio: float = _COUPLING_DEFAULTS["fast_ratio"]
    healthy_ratio: float = _COUPLING_DEFAULTS["healthy_ratio"]
    stall_ratio: float = _COUPLING_DEFAULTS["stall_ratio"]
    oscillation_window: int = _COUPLING_DEFAULTS["oscillation_window"]
    electrical_tolerance_fraction: float = _COUPLING_DEFAULTS["electrical_tolerance_fraction"]
    field_scales: dict = field(default_factory=lambda: default_section("coupling")["field_scales"])

    def validate(self):
        if not np.isfinite(self.tolerance) or self.tolerance <= 0:
            raise ValueError('Coupling tolerance must be positive and finite.')
        for name in ('max_iterations', 'prediction_window', 'reasonable_iterations',
                     'check_interval', 'minimum_iterations'):
            value = getattr(self, name)
            if not isinstance(value, int) or value < 1:
                raise ValueError(f'{name} must be a positive integer.')
        if self.prediction_window < 2 or self.minimum_iterations > self.max_iterations:
            raise ValueError('Invalid coupling iteration limits.')
        if self.oscillation_window < 4:
            raise ValueError('Coupling oscillation_window must be at least 4.')
        if not 0.0 < self.electrical_tolerance_fraction <= 1.0:
            raise ValueError('electrical_tolerance_fraction must lie in (0, 1].')
        if not 0.0 < self.fast_ratio < self.healthy_ratio < self.stall_ratio <= 1.0:
            raise ValueError('Coupling residual ratios must increase from fast to stall.')
        expected = {'temperature', 'psi', 'voltage', 'vector_potential_x', 'vector_potential_y'}
        if set(self.field_scales) != expected or any(
                not np.isfinite(v) or v <= 0 for v in self.field_scales.values()):
            raise ValueError('Supply positive finite scales for all primary fields.')


@dataclass
class NumericsConfig:
    """Shared numerical safeguards used by multiple solvers."""

    max_internal_substeps: int = _NUMERICS_DEFAULTS["max_internal_substeps"]
    contact_roundoff_ulps: int = _NUMERICS_DEFAULTS["contact_roundoff_ulps"]

    def validate(self):
        if self.max_internal_substeps < 1:
            raise ValueError('max_internal_substeps must be positive.')
        if self.contact_roundoff_ulps < 0:
            raise ValueError('contact_roundoff_ulps cannot be negative.')


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

    laser: LaserConfig

    pinning: PinningConfig

    josephson: JosephsonConfig

    electrical: ElectricalConfig

    electromagnetic: ElectromagneticConfig

    tdgl: TDGLConfig

    coupling: CouplingConfig = field(default_factory=CouplingConfig)

    numerics: NumericsConfig = field(default_factory=NumericsConfig)

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

        self.laser.validate()

        self.pinning.validate()

        self.josephson.validate()

        self.electrical.validate()

        self.electromagnetic.validate()

        self.tdgl.validate()

        self.coupling.validate()

        self.numerics.validate()

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
        model=thermal_data.get("model", _THERMAL_DEFAULTS["model"]),
        electron_heat_capacity_fraction=float(thermal_data.get("electron_heat_capacity_fraction", _THERMAL_DEFAULTS["electron_heat_capacity_fraction"])),
        electron_thermal_conductivity_fraction=float(thermal_data.get("electron_thermal_conductivity_fraction", _THERMAL_DEFAULTS["electron_thermal_conductivity_fraction"])),
        electron_phonon_coupling_W_m3_K=float(thermal_data.get("electron_phonon_coupling_W_m3_K", _THERMAL_DEFAULTS["electron_phonon_coupling_W_m3_K"])),
        phonon_escape_rate_per_s=float(thermal_data.get("phonon_escape_rate_per_s", _THERMAL_DEFAULTS["phonon_escape_rate_per_s"])),
        bath_temperature=float(
            thermal_data.get(
                "bath_temperature",
                data["temperature"],
            )
        ),

        thermal_relaxation_rate=float(
            thermal_data.get(
                "thermal_relaxation_rate",
                _THERMAL_DEFAULTS["thermal_relaxation_rate"],
            )
        ),

        max_substep=float(
            thermal_data.get(
                "max_substep",
                _THERMAL_DEFAULTS["max_substep"],
            )
        ),

        stability_safety_factor=float(thermal_data.get(
            "stability_safety_factor", _THERMAL_DEFAULTS["stability_safety_factor"])),
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
        backend=electrical_solver_data.get(
            "backend", _ELECTRICAL_DEFAULTS["solver"]["backend"]
        ),

        tolerance=float(
            electrical_solver_data.get(
                "tolerance",
                _ELECTRICAL_DEFAULTS["solver"]["tolerance"],
            )
        ),

        max_iterations=int(
            electrical_solver_data.get(
                "max_iterations",
                _ELECTRICAL_DEFAULTS["solver"]["max_iterations"],
            )
        ),

        omega=float(
            electrical_solver_data.get(
                "omega",
                _ELECTRICAL_DEFAULTS["solver"]["omega"],
            )
        ),

        residual_check_interval=int(electrical_solver_data.get(
            "residual_check_interval",
            _ELECTRICAL_DEFAULTS["solver"]["residual_check_interval"],
        )),
    )

    laser_data = data.get("laser", {})
    laser = LaserConfig(
        enabled=laser_data.get("enabled", _LASER_DEFAULTS["enabled"]),
        absorbed_power_W=float(laser_data.get(
            "absorbed_power_W", _LASER_DEFAULTS["absorbed_power_W"]
        )),
        sigma_m=float(laser_data.get("sigma_m", _LASER_DEFAULTS["sigma_m"])),
        repeat_path=laser_data.get("repeat_path", _LASER_DEFAULTS["repeat_path"]),
        waypoints=[LaserWaypointConfig(
            time_s=float(point["time_s"]), x_m=float(point["x_m"]),
            y_m=float(point["y_m"]),
            power_fraction=float(point.get("power_fraction", 1.0)),
        ) for point in laser_data.get("waypoints", _LASER_DEFAULTS["waypoints"])],
    )

    pinning_data = data.get("pinning", {})
    pinning = PinningConfig(
        enabled=pinning_data.get("enabled", _PINNING_DEFAULTS["enabled"]),
        maximum_suppression=float(pinning_data.get(
            "maximum_suppression", _PINNING_DEFAULTS["maximum_suppression"]
        )),
        sites=[PinningSiteConfig(
            x_m=float(site["x_m"]), y_m=float(site["y_m"]),
            sigma_m=float(site["sigma_m"]), strength=float(site["strength"]),
        ) for site in pinning_data.get("sites", _PINNING_DEFAULTS["sites"])],
    )

    josephson_data = data.get("josephson", {})
    josephson = JosephsonConfig(
        enabled=josephson_data.get("enabled", _JOSEPHSON_DEFAULTS["enabled"]),
        center_x_m=float(josephson_data.get("center_x_m", _JOSEPHSON_DEFAULTS["center_x_m"])),
        width_m=float(josephson_data.get("width_m", _JOSEPHSON_DEFAULTS["width_m"])),
        coefficient_suppression=float(josephson_data.get(
            "coefficient_suppression", _JOSEPHSON_DEFAULTS["coefficient_suppression"])),
    )

    electrical = ElectricalConfig(
        drive_mode=electrical_data.get(
            "drive_mode", _ELECTRICAL_DEFAULTS["drive_mode"]
        ),

        source_contact=electrical_data.get(
            "source_contact", _ELECTRICAL_DEFAULTS["source_contact"]
        ),

        sink_contact=electrical_data.get(
            "sink_contact", _ELECTRICAL_DEFAULTS["sink_contact"]
        ),

        reference_voltage=float(electrical_data.get(
            "reference_voltage", _ELECTRICAL_DEFAULTS["reference_voltage"]
        )),

        voltage_left=float(
            electrical_data.get(
                "voltage_left",
                _ELECTRICAL_DEFAULTS["voltage_left"],
            )
        ),

        voltage_right=float(
            electrical_data.get(
                "voltage_right",
                _ELECTRICAL_DEFAULTS["voltage_right"],
            )
        ),

        normal_conductivity_model=electrical_data.get(
            "normal_conductivity_model",
            _ELECTRICAL_DEFAULTS["normal_conductivity_model"],
        ),

        solver=electrical_solver,
    )

    electromagnetic_data = data.get("electromagnetic", {})
    electromagnetic = ElectromagneticConfig(
        include_self_field=electromagnetic_data.get(
            "include_self_field", _ELECTROMAGNETIC_DEFAULTS["include_self_field"]
        ),
        include_displacement_current=electromagnetic_data.get(
            "include_displacement_current",
            _ELECTROMAGNETIC_DEFAULTS["include_displacement_current"],
        ),
        screening_tolerance=float(electromagnetic_data.get(
            "screening_tolerance", _ELECTROMAGNETIC_DEFAULTS["screening_tolerance"]
        )),
        screening_max_iterations=int(electromagnetic_data.get(
            "screening_max_iterations",
            _ELECTROMAGNETIC_DEFAULTS["screening_max_iterations"],
        )),
        screening_step_size=float(electromagnetic_data.get(
            "screening_step_size", _ELECTROMAGNETIC_DEFAULTS["screening_step_size"]
        )),
        screening_source_stride=int(electromagnetic_data.get(
            "screening_source_stride",
            _ELECTROMAGNETIC_DEFAULTS["screening_source_stride"],
        )),
    )

    #
    # TDGL configuration.
    #

    tdgl_data = data.get(
        "tdgl",
        {},
    )

    tdgl = TDGLConfig(
        time_integrator=tdgl_data.get("time_integrator", _TDGL_DEFAULTS["time_integrator"]),
        u=float(
            tdgl_data.get(
                "u",
                _TDGL_DEFAULTS["u"],
            )
        ),

        gamma=float(
            tdgl_data.get(
                "gamma",
                _TDGL_DEFAULTS["gamma"],
            )
        ),

        kappa=float(
            tdgl_data.get(
                "kappa",
                _TDGL_DEFAULTS["kappa"],
            )
        ),

        normalization=tdgl_data.get(
            "normalization",
            _TDGL_DEFAULTS["normalization"],
        ),

        temperature_model=tdgl_data.get(
            "temperature_model",
            _TDGL_DEFAULTS["temperature_model"],
        ),

        max_normalized_timestep=float(
            tdgl_data.get(
                "max_normalized_timestep",
                _TDGL_DEFAULTS["max_normalized_timestep"],
            )
        ),

        include_scalar_potential=tdgl_data.get(
            "include_scalar_potential",
            _TDGL_DEFAULTS["include_scalar_potential"],
        ),

        stability_safety_factor=float(tdgl_data.get(
            "stability_safety_factor", _TDGL_DEFAULTS["stability_safety_factor"])),
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

        laser=laser,

        pinning=pinning,

        josephson=josephson,

        electrical=electrical,

        electromagnetic=electromagnetic,

        tdgl=tdgl,

        coupling=CouplingConfig(**data.get('coupling', {})),

        numerics=NumericsConfig(**data.get('numerics', {})),
    )

    config.validate()

    return config
