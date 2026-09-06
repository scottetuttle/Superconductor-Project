import argparse
import csv
import json
import math
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from shs.config.builder import build_simulation
from shs.physics.thermal import ThermalModel
from shs.physics.electrical import ElectricalModel
from shs.tdgl import TDGLModel, TDGLParameters
from shs.solvers.tdgl_solver import tdgl_step
from shs.solvers.electrical_solver import electrical_step
from shs.solvers.thermal_solver import thermal_step


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

DEFAULT_CONFIG = "configs/simulations/nbn_hotspot_test.json"

DEFAULT_DT_VALUES = (
    1.0e-12,
    5.0e-13,
    1.0e-13,
    5.0e-14,
    1.0e-14,
)

DEFAULT_STEPS = 1

DEFAULT_COUPLING_TOLERANCE = 1.0e-6
DEFAULT_MAX_COUPLING_ITERATIONS = 200

DEFAULT_SIGMAS = (
    1.00,
    0.90,
    0.85,
    0.80,
    0.75,
    0.70,
    0.60,
    0.50,
)

DEFAULT_MIN_SIGMA = 0.25
DEFAULT_SIGMA_REDUCTION = 0.80

DEFAULT_OSCILLATION_WINDOW = 6
DEFAULT_OSCILLATION_REQUIRED = 3

DEFAULT_BATH_TEMPERATURE = 3.0
DEFAULT_THERMAL_RELAXATION_RATE = 0.0

DEFAULT_VOLTAGE_LEFT = 1.0e-3
DEFAULT_VOLTAGE_RIGHT = 0.0

DEFAULT_EXTERNAL_HEAT = 0.0
DEFAULT_HEAT_RADIUS = 0.15

DEFAULT_OUTPUT_DIRECTORY = "benchmark_results_v2"
DEFAULT_CSV_NAME = "coupled_convergence_results_v2.csv"
DEFAULT_JSON_NAME = "coupled_convergence_results_v2.json"


DEFAULT_ORDERINGS = (
    "tdgl_electrical_thermal",
    "tdgl_thermal_electrical",
    "electrical_tdgl_thermal",
    "electrical_thermal_tdgl",
    "thermal_tdgl_electrical",
    "thermal_electrical_tdgl",
)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class FieldSnapshot:
    temperature: np.ndarray
    psi: np.ndarray
    voltage: np.ndarray
    vector_potential_x: np.ndarray
    vector_potential_y: np.ndarray


@dataclass
class CouplingIterationResult:
    iteration: int
    sigma: float
    global_residual: float
    residuals: Dict[str, float]
    oscillation_detected: bool
    sigma_changed: bool


@dataclass
class BenchmarkCaseResult:
    ordering: str
    initial_sigma: float
    final_sigma: float
    dt: float
    steps: int

    converged: bool
    failed_steps: int

    total_runtime: float
    runtime_per_step: float

    average_iterations: float
    maximum_iterations: int

    final_residual: float

    final_temperature: float
    final_temperature_max: float

    final_psi_amplitude_mean: float
    final_psi_amplitude_min: float
    final_psi_amplitude_max: float

    final_current_mean: float
    final_current_max: float

    final_joule_heating_mean: float
    final_joule_heating_max: float

    oscillation_events: int
    sigma_reductions: int

    residual_history: List[float] = field(default_factory=list)
    field_residual_history: List[Dict[str, float]] = field(default_factory=list)
    sigma_history: List[float] = field(default_factory=list)
    oscillation_history: List[bool] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Snapshot / state handling
# ---------------------------------------------------------------------------

def snapshot_state(simulation) -> FieldSnapshot:
    fields = simulation.fields

    return FieldSnapshot(
        temperature=fields.temperature.copy(),
        psi=fields.psi.copy(),
        voltage=fields.voltage.copy(),
        vector_potential_x=fields.vector_potential_x.copy(),
        vector_potential_y=fields.vector_potential_y.copy(),
    )


def restore_state(simulation, snapshot: FieldSnapshot) -> None:
    fields = simulation.fields

    fields.temperature = snapshot.temperature.copy()
    fields.psi = snapshot.psi.copy()
    fields.voltage = snapshot.voltage.copy()
    fields.vector_potential_x = snapshot.vector_potential_x.copy()
    fields.vector_potential_y = snapshot.vector_potential_y.copy()


def copy_coupling_state(simulation) -> FieldSnapshot:
    return snapshot_state(simulation)


def relax_state(
    old_state: FieldSnapshot,
    new_state: FieldSnapshot,
    sigma: float,
) -> FieldSnapshot:
    return FieldSnapshot(
        temperature=(
            old_state.temperature
            + sigma * (new_state.temperature - old_state.temperature)
        ),
        psi=(
            old_state.psi
            + sigma * (new_state.psi - old_state.psi)
        ),
        voltage=(
            old_state.voltage
            + sigma * (new_state.voltage - old_state.voltage)
        ),
        vector_potential_x=(
            old_state.vector_potential_x
            + sigma
            * (
                new_state.vector_potential_x
                - old_state.vector_potential_x
            )
        ),
        vector_potential_y=(
            old_state.vector_potential_y
            + sigma
            * (
                new_state.vector_potential_y
                - old_state.vector_potential_y
            )
        ),
    )


def apply_snapshot(simulation, snapshot: FieldSnapshot) -> None:
    restore_state(simulation, snapshot)


# ---------------------------------------------------------------------------
# Residual calculations
# ---------------------------------------------------------------------------

def relative_residual(
    old: np.ndarray,
    new: np.ndarray,
    floor: float = 1.0e-12,
) -> float:
    numerator = np.linalg.norm(new - old)
    denominator = max(np.linalg.norm(old), floor)

    return float(numerator / denominator)


def field_residuals(
    old: FieldSnapshot,
    new: FieldSnapshot,
) -> Dict[str, float]:
    return {
        "temperature": relative_residual(
            old.temperature,
            new.temperature,
        ),
        "psi": relative_residual(
            old.psi,
            new.psi,
        ),
        "voltage": relative_residual(
            old.voltage,
            new.voltage,
        ),
        "vector_potential_x": relative_residual(
            old.vector_potential_x,
            new.vector_potential_x,
        ),
        "vector_potential_y": relative_residual(
            old.vector_potential_y,
            new.vector_potential_y,
        ),
    }


def global_residual(residuals: Dict[str, float]) -> float:
    finite_values = [
        value
        for value in residuals.values()
        if np.isfinite(value)
    ]

    if not finite_values:
        return float("inf")

    return max(finite_values)


# ---------------------------------------------------------------------------
# Oscillation detection
# ---------------------------------------------------------------------------

def update_direction(
    previous: FieldSnapshot,
    current: FieldSnapshot,
) -> Dict[str, np.ndarray]:
    return {
        "temperature": current.temperature - previous.temperature,
        "psi_real": current.psi.real - previous.psi.real,
        "psi_imag": current.psi.imag - previous.psi.imag,
        "voltage": current.voltage - previous.voltage,
        "vector_potential_x": (
            current.vector_potential_x
            - previous.vector_potential_x
        ),
        "vector_potential_y": (
            current.vector_potential_y
            - previous.vector_potential_y
        ),
    }


def cosine_between(
    a: np.ndarray,
    b: np.ndarray,
    floor: float = 1.0e-30,
) -> float:
    a_norm = np.linalg.norm(a)
    b_norm = np.linalg.norm(b)

    if a_norm <= floor or b_norm <= floor:
        return 1.0

    return float(
        np.vdot(a, b).real
        / (a_norm * b_norm)
    )


def detect_field_oscillation(
    previous_direction: Optional[Dict[str, np.ndarray]],
    current_direction: Dict[str, np.ndarray],
    cosine_threshold: float = -0.2,
) -> Tuple[bool, Dict[str, float]]:
    """
    Detect reversal of the update direction.

    A strongly negative cosine means that the new update is pointing
    opposite to the previous update, which is a useful indicator of
    oscillatory fixed-point behavior.
    """

    if previous_direction is None:
        return False, {}

    directional_cosines = {}

    for name in current_direction:
        if name not in previous_direction:
            continue

        directional_cosines[name] = cosine_between(
            previous_direction[name],
            current_direction[name],
        )

    oscillating_fields = [
        name
        for name, cosine in directional_cosines.items()
        if cosine < cosine_threshold
    ]

    return bool(oscillating_fields), directional_cosines


# ---------------------------------------------------------------------------
# Physics model creation
# ---------------------------------------------------------------------------

def create_models(
    bath_temperature: float,
    thermal_relaxation_rate: float,
) -> Tuple[TDGLModel, ThermalModel, ElectricalModel]:

    tdgl_model = TDGLModel(
        TDGLParameters()
    )

    thermal_model = ThermalModel(
        bath_temperature=bath_temperature,
        thermal_relaxation_rate=thermal_relaxation_rate,
    )

    electrical_model = ElectricalModel()

    return tdgl_model, thermal_model, electrical_model


# ---------------------------------------------------------------------------
# External thermal perturbation
# ---------------------------------------------------------------------------

def create_external_heat_source(
    simulation,
    amplitude: float,
    radius_fraction: float,
) -> np.ndarray:
    """
    Create a smooth localized heat source centered on the mesh.

    This is deliberately implemented at the benchmark level. It does
    not modify the production thermal solver or introduce a direct
    perturbation to psi.
    """

    if amplitude == 0.0:
        return np.zeros_like(
            simulation.fields.temperature
        )

    mesh = simulation.mesh

    x = np.arange(mesh.nx) * mesh.dx
    y = np.arange(mesh.ny) * mesh.dy

    xx, yy = np.meshgrid(x, y)

    x_center = 0.5 * (x[0] + x[-1])
    y_center = 0.5 * (y[0] + y[-1])

    domain_x = max(x[-1] - x[0], mesh.dx)
    domain_y = max(y[-1] - y[0], mesh.dy)

    radius = radius_fraction * min(domain_x, domain_y)

    gaussian = np.exp(
        -(
            (xx - x_center) ** 2
            + (yy - y_center) ** 2
        )
        / (2.0 * radius**2)
    )

    return amplitude * gaussian


# ---------------------------------------------------------------------------
# Individual physics operations
# ---------------------------------------------------------------------------

def run_tdgl(
    simulation,
    accepted_state: FieldSnapshot,
    dt: float,
    tdgl_model: TDGLModel,
) -> None:
    """
    Advance TDGL by one physical timestep starting from the accepted
    state at t_n, while retaining the current coupling guess for the
    other fields.
    """

    fields = simulation.fields

    trial_temperature = fields.temperature.copy()
    trial_Ax = fields.vector_potential_x.copy()
    trial_Ay = fields.vector_potential_y.copy()

    fields.temperature = trial_temperature
    fields.vector_potential_x = trial_Ax
    fields.vector_potential_y = trial_Ay

    # TDGL advances psi from the accepted physical state.
    fields.psi = accepted_state.psi.copy()

    tdgl_step(
        simulation,
        dt,
        tdgl_model,
    )


def run_electrical(
    simulation,
    electrical_model: ElectricalModel,
) -> None:
    fields = simulation.fields

    superconducting_fraction = np.abs(fields.psi) ** 2

    electrical_step(
        fields,
        simulation.mesh,
        simulation.material_map,
        simulation.contact_map,
        voltage_left=simulation._benchmark_voltage_left,
        voltage_right=simulation._benchmark_voltage_right,
        superconducting_fraction=superconducting_fraction,
        superconducting_current_x=fields.supercurrent_density_x,
        superconducting_current_y=fields.supercurrent_density_y,
    )


def run_thermal(
    simulation,
    dt: float,
    thermal_model: ThermalModel,
    external_heat_source: Optional[np.ndarray],
) -> None:
    fields = simulation.fields

    if external_heat_source is not None:
        fields.heat_source = (
            fields.heat_source
            + external_heat_source
        )

    thermal_step(
        simulation,
        dt,
        thermal_model,
    )


def run_operation(
    operation: str,
    simulation,
    accepted_state: FieldSnapshot,
    dt: float,
    tdgl_model: TDGLModel,
    thermal_model: ThermalModel,
    electrical_model: ElectricalModel,
    external_heat_source: Optional[np.ndarray],
) -> None:

    if operation == "tdgl":
        run_tdgl(
            simulation,
            accepted_state,
            dt,
            tdgl_model,
        )

    elif operation == "electrical":
        run_electrical(
            simulation,
            electrical_model,
        )

    elif operation == "thermal":
        run_thermal(
            simulation,
            dt,
            thermal_model,
            external_heat_source,
        )

    else:
        raise ValueError(
            f"Unknown coupling operation: {operation}"
        )


# ---------------------------------------------------------------------------
# Ordering utilities
# ---------------------------------------------------------------------------

def parse_ordering(ordering: str) -> List[str]:
    pieces = ordering.split("_")

    expected = {"tdgl", "electrical", "thermal"}

    if len(pieces) != 3:
        raise ValueError(
            f"Invalid ordering: {ordering}"
        )

    if set(pieces) != expected:
        raise ValueError(
            f"Ordering must contain tdgl, electrical, thermal: {ordering}"
        )

    return pieces


# ---------------------------------------------------------------------------
# Candidate coupled state
# ---------------------------------------------------------------------------

def calculate_candidate(
    simulation,
    accepted_state: FieldSnapshot,
    previous_guess: FieldSnapshot,
    ordering: Sequence[str],
    dt: float,
    tdgl_model: TDGLModel,
    thermal_model: ThermalModel,
    electrical_model: ElectricalModel,
    external_heat_source: Optional[np.ndarray],
) -> FieldSnapshot:

    # Begin from the current Picard guess.
    apply_snapshot(
        simulation,
        previous_guess,
    )

    for operation in ordering:
        run_operation(
            operation,
            simulation,
            accepted_state,
            dt,
            tdgl_model,
            thermal_model,
            electrical_model,
            external_heat_source,
        )

    return copy_coupling_state(simulation)


# ---------------------------------------------------------------------------
# One physical timestep
# ---------------------------------------------------------------------------

def coupled_physical_timestep(
    simulation,
    accepted_state: FieldSnapshot,
    dt: float,
    ordering: Sequence[str],
    initial_sigma: float,
    tolerance: float,
    max_iterations: int,
    min_sigma: float,
    sigma_reduction: float,
    oscillation_window: int,
    oscillation_required: int,
    tdgl_model: TDGLModel,
    thermal_model: ThermalModel,
    electrical_model: ElectricalModel,
    external_heat_source: Optional[np.ndarray],
) -> Tuple[
    bool,
    FieldSnapshot,
    float,
    int,
    List[CouplingIterationResult],
]:

    sigma = initial_sigma

    previous_guess = accepted_state

    previous_direction = None
    oscillation_count = 0

    iteration_results = []

    for iteration in range(1, max_iterations + 1):

        candidate = calculate_candidate(
            simulation,
            accepted_state,
            previous_guess,
            ordering,
            dt,
            tdgl_model,
            thermal_model,
            electrical_model,
            external_heat_source,
        )

        residuals = field_residuals(
            previous_guess,
            candidate,
        )

        residual = global_residual(residuals)

        current_direction = update_direction(
            previous_guess,
            candidate,
        )

        oscillating, directional_cosines = detect_field_oscillation(
            previous_direction,
            current_direction,
        )

        if oscillating:
            oscillation_count += 1
        else:
            oscillation_count = 0

        sigma_changed = False

        # Require several consecutive oscillatory iterations before
        # changing sigma. This prevents one noisy iteration from
        # triggering unnecessary relaxation.
        if (
            oscillation_count >= oscillation_required
            and sigma > min_sigma
        ):
            new_sigma = max(
                min_sigma,
                sigma * sigma_reduction,
            )

            if new_sigma < sigma:
                sigma = new_sigma
                sigma_changed = True
                oscillation_count = 0

        iteration_results.append(
            CouplingIterationResult(
                iteration=iteration,
                sigma=sigma,
                global_residual=residual,
                residuals=residuals,
                oscillation_detected=oscillating,
                sigma_changed=sigma_changed,
            )
        )

        # Convergence is based on the actual Picard-map residual
        # before relaxation.
        if residual <= tolerance:
            apply_snapshot(
                simulation,
                candidate,
            )

            return (
                True,
                candidate,
                sigma,
                iteration,
                iteration_results,
            )

        # Relax toward the candidate.
        relaxed = relax_state(
            previous_guess,
            candidate,
            sigma,
        )

        previous_direction = current_direction
        previous_guess = relaxed

    # Failed timestep: restore accepted physical state.
    apply_snapshot(
        simulation,
        accepted_state,
    )

    final_state = accepted_state

    return (
        False,
        final_state,
        sigma,
        max_iterations,
        iteration_results,
    )


# ---------------------------------------------------------------------------
# Final observables
# ---------------------------------------------------------------------------

def final_observables(simulation) -> Dict[str, float]:
    fields = simulation.fields

    psi_amplitude = np.abs(fields.psi)

    current_magnitude = np.sqrt(
        fields.current_density_x**2
        + fields.current_density_y**2
    )

    heating = fields.heat_source

    return {
        "final_temperature": float(
            np.mean(fields.temperature)
        ),
        "final_temperature_max": float(
            np.max(fields.temperature)
        ),
        "final_psi_amplitude_mean": float(
            np.mean(psi_amplitude)
        ),
        "final_psi_amplitude_min": float(
            np.min(psi_amplitude)
        ),
        "final_psi_amplitude_max": float(
            np.max(psi_amplitude)
        ),
        "final_current_mean": float(
            np.mean(current_magnitude)
        ),
        "final_current_max": float(
            np.max(current_magnitude)
        ),
        "final_joule_heating_mean": float(
            np.mean(heating)
        ),
        "final_joule_heating_max": float(
            np.max(heating)
        ),
    }


# ---------------------------------------------------------------------------
# One benchmark case
# ---------------------------------------------------------------------------

def run_benchmark_case(
    config_path: str,
    dt: float,
    steps: int,
    ordering: str,
    sigma: float,
    tolerance: float,
    max_iterations: int,
    min_sigma: float,
    sigma_reduction: float,
    oscillation_window: int,
    oscillation_required: int,
    bath_temperature: float,
    thermal_relaxation_rate: float,
    voltage_left: float,
    voltage_right: float,
    external_heat: float,
    heat_radius: float,
) -> BenchmarkCaseResult:

    simulation = build_simulation(
        config_path
    )

    simulation._benchmark_voltage_left = voltage_left
    simulation._benchmark_voltage_right = voltage_right

    tdgl_model, thermal_model, electrical_model = create_models(
        bath_temperature,
        thermal_relaxation_rate,
    )

    external_heat_source = create_external_heat_source(
        simulation,
        external_heat,
        heat_radius,
    )

    ordering_sequence = parse_ordering(ordering)

    accepted_state = snapshot_state(
        simulation
    )

    total_iterations = 0
    maximum_iterations = 0

    failed_steps = 0

    all_residual_history = []
    all_field_residual_history = []
    all_sigma_history = []
    all_oscillation_history = []

    oscillation_events = 0
    sigma_reductions = 0

    current_sigma = sigma

    start_time = time.perf_counter()

    final_residual = float("inf")

    for step in range(steps):

        (
            converged,
            final_state,
            final_sigma,
            iterations,
            iteration_results,
        ) = coupled_physical_timestep(
            simulation=simulation,
            accepted_state=accepted_state,
            dt=dt,
            ordering=ordering_sequence,
            initial_sigma=current_sigma,
            tolerance=tolerance,
            max_iterations=max_iterations,
            min_sigma=min_sigma,
            sigma_reduction=sigma_reduction,
            oscillation_window=oscillation_window,
            oscillation_required=oscillation_required,
            tdgl_model=tdgl_model,
            thermal_model=thermal_model,
            electrical_model=electrical_model,
            external_heat_source=external_heat_source,
        )

        total_iterations += iterations
        maximum_iterations = max(
            maximum_iterations,
            iterations,
        )

        for result in iteration_results:

            all_residual_history.append(
                result.global_residual
            )

            all_field_residual_history.append(
                result.residuals
            )

            all_sigma_history.append(
                result.sigma
            )

            all_oscillation_history.append(
                result.oscillation_detected
            )

            if result.oscillation_detected:
                oscillation_events += 1

            if result.sigma_changed:
                sigma_reductions += 1

        if iteration_results:
            final_residual = (
                iteration_results[-1].global_residual
            )

        if not converged:
            failed_steps += 1
            break

        accepted_state = snapshot_state(
            simulation
        )

        current_sigma = final_sigma

    total_runtime = (
        time.perf_counter()
        - start_time
    )

    successful_steps = steps - failed_steps

    if successful_steps > 0:
        runtime_per_step = (
            total_runtime
            / successful_steps
        )
    else:
        runtime_per_step = float("nan")

    average_iterations = (
        total_iterations / steps
        if steps > 0
        else float("nan")
    )

    observables = final_observables(
        simulation
    )

    return BenchmarkCaseResult(
        ordering=ordering,
        initial_sigma=sigma,
        final_sigma=current_sigma,
        dt=dt,
        steps=steps,

        converged=(
            failed_steps == 0
            and final_residual <= tolerance
        ),
        failed_steps=failed_steps,

        total_runtime=total_runtime,
        runtime_per_step=runtime_per_step,

        average_iterations=average_iterations,
        maximum_iterations=maximum_iterations,

        final_residual=final_residual,

        final_temperature=observables[
            "final_temperature"
        ],
        final_temperature_max=observables[
            "final_temperature_max"
        ],

        final_psi_amplitude_mean=observables[
            "final_psi_amplitude_mean"
        ],
        final_psi_amplitude_min=observables[
            "final_psi_amplitude_min"
        ],
        final_psi_amplitude_max=observables[
            "final_psi_amplitude_max"
        ],

        final_current_mean=observables[
            "final_current_mean"
        ],
        final_current_max=observables[
            "final_current_max"
        ],

        final_joule_heating_mean=observables[
            "final_joule_heating_mean"
        ],
        final_joule_heating_max=observables[
            "final_joule_heating_max"
        ],

        oscillation_events=oscillation_events,
        sigma_reductions=sigma_reductions,

        residual_history=all_residual_history,
        field_residual_history=all_field_residual_history,
        sigma_history=all_sigma_history,
        oscillation_history=all_oscillation_history,
    )


# ---------------------------------------------------------------------------
# Sweep
# ---------------------------------------------------------------------------

def run_sweep(
    config_path: str,
    dt_values: Sequence[float],
    steps: int,
    orderings: Sequence[str],
    sigmas: Sequence[float],
    tolerance: float,
    max_iterations: int,
    min_sigma: float,
    sigma_reduction: float,
    oscillation_window: int,
    oscillation_required: int,
    bath_temperature: float,
    thermal_relaxation_rate: float,
    voltage_left: float,
    voltage_right: float,
    external_heat: float,
    heat_radius: float,
) -> List[BenchmarkCaseResult]:

    results = []

    total_cases = (
        len(dt_values)
        * len(orderings)
        * len(sigmas)
    )

    case_number = 0

    for dt in dt_values:
        for ordering in orderings:
            for sigma in sigmas:

                case_number += 1

                print(
                    f"\n[{case_number}/{total_cases}] "
                    f"dt={dt:.3e} "
                    f"{ordering} "
                    f"sigma={sigma:.2f}"
                )

                result = run_benchmark_case(
                    config_path=config_path,
                    dt=dt,
                    steps=steps,
                    ordering=ordering,
                    sigma=sigma,
                    tolerance=tolerance,
                    max_iterations=max_iterations,
                    min_sigma=min_sigma,
                    sigma_reduction=sigma_reduction,
                    oscillation_window=oscillation_window,
                    oscillation_required=oscillation_required,
                    bath_temperature=bath_temperature,
                    thermal_relaxation_rate=thermal_relaxation_rate,
                    voltage_left=voltage_left,
                    voltage_right=voltage_right,
                    external_heat=external_heat,
                    heat_radius=heat_radius,
                )

                results.append(result)

                print(
                    f"    status={
                    "CONVERGED"
                    if result.converged
                    else "FAILED"
                    } "
                    f"iterations={result.maximum_iterations} "
                    f"residual={result.final_residual:.3e} "
                    f"final_sigma={result.final_sigma:.3f} "
                    f"oscillations={result.oscillation_events} "
                    f"runtime={result.total_runtime:.2f}s"
                )

    return results


# ---------------------------------------------------------------------------
# Serialization
# ---------------------------------------------------------------------------

def result_to_dict(
    result: BenchmarkCaseResult,
) -> Dict:

    return {
        "ordering": result.ordering,
        "initial_sigma": result.initial_sigma,
        "final_sigma": result.final_sigma,

        "dt": result.dt,
        "steps": result.steps,

        "converged": result.converged,
        "failed_steps": result.failed_steps,

        "total_runtime": result.total_runtime,
        "runtime_per_step": result.runtime_per_step,

        "average_iterations": result.average_iterations,
        "maximum_iterations": result.maximum_iterations,

        "final_residual": result.final_residual,

        "final_temperature": result.final_temperature,
        "final_temperature_max": result.final_temperature_max,

        "final_psi_amplitude_mean":
            result.final_psi_amplitude_mean,
        "final_psi_amplitude_min":
            result.final_psi_amplitude_min,
        "final_psi_amplitude_max":
            result.final_psi_amplitude_max,

        "final_current_mean":
            result.final_current_mean,
        "final_current_max":
            result.final_current_max,

        "final_joule_heating_mean":
            result.final_joule_heating_mean,
        "final_joule_heating_max":
            result.final_joule_heating_max,

        "oscillation_events":
            result.oscillation_events,
        "sigma_reductions":
            result.sigma_reductions,

        "residual_history":
            result.residual_history,
        "field_residual_history":
            result.field_residual_history,
        "sigma_history":
            result.sigma_history,
        "oscillation_history":
            result.oscillation_history,
    }


def write_json(
    results: Sequence[BenchmarkCaseResult],
    output_path: Path,
    metadata: Dict,
) -> None:

    payload = {
        "metadata": metadata,
        "results": [
            result_to_dict(result)
            for result in results
        ],
    }

    with output_path.open(
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(
            payload,
            handle,
            indent=2,
            allow_nan=True,
        )


def write_csv(
    results: Sequence[BenchmarkCaseResult],
    output_path: Path,
) -> None:

    fieldnames = [
        "ordering",
        "initial_sigma",
        "final_sigma",
        "dt",
        "steps",
        "converged",
        "failed_steps",
        "total_runtime",
        "runtime_per_step",
        "average_iterations",
        "maximum_iterations",
        "final_residual",
        "final_temperature",
        "final_temperature_max",
        "final_psi_amplitude_mean",
        "final_psi_amplitude_min",
        "final_psi_amplitude_max",
        "final_current_mean",
        "final_current_max",
        "final_joule_heating_mean",
        "final_joule_heating_max",
        "oscillation_events",
        "sigma_reductions",
    ]

    with output_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )

        writer.writeheader()

        for result in results:
            row = result_to_dict(result)

            writer.writerow({
                key: row[key]
                for key in fieldnames
            })


# ---------------------------------------------------------------------------
# Terminal summary
# ---------------------------------------------------------------------------

def print_summary(
    results: Sequence[BenchmarkCaseResult],
) -> None:

    print()
    print("=" * 125)
    print(
        f"{'dt':<12}"
        f"{'Ordering':<35}"
        f"{'Sigma':>8}"
        f"{'Final σ':>9}"
        f"{'Status':>11}"
        f"{'Avg Iter':>11}"
        f"{'Runtime':>12}"
        f"{'Residual':>14}"
        f"{'Oscill.':>9}"
    )
    print("-" * 125)

    for result in results:

        status = (
            "CONVERGED"
            if result.converged
            else "FAILED"
        )

        print(
            f"{result.dt:<12.3e}"
            f"{result.ordering:<35}"
            f"{result.initial_sigma:>8.2f}"
            f"{result.final_sigma:>9.3f}"
            f"{status:>11}"
            f"{result.average_iterations:>11.2f}"
            f"{result.total_runtime:>12.3f}"
            f"{result.final_residual:>14.3e}"
            f"{result.oscillation_events:>9}"
        )

    print("=" * 125)


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_float_list(values: Sequence[str]) -> List[float]:
    result = []

    for value in values:
        for piece in value.split(","):
            result.append(float(piece))

    return result


def build_parser() -> argparse.ArgumentParser:

    parser = argparse.ArgumentParser(
        description=(
            "Coupled TDGL/electrical/thermal convergence "
            "benchmark with timestep and adaptive-relaxation sweeps."
        )
    )

    parser.add_argument(
        "--config",
        default=DEFAULT_CONFIG,
    )

    parser.add_argument(
        "--dt",
        nargs="+",
        default=[
            str(value)
            for value in DEFAULT_DT_VALUES
        ],
        help="One or more physical timesteps.",
    )

    parser.add_argument(
        "--steps",
        type=int,
        default=DEFAULT_STEPS,
    )

    parser.add_argument(
        "--sigma",
        nargs="+",
        default=[
            str(value)
            for value in DEFAULT_SIGMAS
        ],
        help="Initial coupling relaxation values.",
    )

    parser.add_argument(
        "--ordering",
        nargs="+",
        default=list(DEFAULT_ORDERINGS),
    )

    parser.add_argument(
        "--tolerance",
        type=float,
        default=DEFAULT_COUPLING_TOLERANCE,
    )

    parser.add_argument(
        "--max-iterations",
        type=int,
        default=DEFAULT_MAX_COUPLING_ITERATIONS,
    )

    parser.add_argument(
        "--min-sigma",
        type=float,
        default=DEFAULT_MIN_SIGMA,
    )

    parser.add_argument(
        "--sigma-reduction",
        type=float,
        default=DEFAULT_SIGMA_REDUCTION,
        help=(
            "Multiplier applied when oscillation causes "
            "automatic sigma reduction."
        ),
    )

    parser.add_argument(
        "--oscillation-window",
        type=int,
        default=DEFAULT_OSCILLATION_WINDOW,
        help=(
            "Reserved diagnostic parameter controlling "
            "the intended observation window."
        ),
    )

    parser.add_argument(
        "--oscillation-required",
        type=int,
        default=DEFAULT_OSCILLATION_REQUIRED,
        help=(
            "Consecutive oscillatory iterations required "
            "before sigma is reduced."
        ),
    )

    parser.add_argument(
        "--bath-temperature",
        type=float,
        default=DEFAULT_BATH_TEMPERATURE,
    )

    parser.add_argument(
        "--thermal-relaxation-rate",
        type=float,
        default=DEFAULT_THERMAL_RELAXATION_RATE,
    )

    parser.add_argument(
        "--voltage-left",
        type=float,
        default=DEFAULT_VOLTAGE_LEFT,
    )

    parser.add_argument(
        "--voltage-right",
        type=float,
        default=DEFAULT_VOLTAGE_RIGHT,
    )

    parser.add_argument(
        "--external-heat",
        type=float,
        default=DEFAULT_EXTERNAL_HEAT,
        help=(
            "Localized benchmark heat-source amplitude. "
            "Zero preserves the baseline case."
        ),
    )

    parser.add_argument(
        "--heat-radius",
        type=float,
        default=DEFAULT_HEAT_RADIUS,
        help=(
            "Localized heat radius as a fraction of the "
            "smaller domain dimension."
        ),
    )

    parser.add_argument(
        "--output-dir",
        default=DEFAULT_OUTPUT_DIRECTORY,
    )

    parser.add_argument(
        "--csv-name",
        default=DEFAULT_CSV_NAME,
    )

    parser.add_argument(
        "--json-name",
        default=DEFAULT_JSON_NAME,
    )

    parser.add_argument(
        "--no-output",
        action="store_true",
    )

    return parser


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:

    parser = build_parser()
    args = parser.parse_args()

    dt_values = parse_float_list(
        args.dt
    )

    sigmas = parse_float_list(
        args.sigma
    )

    if args.steps <= 0:
        parser.error(
            "--steps must be positive."
        )

    if args.tolerance <= 0:
        parser.error(
            "--tolerance must be positive."
        )

    if args.max_iterations <= 0:
        parser.error(
            "--max-iterations must be positive."
        )

    if not 0.0 < args.min_sigma <= 1.0:
        parser.error(
            "--min-sigma must be in (0, 1]."
        )

    if not 0.0 < args.sigma_reduction < 1.0:
        parser.error(
            "--sigma-reduction must be between 0 and 1."
        )

    for sigma in sigmas:
        if not 0.0 < sigma <= 1.0:
            parser.error(
                "All sigma values must be in (0, 1]."
            )

    for dt in dt_values:
        if dt <= 0:
            parser.error(
                "All dt values must be positive."
            )

    results = run_sweep(
        config_path=args.config,
        dt_values=dt_values,
        steps=args.steps,
        orderings=args.ordering,
        sigmas=sigmas,
        tolerance=args.tolerance,
        max_iterations=args.max_iterations,
        min_sigma=args.min_sigma,
        sigma_reduction=args.sigma_reduction,
        oscillation_window=args.oscillation_window,
        oscillation_required=args.oscillation_required,
        bath_temperature=args.bath_temperature,
        thermal_relaxation_rate=args.thermal_relaxation_rate,
        voltage_left=args.voltage_left,
        voltage_right=args.voltage_right,
        external_heat=args.external_heat,
        heat_radius=args.heat_radius,
    )

    print_summary(
        results
    )

    if not args.no_output:

        output_directory = Path(
            args.output_dir
        )

        output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        metadata = {
            "configuration": args.config,
            "dt_values": dt_values,
            "steps": args.steps,
            "coupling_tolerance": args.tolerance,
            "max_coupling_iterations":
                args.max_iterations,
            "initial_sigmas": sigmas,
            "min_sigma": args.min_sigma,
            "sigma_reduction":
                args.sigma_reduction,
            "oscillation_window":
                args.oscillation_window,
            "oscillation_required":
                args.oscillation_required,
            "orderings": args.ordering,
            "voltage_left":
                args.voltage_left,
            "voltage_right":
                args.voltage_right,
            "bath_temperature":
                args.bath_temperature,
            "thermal_relaxation_rate":
                args.thermal_relaxation_rate,
            "external_heat":
                args.external_heat,
            "heat_radius":
                args.heat_radius,
        }

        write_csv(
            results,
            output_directory / args.csv_name,
        )

        write_json(
            results,
            output_directory / args.json_name,
            metadata,
        )

        print()
        print(
            f"Results written to: "
            f"{output_directory}"
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())