"""
Coupled convergence benchmark for the Superconducting Hotspot Simulator.

This benchmark investigates fixed-point convergence of the coupled:

    TDGL
    Electrical transport
    Thermal evolution

system.

The benchmark is intentionally separate from the production coupled
solvers. It is a developer/research tool for determining:

    1. Which physics ordering converges most reliably.
    2. How the coupling relaxation parameter affects convergence.
    3. How many coupling iterations are required.
    4. How convergence behavior affects runtime.
    5. Whether different successful orderings produce comparable states.

IMPORTANT NUMERICAL DISTINCTION
--------------------------------

A physical timestep is:

    state_n -> state_(n+1)

The coupling iteration must NOT repeatedly advance physical time.

Instead, for one physical timestep:

    accepted state at t_n
            |
            v
       initial guess
            |
            v
       physics sweep
            |
            v
       candidate t_(n+1)
            |
            v
       relaxation
            |
            v
       next coupling guess
            |
           ...
            |
            v
       converged t_(n+1)

The time-evolving variables are therefore always anchored to their accepted
values at t_n when their respective solver is called:

    TDGL:
        psi_n -> psi_(n+1)

    Thermal:
        T_n -> T_(n+1)

Meanwhile, the other coupled fields are supplied from the current Picard
iteration guess.

RELAXATION
----------

The coupling relaxation parameter is:

    x_new = (1 - sigma) * x_old + sigma * x_candidate

where:

    sigma = 1.00
        Full update.

    sigma = 0.50
        Halfway between the previous guess and candidate.

    sigma = 0.10
        Strong damping.

This is the coupling relaxation parameter. It is NOT electrical
conductivity.

DEFAULT BENCHMARK
-----------------

The default sweep contains:

    6 physics orderings
    x
    6 relaxation values

for:

    36 benchmark cases.

Default relaxation values:

    1.00
    0.75
    0.50
    0.25
    0.10
    0.05

Default orderings:

    tdgl_electrical_thermal
    tdgl_thermal_electrical
    electrical_tdgl_thermal
    electrical_thermal_tdgl
    thermal_tdgl_electrical
    thermal_electrical_tdgl

EXAMPLES
--------

Run the complete sweep:

    python tools/coupled_convergence_benchmark.py

Run only sigma = 0.50:

    python tools/coupled_convergence_benchmark.py --sigma 0.50

Run two sigma values:

    python tools/coupled_convergence_benchmark.py --sigma 0.25 0.50

Run one ordering:

    python tools/coupled_convergence_benchmark.py \
        --ordering tdgl_electrical_thermal

Run a small test:

    python tools/coupled_convergence_benchmark.py \
        --steps 2 \
        --max-iterations 25

Change timestep:

    python tools/coupled_convergence_benchmark.py \
        --dt 1e-13

Change convergence tolerance:

    python tools/coupled_convergence_benchmark.py \
        --tolerance 1e-6

Change simulation configuration:

    python tools/coupled_convergence_benchmark.py \
        --config configs/simulations/nbn_hotspot_test.json

Results are written to:

    benchmark_results/coupled_convergence_results.csv
    benchmark_results/coupled_convergence_results.json
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from shs.config.builder import build_simulation

from shs.physics.thermal import ThermalModel
from shs.physics.electrical import ElectricalModel

from shs.tdgl import (
    TDGLModel,
    TDGLParameters,
)

from shs.solvers.tdgl_solver import tdgl_step
from shs.solvers.electrical_solver import electrical_step
from shs.solvers.thermal_solver import thermal_step


# ============================================================================
# DEFAULT CONFIGURATION
# ============================================================================

DEFAULT_CONFIG = (
    "configs/simulations/nbn_hotspot_test.json"
)

DEFAULT_DT = 1.0e-13

DEFAULT_STEPS = 5

DEFAULT_COUPLING_TOLERANCE = 1.0e-6

DEFAULT_MAX_COUPLING_ITERATIONS = 100

DEFAULT_SIGMAS = (
    1.00,
    0.75,
    0.50,
    0.25,
    0.10,
    0.05,
)

DEFAULT_ORDERINGS = (
    "tdgl_electrical_thermal",
    "tdgl_thermal_electrical",
    "electrical_tdgl_thermal",
    "electrical_thermal_tdgl",
    "thermal_tdgl_electrical",
    "thermal_electrical_tdgl",
)

DEFAULT_BATH_TEMPERATURE = 3.0

DEFAULT_THERMAL_RELAXATION_RATE = 0.0

DEFAULT_VOLTAGE_LEFT = 1.0e-3

DEFAULT_VOLTAGE_RIGHT = 0.0

DEFAULT_OUTPUT_DIRECTORY = (
    "benchmark_results"
)

DEFAULT_CSV_NAME = (
    "coupled_convergence_results.csv"
)

DEFAULT_JSON_NAME = (
    "coupled_convergence_results.json"
)


# ============================================================================
# DATA STRUCTURES
# ============================================================================


@dataclass
class FieldSnapshot:
    """
    Snapshot of all ndarray fields in a Simulation.
    """

    arrays: Dict[str, np.ndarray]


@dataclass
class CouplingIterationResult:
    """
    Result of one physical timestep.
    """

    converged: bool
    iterations: int
    final_residual: float
    residual_history: List[float]
    failed: bool


@dataclass
class BenchmarkCaseResult:
    """
    Complete result for one ordering / sigma combination.
    """

    ordering: str
    sigma: float

    steps: int
    dt: float

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

    residual_history: List[float]


# ============================================================================
# FIELD SNAPSHOT UTILITIES
# ============================================================================


def snapshot_fields(fields) -> FieldSnapshot:
    """
    Copy every ndarray contained in the Fields object.
    """

    arrays = {}

    for name, value in vars(fields).items():

        if isinstance(value, np.ndarray):
            arrays[name] = value.copy()

    return FieldSnapshot(
        arrays=arrays
    )


def restore_fields(
    fields,
    snapshot: FieldSnapshot,
) -> None:
    """
    Restore every field from a snapshot.
    """

    for name, array in snapshot.arrays.items():

        setattr(
            fields,
            name,
            array.copy(),
        )


def copy_snapshot(
    snapshot: FieldSnapshot,
) -> FieldSnapshot:
    """
    Make a completely independent field snapshot.
    """

    return FieldSnapshot(
        arrays={
            name: array.copy()
            for name, array in snapshot.arrays.items()
        }
    )


# ============================================================================
# RELAXATION
# ============================================================================


def relax_array(
    previous: np.ndarray,
    candidate: np.ndarray,
    sigma: float,
) -> np.ndarray:
    """
    Relax a candidate field toward the previous coupling guess.

        x_new = (1 - sigma) x_old + sigma x_candidate
    """

    return (
        (1.0 - sigma) * previous
        + sigma * candidate
    )


def relaxed_snapshot(
    previous: FieldSnapshot,
    candidate: FieldSnapshot,
    sigma: float,
) -> FieldSnapshot:
    """
    Relax all fields in a candidate state.
    """

    arrays = {}

    field_names = (
        set(previous.arrays)
        | set(candidate.arrays)
    )

    for name in field_names:

        if name not in candidate.arrays:

            arrays[name] = (
                previous.arrays[name].copy()
            )

            continue

        if name not in previous.arrays:

            arrays[name] = (
                candidate.arrays[name].copy()
            )

            continue

        arrays[name] = relax_array(
            previous.arrays[name],
            candidate.arrays[name],
            sigma,
        )

    return FieldSnapshot(
        arrays=arrays
    )


# ============================================================================
# CONVERGENCE METRICS
# ============================================================================


def relative_norm(
    new: np.ndarray,
    old: np.ndarray,
    floor: float = 1.0e-12,
) -> float:
    """
    Calculate a relative L2 update norm.

        ||new - old|| / max(||old||, floor)
    """

    difference = new - old

    numerator = np.linalg.norm(
        difference.ravel()
    )

    denominator = max(
        np.linalg.norm(old.ravel()),
        floor,
    )

    value = numerator / denominator

    if not np.isfinite(value):
        return math.inf

    return float(value)


def state_residual(
    previous: FieldSnapshot,
    candidate: FieldSnapshot,
) -> float:
    """
    Calculate the overall Picard update residual.

    The maximum relative change among the important coupled fields is used.

    This gives each field a dimensionless relative measure instead of
    allowing a dimensional quantity such as voltage or temperature to
    dominate solely because of its units.
    """

    fields_to_check = (
        "temperature",
        "psi",
        "voltage",
        "electric_field_x",
        "electric_field_y",
        "current_density_x",
        "current_density_y",
        "supercurrent_density_x",
        "supercurrent_density_y",
        "normal_current_density_x",
        "normal_current_density_y",
        "heat_source",
    )

    residuals = []

    for name in fields_to_check:

        if name not in previous.arrays:
            continue

        if name not in candidate.arrays:
            continue

        residuals.append(
            relative_norm(
                candidate.arrays[name],
                previous.arrays[name],
            )
        )

    if not residuals:
        return math.inf

    value = max(residuals)

    if not np.isfinite(value):
        return math.inf

    return float(value)


# ============================================================================
# MODEL CREATION
# ============================================================================


def create_tdgl_model() -> TDGLModel:
    """
    Construct the current SHS TDGL model.

    TDGLParameters currently provides defaults for all parameters, so the
    benchmark deliberately uses the same default construction used by the
    existing full-coupling test.
    """

    parameters = TDGLParameters()

    return TDGLModel(
        parameters
    )


def create_thermal_model() -> ThermalModel:
    """
    Construct the thermal model using the same bath configuration as the
    existing full-coupling test.
    """

    return ThermalModel(
        bath_temperature=DEFAULT_BATH_TEMPERATURE,
        thermal_relaxation_rate=(
            DEFAULT_THERMAL_RELAXATION_RATE
        ),
    )


def create_electrical_model() -> ElectricalModel:
    """
    Construct the electrical model using its current defaults.
    """

    return ElectricalModel()


# ============================================================================
# PHYSICS OPERATIONS
# ============================================================================


def run_tdgl(
    simulation,
    dt: float,
    tdgl_model: TDGLModel,
    accepted_state: FieldSnapshot,
) -> None:
    """
    Advance TDGL for one physical timestep.

    The order parameter must begin this physical timestep from psi_n.

    Temperature and electromagnetic/electrical coupling quantities remain
    available from the current Picard guess.
    """

    fields = simulation.fields

    psi_guess = fields.psi.copy()

    fields.psi = (
        accepted_state.arrays["psi"].copy()
    )

    tdgl_step(
        simulation,
        dt,
        tdgl_model,
    )

    # tdgl_step computes the new psi and supercurrent.
    #
    # The local variable is retained above only to make the temporal-anchor
    # logic explicit.
    del psi_guess


def run_electrical(
    simulation,
    electrical_model: ElectricalModel,
    voltage_left: float,
    voltage_right: float,
) -> None:
    """
    Perform the electrical solve.

    The electrical solve is algebraic rather than a time integration, so
    the current Picard values of psi/supercurrent are used directly.
    """

    fields = simulation.fields

    superconducting_fraction = (
        np.abs(fields.psi) ** 2
    )

    electrical_step(
        fields=fields,
        mesh=simulation.mesh,
        material_map=simulation.material_map,
        contact_map=simulation.contact_map,
        voltage_left=voltage_left,
        voltage_right=voltage_right,
        superconducting_fraction=(
            superconducting_fraction
        ),
        superconducting_current_x=(
            fields.supercurrent_density_x
        ),
        superconducting_current_y=(
            fields.supercurrent_density_y
        ),
    )


def run_thermal(
    simulation,
    dt: float,
    thermal_model: ThermalModel,
    accepted_state: FieldSnapshot,
) -> None:
    """
    Advance thermal evolution for one physical timestep.

    Temperature is anchored to T_n for every Picard iteration.

    The current heat source is supplied by the current coupling guess.
    """

    fields = simulation.fields

    fields.temperature = (
        accepted_state.arrays[
            "temperature"
        ].copy()
    )

    thermal_step(
        simulation,
        dt,
        thermal_model,
    )


# ============================================================================
# PHYSICS ORDERING
# ============================================================================


def parse_ordering(
    ordering: str,
) -> Tuple[str, str, str]:
    """
    Validate and split a physics ordering.
    """

    parts = (
        ordering
        .lower()
        .split("_")
    )

    expected = {
        "tdgl",
        "electrical",
        "thermal",
    }

    if len(parts) != 3:
        raise ValueError(
            f"Invalid ordering '{ordering}'. "
            "Expected three physics names."
        )

    if set(parts) != expected:
        raise ValueError(
            f"Invalid ordering '{ordering}'. "
            "Expected exactly one each of "
            "tdgl, electrical, and thermal."
        )

    return (
        parts[0],
        parts[1],
        parts[2],
    )


def run_physics_operation(
    operation: str,
    simulation,
    dt: float,
    tdgl_model: TDGLModel,
    thermal_model: ThermalModel,
    electrical_model: ElectricalModel,
    accepted_state: FieldSnapshot,
    voltage_left: float,
    voltage_right: float,
) -> None:
    """
    Run one subsystem in the requested Picard ordering.
    """

    if operation == "tdgl":

        run_tdgl(
            simulation=simulation,
            dt=dt,
            tdgl_model=tdgl_model,
            accepted_state=accepted_state,
        )

    elif operation == "electrical":

        run_electrical(
            simulation=simulation,
            electrical_model=electrical_model,
            voltage_left=voltage_left,
            voltage_right=voltage_right,
        )

    elif operation == "thermal":

        run_thermal(
            simulation=simulation,
            dt=dt,
            thermal_model=thermal_model,
            accepted_state=accepted_state,
        )

    else:

        raise ValueError(
            f"Unknown physics operation: {operation}"
        )


# ============================================================================
# PICARD CANDIDATE
# ============================================================================


def calculate_candidate_state(
    simulation,
    accepted_state: FieldSnapshot,
    previous_guess: FieldSnapshot,
    ordering: Sequence[str],
    dt: float,
    tdgl_model: TDGLModel,
    thermal_model: ThermalModel,
    electrical_model: ElectricalModel,
    voltage_left: float,
    voltage_right: float,
) -> FieldSnapshot:
    """
    Calculate one candidate state for a fixed physical timestep.

    The important distinction here is:

        accepted_state
            = actual state at t_n

        previous_guess
            = current estimate of state at t_(n+1)

        candidate
            = result of applying one sequential Picard map to that guess

    TDGL and thermal evolution are anchored to the accepted t_n state each
    time they are called. Algebraic quantities remain available from the
    current coupling guess.

    This allows different orderings to represent different Gauss-Seidel-like
    fixed-point maps without accidentally advancing physical time multiple
    times.
    """

    fields = simulation.fields

    # ----------------------------------------------------------------------
    # Start from the current coupling guess.
    # ----------------------------------------------------------------------

    restore_fields(
        fields,
        previous_guess,
    )

    # ----------------------------------------------------------------------
    # Execute requested ordering.
    # ----------------------------------------------------------------------

    for operation in ordering:

        run_physics_operation(
            operation=operation,
            simulation=simulation,
            dt=dt,
            tdgl_model=tdgl_model,
            thermal_model=thermal_model,
            electrical_model=electrical_model,
            accepted_state=accepted_state,
            voltage_left=voltage_left,
            voltage_right=voltage_right,
        )

    # ----------------------------------------------------------------------
    # Candidate state represents one estimate of t_(n+1).
    # ----------------------------------------------------------------------

    return snapshot_fields(
        fields
    )


# ============================================================================
# SINGLE PHYSICAL TIMESTEP
# ============================================================================


def solve_coupled_timestep(
    simulation,
    dt: float,
    ordering: str,
    sigma: float,
    tolerance: float,
    max_iterations: int,
    tdgl_model: TDGLModel,
    thermal_model: ThermalModel,
    electrical_model: ElectricalModel,
    voltage_left: float,
    voltage_right: float,
) -> CouplingIterationResult:
    """
    Solve one physical timestep using Picard coupling.

    Physical time advances exactly once:

        t_n -> t_(n+1)

    regardless of how many coupling iterations are required.
    """

    if not (
        0.0 < sigma <= 1.0
    ):
        raise ValueError(
            "sigma must satisfy "
            "0 < sigma <= 1."
        )

    if tolerance <= 0.0:
        raise ValueError(
            "Coupling tolerance must be positive."
        )

    if max_iterations < 1:
        raise ValueError(
            "max_iterations must be >= 1."
        )

    parsed_ordering = parse_ordering(
        ordering
    )

    # ----------------------------------------------------------------------
    # Accepted state at t_n.
    # ----------------------------------------------------------------------

    accepted_state = snapshot_fields(
        simulation.fields
    )

    # ----------------------------------------------------------------------
    # Initial Picard guess.
    #
    # The natural initial guess is the previous physical state.
    # ----------------------------------------------------------------------

    previous_guess = copy_snapshot(
        accepted_state
    )

    residual_history = []

    # ----------------------------------------------------------------------
    # Picard iterations.
    # ----------------------------------------------------------------------

    for iteration in range(
        1,
        max_iterations + 1,
    ):

        candidate = calculate_candidate_state(
            simulation=simulation,
            accepted_state=accepted_state,
            previous_guess=previous_guess,
            ordering=parsed_ordering,
            dt=dt,
            tdgl_model=tdgl_model,
            thermal_model=thermal_model,
            electrical_model=electrical_model,
            voltage_left=voltage_left,
            voltage_right=voltage_right,
        )

        residual = state_residual(
            previous_guess,
            candidate,
        )

        residual_history.append(
            residual
        )

        # --------------------------------------------------------------
        # Numerical failure.
        # --------------------------------------------------------------

        if not np.isfinite(residual):

            restore_fields(
                simulation.fields,
                accepted_state,
            )

            return CouplingIterationResult(
                converged=False,
                iterations=iteration,
                final_residual=math.inf,
                residual_history=(
                    residual_history
                ),
                failed=True,
            )

        # --------------------------------------------------------------
        # Fixed-point convergence.
        # --------------------------------------------------------------

        if residual <= tolerance:

            restore_fields(
                simulation.fields,
                candidate,
            )

            return CouplingIterationResult(
                converged=True,
                iterations=iteration,
                final_residual=residual,
                residual_history=(
                    residual_history
                ),
                failed=False,
            )

        # --------------------------------------------------------------
        # Relax candidate to form next coupling guess.
        # --------------------------------------------------------------

        previous_guess = relaxed_snapshot(
            previous=previous_guess,
            candidate=candidate,
            sigma=sigma,
        )

    # ----------------------------------------------------------------------
    # Coupling failed to converge.
    #
    # Do not leave a partially converged state in the Simulation.
    # ----------------------------------------------------------------------

    restore_fields(
        simulation.fields,
        accepted_state,
    )

    return CouplingIterationResult(
        converged=False,
        iterations=max_iterations,
        final_residual=(
            residual_history[-1]
            if residual_history
            else math.inf
        ),
        residual_history=(
            residual_history
        ),
        failed=False,
    )


# ============================================================================
# FINAL OBSERVABLES
# ============================================================================


def finite_or_nan(
    value: float,
) -> float:
    """
    Return NaN instead of an invalid numerical value.
    """

    if np.isfinite(value):
        return float(value)

    return float("nan")


def final_observables(
    simulation,
) -> Dict[str, float]:
    """
    Extract useful physical quantities from the final state.
    """

    fields = simulation.fields

    temperature = np.asarray(
        fields.temperature
    )

    psi_amplitude = np.abs(
        fields.psi
    )

    current_x = np.asarray(
        fields.current_density_x
    )

    current_y = np.asarray(
        fields.current_density_y
    )

    current_magnitude = np.sqrt(
        current_x ** 2
        + current_y ** 2
    )

    heat_source = np.asarray(
        fields.heat_source
    )

    return {
        "final_temperature": finite_or_nan(
            np.mean(temperature)
        ),

        "final_temperature_max": finite_or_nan(
            np.max(temperature)
        ),

        "final_psi_amplitude_mean": (
            finite_or_nan(
                np.mean(psi_amplitude)
            )
        ),

        "final_psi_amplitude_min": (
            finite_or_nan(
                np.min(psi_amplitude)
            )
        ),

        "final_psi_amplitude_max": (
            finite_or_nan(
                np.max(psi_amplitude)
            )
        ),

        "final_current_mean": finite_or_nan(
            np.mean(current_magnitude)
        ),

        "final_current_max": finite_or_nan(
            np.max(current_magnitude)
        ),

        "final_joule_heating_mean": (
            finite_or_nan(
                np.mean(heat_source)
            )
        ),

        "final_joule_heating_max": (
            finite_or_nan(
                np.max(heat_source)
            )
        ),
    }


# ============================================================================
# SINGLE BENCHMARK CASE
# ============================================================================


def run_benchmark_case(
    config_path: str,
    ordering: str,
    sigma: float,
    dt: float,
    steps: int,
    tolerance: float,
    max_iterations: int,
    voltage_left: float,
    voltage_right: float,
) -> BenchmarkCaseResult:
    """
    Run one complete ordering / sigma benchmark case.

    A fresh Simulation is constructed for every case.
    """

    simulation = build_simulation(
        config_path
    )

    tdgl_model = create_tdgl_model()

    thermal_model = create_thermal_model()

    electrical_model = create_electrical_model()

    iteration_counts = []

    residual_history = []

    failed_steps = 0

    start_time = time.perf_counter()

    for step in range(
        steps
    ):

        result = solve_coupled_timestep(
            simulation=simulation,
            dt=dt,
            ordering=ordering,
            sigma=sigma,
            tolerance=tolerance,
            max_iterations=(
                max_iterations
            ),
            tdgl_model=tdgl_model,
            thermal_model=thermal_model,
            electrical_model=electrical_model,
            voltage_left=voltage_left,
            voltage_right=voltage_right,
        )

        iteration_counts.append(
            result.iterations
        )

        residual_history.extend(
            result.residual_history
        )

        if not result.converged:

            failed_steps += 1

            break

    total_runtime = (
        time.perf_counter()
        - start_time
    )

    successful_steps = (
        steps - failed_steps
    )

    if successful_steps > 0:

        runtime_per_step = (
            total_runtime
            / successful_steps
        )

    else:

        runtime_per_step = float("nan")

    converged = (
        failed_steps == 0
        and len(iteration_counts) == steps
    )

    if iteration_counts:

        average_iterations = float(
            np.mean(iteration_counts)
        )

        maximum_iterations = int(
            np.max(iteration_counts)
        )

    else:

        average_iterations = float("nan")

        maximum_iterations = 0

    if residual_history:

        final_residual = (
            residual_history[-1]
        )

    else:

        final_residual = math.inf

    observables = final_observables(
        simulation
    )

    return BenchmarkCaseResult(
        ordering=ordering,
        sigma=float(sigma),

        steps=steps,
        dt=float(dt),

        converged=converged,
        failed_steps=failed_steps,

        total_runtime=float(
            total_runtime
        ),

        runtime_per_step=float(
            runtime_per_step
        ),

        average_iterations=float(
            average_iterations
        ),

        maximum_iterations=(
            maximum_iterations
        ),

        final_residual=float(
            final_residual
        ),

        final_temperature=(
            observables[
                "final_temperature"
            ]
        ),

        final_temperature_max=(
            observables[
                "final_temperature_max"
            ]
        ),

        final_psi_amplitude_mean=(
            observables[
                "final_psi_amplitude_mean"
            ]
        ),

        final_psi_amplitude_min=(
            observables[
                "final_psi_amplitude_min"
            ]
        ),

        final_psi_amplitude_max=(
            observables[
                "final_psi_amplitude_max"
            ]
        ),

        final_current_mean=(
            observables[
                "final_current_mean"
            ]
        ),

        final_current_max=(
            observables[
                "final_current_max"
            ]
        ),

        final_joule_heating_mean=(
            observables[
                "final_joule_heating_mean"
            ]
        ),

        final_joule_heating_max=(
            observables[
                "final_joule_heating_max"
            ]
        ),

        residual_history=(
            residual_history
        ),
    )


# ============================================================================
# BENCHMARK SWEEP
# ============================================================================


def run_benchmark_sweep(
    config_path: str,
    orderings: Sequence[str],
    sigmas: Sequence[float],
    dt: float,
    steps: int,
    tolerance: float,
    max_iterations: int,
    voltage_left: float,
    voltage_right: float,
) -> List[BenchmarkCaseResult]:
    """
    Run the complete requested benchmark matrix.
    """

    results = []

    total_cases = (
        len(orderings)
        * len(sigmas)
    )

    case_number = 0

    for ordering in orderings:

        for sigma in sigmas:

            case_number += 1

            print()
            print("=" * 72)
            print(
                f"CASE {case_number}/{total_cases}"
            )
            print("=" * 72)

            print(
                f"Ordering : {ordering}"
            )

            print(
                f"Sigma    : {sigma:.3f}"
            )

            try:

                result = run_benchmark_case(
                    config_path=config_path,
                    ordering=ordering,
                    sigma=sigma,
                    dt=dt,
                    steps=steps,
                    tolerance=tolerance,
                    max_iterations=(
                        max_iterations
                    ),
                    voltage_left=(
                        voltage_left
                    ),
                    voltage_right=(
                        voltage_right
                    ),
                )

                results.append(
                    result
                )

                print_case_result(
                    result
                )

            except Exception as exc:

                print()
                print(
                    "CASE FAILED"
                )

                print(
                    f"{type(exc).__name__}: "
                    f"{exc}"
                )

                results.append(
                    BenchmarkCaseResult(
                        ordering=ordering,
                        sigma=float(sigma),
                        steps=steps,
                        dt=float(dt),
                        converged=False,
                        failed_steps=steps,
                        total_runtime=float(
                            "nan"
                        ),
                        runtime_per_step=float(
                            "nan"
                        ),
                        average_iterations=float(
                            "nan"
                        ),
                        maximum_iterations=0,
                        final_residual=float(
                            "nan"
                        ),
                        final_temperature=float(
                            "nan"
                        ),
                        final_temperature_max=float(
                            "nan"
                        ),
                        final_psi_amplitude_mean=float(
                            "nan"
                        ),
                        final_psi_amplitude_min=float(
                            "nan"
                        ),
                        final_psi_amplitude_max=float(
                            "nan"
                        ),
                        final_current_mean=float(
                            "nan"
                        ),
                        final_current_max=float(
                            "nan"
                        ),
                        final_joule_heating_mean=float(
                            "nan"
                        ),
                        final_joule_heating_max=float(
                            "nan"
                        ),
                        residual_history=[],
                    )
                )

    return results


# ============================================================================
# RESULT SERIALIZATION
# ============================================================================


def result_to_dict(
    result: BenchmarkCaseResult,
) -> Dict:
    """
    Convert one result to a serializable dictionary.
    """

    return {
        "ordering": result.ordering,
        "sigma": result.sigma,

        "steps": result.steps,
        "dt": result.dt,

        "converged": result.converged,
        "failed_steps": result.failed_steps,

        "total_runtime": (
            result.total_runtime
        ),

        "runtime_per_step": (
            result.runtime_per_step
        ),

        "average_iterations": (
            result.average_iterations
        ),

        "maximum_iterations": (
            result.maximum_iterations
        ),

        "final_residual": (
            result.final_residual
        ),

        "final_temperature": (
            result.final_temperature
        ),

        "final_temperature_max": (
            result.final_temperature_max
        ),

        "final_psi_amplitude_mean": (
            result.final_psi_amplitude_mean
        ),

        "final_psi_amplitude_min": (
            result.final_psi_amplitude_min
        ),

        "final_psi_amplitude_max": (
            result.final_psi_amplitude_max
        ),

        "final_current_mean": (
            result.final_current_mean
        ),

        "final_current_max": (
            result.final_current_max
        ),

        "final_joule_heating_mean": (
            result.final_joule_heating_mean
        ),

        "final_joule_heating_max": (
            result.final_joule_heating_max
        ),

        "residual_history": (
            result.residual_history
        ),
    }


def write_json(
    results: Sequence[BenchmarkCaseResult],
    path: Path,
    metadata: Dict,
) -> None:
    """
    Write complete benchmark results to JSON.
    """

    payload = {
        "metadata": metadata,

        "results": [
            result_to_dict(result)
            for result in results
        ],
    }

    with path.open(
        "w",
        encoding="utf-8",
    ) as handle:

        json.dump(
            payload,
            handle,
            indent=2,
        )


def write_csv(
    results: Sequence[BenchmarkCaseResult],
    path: Path,
) -> None:
    """
    Write summary results to CSV.

    Residual histories remain in JSON because they have variable length.
    """

    if not results:
        return

    rows = [
        result_to_dict(result)
        for result in results
    ]

    fieldnames = [
        key
        for key in rows[0]
        if key != "residual_history"
    ]

    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )

        writer.writeheader()

        for row in rows:

            writer.writerow(
                {
                    key: row[key]
                    for key in fieldnames
                }
            )


# ============================================================================
# TERMINAL REPORTING
# ============================================================================


def print_case_result(
    result: BenchmarkCaseResult,
) -> None:
    """
    Print one case's summary.
    """

    status = (
        "CONVERGED"
        if result.converged
        else "FAILED"
    )

    print()
    print(
        f"Status             : {status}"
    )

    print(
        f"Runtime            : "
        f"{result.total_runtime:.4f} s"
    )

    print(
        f"Runtime / step     : "
        f"{result.runtime_per_step:.4f} s"
    )

    print(
        f"Average iterations : "
        f"{result.average_iterations:.2f}"
    )

    print(
        f"Maximum iterations : "
        f"{result.maximum_iterations}"
    )

    print(
        f"Final residual     : "
        f"{result.final_residual:.6e}"
    )

    print(
        f"Failed steps       : "
        f"{result.failed_steps}"
    )


def print_summary(
    results: Sequence[BenchmarkCaseResult],
) -> None:
    """
    Print the complete benchmark comparison table.
    """

    print()
    print()
    print("=" * 116)
    print(
        "COUPLED CONVERGENCE BENCHMARK SUMMARY"
    )
    print("=" * 116)

    print(
        f"{'Ordering':<32}"
        f"{'Sigma':>8}"
        f"{'Status':>13}"
        f"{'Avg Iter':>12}"
        f"{'Max Iter':>12}"
        f"{'Runtime':>14}"
        f"{'Residual':>15}"
    )

    print("-" * 116)

    for result in results:

        status = (
            "CONVERGED"
            if result.converged
            else "FAILED"
        )

        print(
            f"{result.ordering:<32}"
            f"{result.sigma:>8.2f}"
            f"{status:>13}"
            f"{result.average_iterations:>12.2f}"
            f"{result.maximum_iterations:>12}"
            f"{result.total_runtime:>14.4f}"
            f"{result.final_residual:>15.3e}"
        )

    print("=" * 116)


def print_best_cases(
    results: Sequence[BenchmarkCaseResult],
) -> None:
    """
    Report useful numerical observations.

    This does not claim that the fastest case is physically best.
    """

    successful = [
        result
        for result in results
        if result.converged
        and np.isfinite(
            result.total_runtime
        )
    ]

    if not successful:

        print()
        print(
            "No successful benchmark cases."
        )

        return

    fastest = min(
        successful,
        key=lambda result:
            result.total_runtime,
    )

    fewest_iterations = min(
        successful,
        key=lambda result:
            result.average_iterations,
    )

    print()
    print(
        "BENCHMARK OBSERVATIONS"
    )
    print("-" * 72)

    print(
        "Fastest successful case:"
    )

    print(
        f"  {fastest.ordering}"
        f", sigma={fastest.sigma:.2f}"
        f", runtime={fastest.total_runtime:.4f} s"
    )

    print()
    print(
        "Fewest average coupling iterations:"
    )

    print(
        f"  {fewest_iterations.ordering}"
        f", sigma={fewest_iterations.sigma:.2f}"
        f", average="
        f"{fewest_iterations.average_iterations:.2f}"
    )

    print()
    print(
        "These observations characterize numerical "
        "performance only."
    )

    print(
        "They are not by themselves physical validation."
    )


# ============================================================================
# COMMAND LINE
# ============================================================================


def build_argument_parser():
    """
    Construct the command-line interface.
    """

    parser = argparse.ArgumentParser(
        description=(
            "Benchmark Picard convergence of the "
            "SHS TDGL/electrical/thermal system."
        )
    )

    parser.add_argument(
        "--config",
        default=DEFAULT_CONFIG,
        help=(
            "Simulation configuration file."
        ),
    )

    parser.add_argument(
        "--dt",
        type=float,
        default=DEFAULT_DT,
        help=(
            "Physical timestep."
        ),
    )

    parser.add_argument(
        "--steps",
        type=int,
        default=DEFAULT_STEPS,
        help=(
            "Number of physical timesteps per case."
        ),
    )

    parser.add_argument(
        "--tolerance",
        type=float,
        default=DEFAULT_COUPLING_TOLERANCE,
        help=(
            "Coupling convergence tolerance."
        ),
    )

    parser.add_argument(
        "--max-iterations",
        type=int,
        default=DEFAULT_MAX_COUPLING_ITERATIONS,
        help=(
            "Maximum Picard iterations per timestep."
        ),
    )

    parser.add_argument(
        "--sigma",
        type=float,
        nargs="+",
        default=None,
        help=(
            "Relaxation value(s). "
            "If omitted, use the complete default sweep."
        ),
    )

    parser.add_argument(
        "--ordering",
        nargs="+",
        default=None,
        help=(
            "Physics ordering(s). "
            "If omitted, use all six orderings."
        ),
    )

    parser.add_argument(
        "--voltage-left",
        type=float,
        default=DEFAULT_VOLTAGE_LEFT,
        help=(
            "Left contact voltage."
        ),
    )

    parser.add_argument(
        "--voltage-right",
        type=float,
        default=DEFAULT_VOLTAGE_RIGHT,
        help=(
            "Right contact voltage."
        ),
    )

    parser.add_argument(
        "--output-dir",
        default=DEFAULT_OUTPUT_DIRECTORY,
        help=(
            "Directory for benchmark results."
        ),
    )

    parser.add_argument(
        "--csv-name",
        default=DEFAULT_CSV_NAME,
        help=(
            "CSV result filename."
        ),
    )

    parser.add_argument(
        "--json-name",
        default=DEFAULT_JSON_NAME,
        help=(
            "JSON result filename."
        ),
    )

    parser.add_argument(
        "--no-output",
        action="store_true",
        help=(
            "Do not write CSV or JSON results."
        ),
    )

    return parser


# ============================================================================
# ARGUMENT VALIDATION
# ============================================================================


def validate_arguments(
    args,
) -> None:
    """
    Validate command-line configuration.
    """

    if args.dt <= 0.0:

        raise ValueError(
            "--dt must be positive."
        )

    if args.steps < 1:

        raise ValueError(
            "--steps must be >= 1."
        )

    if args.tolerance <= 0.0:

        raise ValueError(
            "--tolerance must be positive."
        )

    if args.max_iterations < 1:

        raise ValueError(
            "--max-iterations must be >= 1."
        )

    sigmas = (
        args.sigma
        if args.sigma is not None
        else DEFAULT_SIGMAS
    )

    for sigma in sigmas:

        if not (
            0.0 < sigma <= 1.0
        ):

            raise ValueError(
                f"Invalid sigma={sigma}. "
                "Require 0 < sigma <= 1."
            )

    orderings = (
        args.ordering
        if args.ordering is not None
        else DEFAULT_ORDERINGS
    )

    for ordering in orderings:

        parse_ordering(
            ordering
        )


# ============================================================================
# MAIN
# ============================================================================


def main() -> int:
    """
    Benchmark entry point.
    """

    parser = build_argument_parser()

    args = parser.parse_args()

    try:

        validate_arguments(
            args
        )

    except ValueError as exc:

        parser.error(
            str(exc)
        )

    config_path = Path(
        args.config
    )

    if not config_path.exists():

        print(
            "ERROR: configuration file does not exist:"
        )

        print(
            f"  {config_path}"
        )

        return 1

    sigmas = (
        list(args.sigma)
        if args.sigma is not None
        else list(DEFAULT_SIGMAS)
    )

    orderings = (
        list(args.ordering)
        if args.ordering is not None
        else list(DEFAULT_ORDERINGS)
    )

    print()
    print("=" * 72)
    print(
        "SUPERCONDUCTING HOTSPOT SIMULATOR"
    )
    print(
        "COUPLED CONVERGENCE BENCHMARK"
    )
    print("=" * 72)

    print(
        f"Configuration       : {config_path}"
    )

    print(
        f"Physical timestep   : "
        f"{args.dt:.6e}"
    )

    print(
        f"Physical steps      : "
        f"{args.steps}"
    )

    print(
        f"Coupling tolerance  : "
        f"{args.tolerance:.6e}"
    )

    print(
        f"Maximum iterations  : "
        f"{args.max_iterations}"
    )

    print(
        "Relaxation values   : "
        + ", ".join(
            f"{sigma:.2f}"
            for sigma in sigmas
        )
    )

    print(
        f"Physics orderings   : "
        f"{len(orderings)}"
    )

    print(
        "Total benchmark cases: "
        f"{len(orderings) * len(sigmas)}"
    )

    print("=" * 72)

    results = run_benchmark_sweep(
        config_path=str(
            config_path
        ),
        orderings=orderings,
        sigmas=sigmas,
        dt=args.dt,
        steps=args.steps,
        tolerance=args.tolerance,
        max_iterations=(
            args.max_iterations
        ),
        voltage_left=(
            args.voltage_left
        ),
        voltage_right=(
            args.voltage_right
        ),
    )

    print_summary(
        results
    )

    print_best_cases(
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

        csv_path = (
            output_directory
            / args.csv_name
        )

        json_path = (
            output_directory
            / args.json_name
        )

        metadata = {
            "configuration": str(
                config_path
            ),

            "dt": args.dt,

            "steps": args.steps,

            "coupling_tolerance": (
                args.tolerance
            ),

            "max_coupling_iterations": (
                args.max_iterations
            ),

            "sigmas": sigmas,

            "orderings": orderings,

            "voltage_left": (
                args.voltage_left
            ),

            "voltage_right": (
                args.voltage_right
            ),

            "bath_temperature": (
                DEFAULT_BATH_TEMPERATURE
            ),

            "thermal_relaxation_rate": (
                DEFAULT_THERMAL_RELAXATION_RATE
            ),
        }

        write_csv(
            results,
            csv_path,
        )

        write_json(
            results,
            json_path,
            metadata,
        )

        print()
        print(
            "Results written:"
        )

        print(
            f"  CSV : {csv_path}"
        )

        print(
            f"  JSON: {json_path}"
        )

    print()
    print(
        "Benchmark complete."
    )

    return 0


if __name__ == "__main__":
    sys.exit(
        main()
    )