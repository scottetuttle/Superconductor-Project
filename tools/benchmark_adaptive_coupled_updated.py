"""
Adaptive coupled convergence benchmark for SHS.

This benchmark studies convergence of the coupled TDGL/electrical/thermal
system and experimentally evaluates automatic numerical adaptation.

There are two distinct levels of iteration:

    Physical time:
        t_n -> t_(n+1)

    Coupling iterations:
        k = 0, 1, 2, ...

Coupling iterations do NOT advance physical time.

At a physical timestep, the accepted state at t_n remains fixed while
the coupled state at t_(n+1) is iteratively solved.

If convergence is predicted to be poor, the benchmark may automatically
adapt numerical parameters and restart the SAME physical timestep.

Adaptation hierarchy:

    1. Oscillation
        -> reduce coupling relaxation sigma.

    2. One problematic component
        -> adapt only that component.

    3. Two problematic components
        -> adapt both components.

    4. Three problematic components
        -> reduce physical dt and retry the timestep.

    5. Ineffective component adaptation
        -> fall back to reducing the coupled physical dt.

    If physical dt is reduced, the original requested physical-step duration is
    completed using additional accepted substeps so adaptive dt does not alter
    the requested total simulated time.

The accepted physical state is never overwritten by a failed attempt.

This benchmark intentionally keeps the adaptive controller separate from
the production CoupledSolver. Its purpose is to experimentally evaluate
numerical adaptation strategies before transferring them into production.
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
import math
import time
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from shs.config.builder import build_simulation

from shs.solvers.tdgl_solver import tdgl_step
from shs.solvers.electrical_solver import electrical_step
from shs.solvers.thermal_solver import thermal_step

from shs.physics.thermal import ThermalModel

from shs.tdgl import TDGLModel
from shs.tdgl.parameters import TDGLParameters

from shs.numerics.convergence import (
    ConvergenceController,
    ConvergenceMonitor,
    ConvergenceStatus,
)


CONFIG_PATH = Path(
    "tools/config/coupled_convergence_benchmark.json"
)


DEFAULT_ORDERINGS = (
    "tdgl_electrical_thermal",
    "tdgl_thermal_electrical",
    "electrical_tdgl_thermal",
    "electrical_thermal_tdgl",
    "thermal_tdgl_electrical",
    "thermal_electrical_tdgl",
)


# ============================================================================
# Data structures
# ============================================================================


@dataclass
class TrialState:
    temperature: np.ndarray
    psi: np.ndarray
    voltage: np.ndarray
    vector_potential_x: np.ndarray
    vector_potential_y: np.ndarray


@dataclass
class NumericalParameters:
    physical_dt: float
    sigma: float
    electrical_tolerance: float
    tdgl_max_normalized_timestep: float
    thermal_max_substep: float


@dataclass
class AdaptationEvent:
    physical_step: int
    attempt: int
    iteration: int
    reason: str

    parameter: str
    old_value: float
    new_value: float

    residual: float
    status: str

    predicted_iterations: float | None
    prediction_confidence: float

    problematic_components: list[str] = field(
        default_factory=list
    )

    successful_after_restart: bool | None = None


@dataclass
class CouplingResult:
    converged: bool
    iterations: int
    final_residual: float

    residual_history: list[float]

    temperature_residual: list[float]
    psi_residual: list[float]
    voltage_residual: list[float]
    ax_residual: list[float]
    ay_residual: list[float]

    sigma_history: list[float]
    oscillation_history: list[bool]

    residual_ratio_history: list[float | None]
    log_slope_history: list[float | None]
    estimated_remaining_history: list[float | None]
    confidence_history: list[float]
    status_history: list[str]

    checkpoint_history: list[bool]
    adaptation_history: list[bool]

    actual_remaining_history: list[int | None]

    final_sigma: float

    prediction_errors: list[float]

    problematic_components_history: list[list[str]]

    timing_history: dict[str, list[float]]


# ============================================================================
# Basic utilities
# ============================================================================


def relative_residual(
    old,
    new,
    floor=1e-14,
):
    old = np.asarray(old)
    new = np.asarray(new)

    numerator = np.linalg.norm(
        (new - old).ravel()
    )

    denominator = max(
        np.linalg.norm(old.ravel()),
        floor,
    )

    return float(
        numerator / denominator
    )


def copy_trial_state(
    simulation,
):
    fields = simulation.fields

    return TrialState(
        temperature=fields.temperature.copy(),
        psi=fields.psi.copy(),
        voltage=fields.voltage.copy(),
        vector_potential_x=(
            fields.vector_potential_x.copy()
        ),
        vector_potential_y=(
            fields.vector_potential_y.copy()
        ),
    )


def apply_trial_state(
    simulation,
    trial,
):
    fields = simulation.fields

    fields.temperature = trial.temperature.copy()
    fields.psi = trial.psi.copy()
    fields.voltage = trial.voltage.copy()

    fields.vector_potential_x = (
        trial.vector_potential_x.copy()
    )

    fields.vector_potential_y = (
        trial.vector_potential_y.copy()
    )


def calculate_residuals(
    old,
    new,
):
    temperature = relative_residual(
        old.temperature,
        new.temperature,
    )

    psi = relative_residual(
        old.psi,
        new.psi,
    )

    voltage = relative_residual(
        old.voltage,
        new.voltage,
    )

    ax = relative_residual(
        old.vector_potential_x,
        new.vector_potential_x,
    )

    ay = relative_residual(
        old.vector_potential_y,
        new.vector_potential_y,
    )

    total = max(
        temperature,
        psi,
        voltage,
        ax,
        ay,
    )

    return {
        "temperature": temperature,
        "psi": psi,
        "voltage": voltage,
        "ax": ax,
        "ay": ay,
        "total": total,
    }


def relax_trial_state(
    old_trial,
    calculated,
    sigma,
):
    return TrialState(
        temperature=(
            old_trial.temperature
            + sigma
            * (
                calculated.temperature
                - old_trial.temperature
            )
        ),

        psi=(
            old_trial.psi
            + sigma
            * (
                calculated.psi
                - old_trial.psi
            )
        ),

        voltage=(
            old_trial.voltage
            + sigma
            * (
                calculated.voltage
                - old_trial.voltage
            )
        ),

        vector_potential_x=(
            old_trial.vector_potential_x
            + sigma
            * (
                calculated.vector_potential_x
                - old_trial.vector_potential_x
            )
        ),

        vector_potential_y=(
            old_trial.vector_potential_y
            + sigma
            * (
                calculated.vector_potential_y
                - old_trial.vector_potential_y
            )
        ),
    )


def state_update_directions(
    old,
    new,
):
    return {
        "temperature": (
            new.temperature
            - old.temperature
        ).ravel(),

        "psi": (
            new.psi
            - old.psi
        ).view(np.float64).ravel(),

        "voltage": (
            new.voltage
            - old.voltage
        ).ravel(),

        "ax": (
            new.vector_potential_x
            - old.vector_potential_x
        ).ravel(),

        "ay": (
            new.vector_potential_y
            - old.vector_potential_y
        ).ravel(),
    }


def detect_oscillation(
    previous_delta,
    current_delta,
):
    for name in previous_delta:

        old = previous_delta[name]
        new = current_delta[name]

        old_norm = np.linalg.norm(old)
        new_norm = np.linalg.norm(new)

        if (
            old_norm < 1e-14
            or new_norm < 1e-14
        ):
            continue

        dot = np.vdot(
            old,
            new,
        ).real

        if dot < 0.0:
            return True

    return False


# ============================================================================
# Physics operations
# ============================================================================


def run_tdgl(
    simulation,
    dt,
    tdgl_model,
):
    start = time.time()
    tdgl_step(
        simulation,
        dt,
        tdgl_model,
    )
    end = time.time()
    tdgl_step_time = end - start
    return tdgl_step_time


def run_electrical(
    simulation,
    electrical_voltage_left,
    electrical_voltage_right,
    electrical_tolerance,
):
    start = time.time()
    fields = simulation.fields

    electrical_step(
        fields,
        simulation.mesh,
        simulation.material_map,
        simulation.contact_map,
        voltage_left=electrical_voltage_left,
        voltage_right=electrical_voltage_right,
        superconducting_fraction=(
            np.abs(fields.psi) ** 2
        ),
        superconducting_current_x=(
            fields.supercurrent_density_x
        ),
        superconducting_current_y=(
            fields.supercurrent_density_y
        ),
        solver_tolerance=electrical_tolerance,
    )
    end = time.time()
    electrical_step_time = (end - start)
    return electrical_step_time

def run_thermal(
    simulation,
    dt,
    thermal_model,
    external_heat=None,
):
    start = time.time()
    if external_heat is not None:

        original_heat = (
            simulation.fields.heat_source.copy()
        )

        simulation.fields.heat_source = (
            original_heat
            + external_heat
        )

        thermal_step(
            simulation,
            dt,
            thermal_model,
        )

        simulation.fields.heat_source = (
            original_heat
        )

    else:

        thermal_step(
            simulation,
            dt,
            thermal_model,
        )
    end = time.time()
    thermal_step_time = end -start
    return thermal_step_time



def execute_ordering(
    simulation,
    ordering,
    dt,
    tdgl_model,
    thermal_model,
    voltage_left,
    voltage_right,
    external_heat,
):
    operations = ordering.split("_")

    timings = {
        "tdgl": 0.0,
        "electrical": 0.0,
        "thermal": 0.0
    }

    for operation in operations:

        if operation == "tdgl":

            timings["tdgl"] += run_tdgl(
                simulation,
                dt,
                tdgl_model,
            )

        elif operation == "electrical":

            timings["electrical"] += run_electrical(
                simulation,
                voltage_left,
                voltage_right,
                simulation.config.electrical.solver.tolerance
            )

        elif operation == "thermal":

            timings["thermal"] += run_thermal(
                simulation,
                dt,
                thermal_model,
                external_heat=external_heat,
            )

        else:

            raise ValueError(
                f"Unknown coupling operation: "
                f"{operation}"
            )
    return timings

# ============================================================================
# Convergence controller
# ============================================================================


def build_controller(
    tolerance,
    max_iterations,
    check_interval,
    prediction_window,
    reasonable_iterations,
):
    monitor = ConvergenceMonitor(
        tolerance=tolerance,
        prediction_window=prediction_window,
        reasonable_iterations=(
            reasonable_iterations
        ),
    )

    return ConvergenceController(
        monitor=monitor,
        max_iterations=max_iterations,
        check_interval=check_interval,
    )


# ============================================================================
# Component classification
# ============================================================================


def classify_component(
    history,
    tolerance,
    minimum_confidence,
    reasonable_iterations,
    stall_ratio,
    divergence_ratio,
):
    if len(history) < 2:
        return False

    residual = history[-1]

    if residual <= tolerance:
        return False

    previous = history[-2]

    if previous <= 0.0:
        return False

    ratio = residual / previous

    if ratio >= divergence_ratio:
        return True

    if ratio >= stall_ratio:
        return True

    if len(history) >= 3:

        count = min(
            len(history),
            10,
        )

        values = np.asarray(
            history[-count:],
            dtype=float,
        )

        if np.all(values > 0.0):

            iterations = np.arange(
                count,
                dtype=float,
            )

            slope, _ = np.polyfit(
                iterations,
                np.log10(values),
                1,
            )

            if slope >= 0.0:
                return True

            remaining = (
                (
                    np.log10(tolerance)
                    - np.log10(residual)
                )
                / slope
            )

            if (
                np.isfinite(remaining)
                and remaining
                > reasonable_iterations
            ):
                return True

            if (
                not np.isfinite(remaining)
                or remaining < 0.0
            ):
                return True

            predicted = np.polyval(
                np.polyfit(
                    iterations,
                    np.log10(values),
                    1,
                ),
                iterations,
            )

            variance = np.sum(
                (
                    np.log10(values)
                    - np.mean(
                        np.log10(values)
                    )
                )
                ** 2
            )

            if variance > 0.0:

                error = np.sum(
                    (
                        np.log10(values)
                        - predicted
                    )
                    ** 2
                )

                confidence = np.clip(
                    1.0
                    - error / variance,
                    0.0,
                    1.0,
                )
                if (
                    confidence
                    < minimum_confidence
                    and ratio
                    >= stall_ratio
                ):
                    return True
    return False


def identify_problematic_components(
    residual_histories,
    tolerance,
    controller_config,
):
    component_config = controller_config[
        "component_prediction"
    ]

    problematic = []

    for name, history in residual_histories.items():

        if classify_component(
            history=history,
            tolerance=tolerance,
            minimum_confidence=(
                float(
                    component_config[
                        "minimum_confidence"
                    ]
                )
            ),
            reasonable_iterations=(
                int(
                    component_config[
                        "reasonable_iterations"
                    ]
                )
            ),
            stall_ratio=(
                float(
                    component_config[
                        "stall_ratio"
                    ]
                )
            ),
            divergence_ratio=(
                float(
                    component_config[
                        "divergence_ratio"
                    ]
                )
            ),
        ):
            problematic.append(name)

    return problematic


# ============================================================================
# Runtime numerical parameter handling
# ============================================================================


def get_nested_attribute(
    obj,
    path,
    default=None,
):
    current = obj

    for name in path.split("."):

        if current is None:
            return default

        if not hasattr(
            current,
            name,
        ):
            return default

        current = getattr(
            current,
            name,
        )

    return current


def set_nested_attribute(
    obj,
    path,
    value,
):
    parts = path.split(".")

    current = obj

    for name in parts[:-1]:

        if not hasattr(
            current,
            name,
        ):
            return False

        current = getattr(
            current,
            name,
        )

        if current is None:
            return False

    final_name = parts[-1]

    if not hasattr(
        current,
        final_name,
    ):
        return False

    setattr(
        current,
        final_name,
        value,
    )

    return True


def apply_numerical_parameters(
    simulation,
    parameters,
):
    """
    Push the benchmark's numerical state into the simulation configuration.

    TDGL is constructed directly from parameters elsewhere.

    Electrical and thermal solver settings are applied to the simulation
    configuration when those configuration objects expose the expected
    attributes.
    """

    set_nested_attribute(
        simulation,
        "config.electrical.solver.tolerance",
        parameters.electrical_tolerance,
    )

    set_nested_attribute(
        simulation,
        "config.thermal.max_substep",
        parameters.thermal_max_substep,
    )


def build_tdgl_model(
    parameters,
):
    return TDGLModel(
        TDGLParameters(
            max_normalized_timestep=(
                parameters
                .tdgl_max_normalized_timestep
            )
        )
    )


# ============================================================================
# Adaptation
# ============================================================================


def clamp(
    value,
    minimum,
    maximum,
):
    return min(
        max(
            value,
            minimum,
        ),
        maximum,
    )


def adapt_parameter(
    parameters,
    parameter,
    factor,
    minimum,
    maximum,
):
    old_value = getattr(
        parameters,
        parameter,
    )

    new_value = clamp(
        old_value * factor,
        minimum,
        maximum,
    )

    changed = (
        not math.isclose(
            old_value,
            new_value,
            rel_tol=1e-14,
            abs_tol=0.0,
        )
    )

    if changed:
        setattr(
            parameters,
            parameter,
            new_value,
        )

    return (
        changed,
        old_value,
        new_value,
    )


def choose_physical_dt_fallback(
    parameters,
    controller_config,
    bounds,
):
    """Return a bounded physical-dt reduction action."""
    adaptation = controller_config["adaptation"]

    factor = float(adaptation["physical_dt_factor"])
    old_dt = float(parameters.physical_dt)
    minimum = float(bounds["physical_dt_min"])

    if factor >= 1.0:
        return []

    new_dt = clamp(
        old_dt * factor,
        minimum,
        old_dt,
    )

    if new_dt >= old_dt or math.isclose(
        new_dt,
        old_dt,
        rel_tol=1e-14,
        abs_tol=0.0,
    ):
        return []

    print(
        f"physical timestep reduction: "
        f"{old_dt:.3e} - {new_dt:.3e}"
    )
    return [
        (
            "physical_dt",
            old_dt,
            new_dt,
        )
    ]


def choose_adaptation(
    problematic_components,
    parameters,
    controller_config,
    bounds,
):
    """
    Choose the least aggressive adaptation appropriate to the failed attempt.

    The existing hierarchy is preserved:

        0 components -> no component-specific adaptation
        1 component  -> adapt that component
        2 components -> adapt both components
        3+ components -> reduce the coupled physical timestep

    A separate physical-dt fallback is used by the retry loop when a
    component-specific adaptation fails to produce meaningful improvement.
    """
    adaptation = controller_config["adaptation"]
    actions = []

    if len(problematic_components) == 0:
        return actions

    if len(problematic_components) >= 3:
        return choose_physical_dt_fallback(
            parameters,
            controller_config,
            bounds,
        )

    mapping = {
        "voltage": (
            "electrical_tolerance",
            float(adaptation["electrical_tolerance_factor"]),
            float(bounds["electrical_tolerance_min"]),
            float(bounds["electrical_tolerance_max"]),
        ),
        "psi": (
            "tdgl_max_normalized_timestep",
            float(adaptation["tdgl_max_normalized_timestep_factor"]),
            float(bounds["tdgl_max_normalized_timestep_min"]),
            float(bounds["tdgl_max_normalized_timestep_max"]),
        ),
        "temperature": (
            "thermal_max_substep",
            float(adaptation["thermal_max_substep_factor"]),
            float(bounds["thermal_max_substep_min"]),
            float(bounds["thermal_max_substep_max"]),
        ),
    }

    for component in problematic_components:
        if component not in mapping:
            continue

        parameter, factor, minimum, maximum = mapping[component]
        old_value = float(getattr(parameters, parameter))
        new_value = clamp(
            old_value * factor,
            minimum,
            maximum,
        )

        if not math.isclose(
            old_value,
            new_value,
            rel_tol=1e-14,
            abs_tol=0.0,
        ):
            actions.append(
                (
                    parameter,
                    old_value,
                    new_value,
                )
            )
    print(
        f"attempted adjustment"
        f"{actions}"
        )
    return actions


def adaptation_was_effective(
    previous_result,
    current_result,
    minimum_residual_improvement=0.10,
    minimum_iteration_improvement=0.10,
):
    """
    Decide whether an adaptation produced a meaningful improvement.

    An adaptation is considered effective when either the final residual or
    the number of coupling iterations improves by at least the configured
    threshold. This deliberately uses the previous and current failed
    attempts rather than merely checking whether a parameter changed.
    """
    if previous_result is None:
        return True

    previous_residual = float(previous_result.final_residual)
    current_residual = float(current_result.final_residual)

    if not np.isfinite(previous_residual):
        return np.isfinite(current_residual)

    if previous_residual <= 0.0:
        return True

    residual_improvement = (
        previous_residual - current_residual
    ) / previous_residual

    previous_iterations = max(
        int(previous_result.iterations),
        1,
    )
    current_iterations = int(current_result.iterations)

    iteration_improvement = (
        previous_iterations - current_iterations
    ) / previous_iterations

    return bool(
        residual_improvement >= minimum_residual_improvement
        or iteration_improvement >= minimum_iteration_improvement
    )


# ============================================================================
# One coupling attempt
# ============================================================================


def solve_coupled_attempt(
    simulation,
    accepted_state,
    parameters,
    ordering,
    thermal_model,
    voltage_left,
    voltage_right,
    external_heat,
    tolerance,
    max_iterations,
    check_interval,
    prediction_window,
    reasonable_iterations,
    adaptive_relaxation,
    minimum_sigma,
    reduction_factor,
    required_reversals,
    adaptive_controller_config,
):
    trial_state = copy.deepcopy(
        accepted_state
    )

    residual_history = []

    temperature_history = []
    psi_history = []
    voltage_history = []
    ax_history = []
    ay_history = []

    sigma_history = []
    oscillation_history = []

    residual_ratio_history = []
    log_slope_history = []
    estimated_remaining_history = []
    confidence_history = []

    status_history = []

    checkpoint_history = []
    adaptation_history = []

    actual_remaining_history = []

    problematic_components_history = []
    
    sigma = parameters.sigma

    previous_delta = None
    consecutive_reversals = 0

    timing_totals = {
        "tdgl":0.0,
        "electrical": 0.0,
        "thermal": 0.0
    }
    timing_history = {
        "tdgl": [],
        "electrical": [],
        "thermal": []
    }

    controller = build_controller(
        tolerance=tolerance,
        max_iterations=max_iterations,
        check_interval=check_interval,
        prediction_window=prediction_window,
        reasonable_iterations=(
            reasonable_iterations
        ),
    )

    prediction_errors = []

    residual_histories = {
        "temperature": [],
        "psi": [],
        "voltage": [],
    }

    converged = False

    for iteration in range(
        1,
        max_iterations + 1,
    ):

        apply_trial_state(
            simulation,
            accepted_state,
        )

        simulation.fields.temperature = (
            trial_state.temperature.copy()
        )

        simulation.fields.voltage = (
            trial_state.voltage.copy()
        )

        simulation.fields.vector_potential_x = (
            trial_state.vector_potential_x.copy()
        )

        simulation.fields.vector_potential_y = (
            trial_state.vector_potential_y.copy()
        )

        simulation.fields.psi = (
            accepted_state.psi.copy()
        )

        apply_numerical_parameters(
            simulation,
            parameters,
        )

        thermal_model.max_substep = (
            parameters.thermal_max_substep
        )

        tdgl_model = build_tdgl_model(
            parameters
        )

        iteration_timings = execute_ordering(
            simulation=simulation,
            ordering=ordering,
            dt=parameters.physical_dt,
            tdgl_model=tdgl_model,
            thermal_model=thermal_model,
            voltage_left=voltage_left,
            voltage_right=voltage_right,
            external_heat=external_heat,
        )

        for component, elapsed in iteration_timings.items():
            timing_history[component].append(elapsed)

        calculated = copy_trial_state(
            simulation
        )

        residuals = calculate_residuals(
            trial_state,
            calculated,
        )

        residual = residuals["total"]

        residual_history.append(
            residual
        )

        temperature_history.append(
            residuals["temperature"]
        )

        psi_history.append(
            residuals["psi"]
        )

        voltage_history.append(
            residuals["voltage"]
        )

        ax_history.append(
            residuals["ax"]
        )

        ay_history.append(
            residuals["ay"]
        )

        residual_histories[
            "temperature"
        ].append(
            residuals["temperature"]
        )

        residual_histories[
            "psi"
        ].append(
            residuals["psi"]
        )

        residual_histories[
            "voltage"
        ].append(
            residuals["voltage"]
        )

        sigma_history.append(
            sigma
        )

        assessment = controller.record(
            iteration=iteration,
            residual=residual,
        )

        residual_ratio_history.append(
            assessment.residual_ratio
        )

        log_slope_history.append(
            assessment.log_slope
        )

        estimated_remaining_history.append(
            assessment.estimated_iterations_remaining
        )

        confidence_history.append(
            assessment.prediction_confidence
        )

        status_history.append(
            assessment.status.value
        )

        checkpoint = controller.should_check(
            iteration
        )

        checkpoint_history.append(
            checkpoint
        )

        actual_remaining_history.append(
            None
        )

        current_delta = (
            state_update_directions(
                trial_state,
                calculated,
            )
        )

        oscillating = False

        if previous_delta is not None:

            oscillating = detect_oscillation(
                previous_delta,
                current_delta,
            )

        oscillation_history.append(
            oscillating
        )

        if oscillating:
            consecutive_reversals += 1
        else:
            consecutive_reversals = 0

        problematic = []

        if checkpoint:

            problematic = (
                identify_problematic_components(
                    residual_histories=(
                        residual_histories
                    ),
                    tolerance=tolerance,
                    controller_config=(
                        adaptive_controller_config
                    ),
                )
            )
        problematic_components_history.append(
            problematic.copy()
        )
        adaptation_triggered = False

        if (
            oscillating
            and adaptive_relaxation
            and adaptive_controller_config[
                "behavior"
            ][
                "reduce_sigma_only_for_oscillation"
            ]
            and consecutive_reversals
            >= required_reversals
        ):

            new_sigma = max(
                minimum_sigma,
                sigma * reduction_factor,
            )

            if new_sigma < sigma:

                sigma = new_sigma

                parameters.sigma = sigma

                adaptation_triggered = True

            consecutive_reversals = 0

        if adaptation_triggered:
            adaptation_history.append(True)
        else:
            adaptation_history.append(False)

        if (
            assessment.status
            == ConvergenceStatus.CONVERGED
        ):

            converged = True

            trial_state = calculated

            convergence_iteration = iteration

            for index in range(
                len(actual_remaining_history)
            ):

                iteration_number = index + 1

                actual = max(
                    convergence_iteration
                    - iteration_number,
                    0,
                )

                actual_remaining_history[
                    index
                ] = actual

                prediction = (
                    estimated_remaining_history[
                        index
                    ]
                )

                if (
                    prediction is not None
                    and actual > 0
                ):

                    prediction_errors.append(
                        float(
                            abs(
                                prediction
                                - actual
                            )
                            / actual
                        )
                    )

            break

        trial_state = relax_trial_state(
            trial_state,
            calculated,
            sigma,
        )

        previous_delta = current_delta



    apply_trial_state(
        simulation,
        trial_state,
    )

    final_residual = (
        residual_history[-1]
        if residual_history
        else float("inf")
    )
    

    return CouplingResult(
        converged=converged,
        iterations=len(
            residual_history
        ),
        final_residual=final_residual,

        residual_history=residual_history,

        temperature_residual=(
            temperature_history
        ),

        psi_residual=psi_history,

        voltage_residual=voltage_history,

        ax_residual=ax_history,

        ay_residual=ay_history,

        sigma_history=sigma_history,

        oscillation_history=(
            oscillation_history
        ),

        residual_ratio_history=(
            residual_ratio_history
        ),

        log_slope_history=(
            log_slope_history
        ),

        estimated_remaining_history=(
            estimated_remaining_history
        ),

        confidence_history=(
            confidence_history
        ),

        status_history=status_history,

        checkpoint_history=(
            checkpoint_history
        ),

        adaptation_history=(
            adaptation_history
        ),

        actual_remaining_history=(
            actual_remaining_history
        ),

        final_sigma=sigma,

        prediction_errors=prediction_errors,

        problematic_components_history=(
            problematic_components_history
        ),
        timing_history=(timing_history),

    )

# ============================================================================
# External heat
# ============================================================================


def make_gaussian_heat_source(
    simulation,
    amplitude,
    radius_cells,
):
    mesh = simulation.mesh

    y, x = np.indices(
        (
            mesh.ny,
            mesh.nx,
        )
    )

    cx = (
        mesh.nx - radius_cells
    ) / 2.0

    cy = (
        mesh.ny - radius_cells
    ) / 2.0

    r2 = (
        (x - cx) ** 2
        + (y - cy) ** 2
    )

    return amplitude * np.exp(
        -r2
        / (
            2.0
            * radius_cells ** 2
        )
    )


# ============================================================================
# Attempt history
# ============================================================================


def build_attempt_history(
    result,
    physical_step,
    attempt,
    parameters,
    restart_reason=None,
):
    """
    Convert one coupling attempt into iteration-level diagnostic history.

    Each attempt is stored independently so that failed attempts that trigger
    numerical adaptation or timestep reduction are preserved for analysis.
    """

    history = []

    for index, residual in enumerate(
        result.residual_history
    ):

        actual_remaining = (
            result.actual_remaining_history[index]
        )

        predicted_remaining = (
            result.estimated_remaining_history[index]
        )

        prediction_error = None

        if (
            predicted_remaining is not None
            and actual_remaining is not None
            and actual_remaining > 0
        ):
            prediction_error = (
                abs(
                    predicted_remaining
                    - actual_remaining
                )
                / actual_remaining
            )

        history.append(
            {
                "physical_step": physical_step,

                "attempt": attempt,

                "iteration": index + 1,

                "residual": residual,

                "temperature_residual": (
                    result.temperature_residual[index]
                ),

                "psi_residual": (
                    result.psi_residual[index]
                ),

                "voltage_residual": (
                    result.voltage_residual[index]
                ),

                "ax_residual": (
                    result.ax_residual[index]
                ),

                "ay_residual": (
                    result.ay_residual[index]
                ),

                "sigma": (
                    result.sigma_history[index]
                ),

                "physical_dt": (
                    parameters.physical_dt
                ),

                "electrical_tolerance": (
                    parameters.electrical_tolerance
                ),

                "tdgl_max_normalized_timestep": (
                    parameters.tdgl_max_normalized_timestep
                ),

                "thermal_max_substep": (
                    parameters.thermal_max_substep
                ),

                "oscillation": (
                    result.oscillation_history[index]
                ),

                "residual_ratio": (
                    result.residual_ratio_history[index]
                ),

                "log_slope": (
                    result.log_slope_history[index]
                ),

                "estimated_iterations_remaining": (
                    predicted_remaining
                ),

                "actual_iterations_remaining": (
                    actual_remaining
                ),

                "prediction_error": (
                    prediction_error
                ),

                "prediction_confidence": (
                    result.confidence_history[index]
                ),

                "status": (
                    result.status_history[index]
                ),

                "checkpoint": (
                    result.checkpoint_history[index]
                ),

                "adaptation_recommended": (
                    result.adaptation_history[index]
                ),

                "restart_reason": restart_reason,
            }
        )

    return history


def build_attempt_summary(attempt):
    """Return only start/end iteration information for compact JSON output."""
    history = attempt.get("history", [])

    if not history:
        return {
            "physical_step": attempt["physical_step"],
            "attempt": attempt["attempt"],
            "converged": attempt["converged"],
            "restart_reason": attempt["restart_reason"],
            "start": None,
            "end": None,
        }

    fields = (
        "iteration",
        "residual",
        "temperature_residual",
        "psi_residual",
        "voltage_residual",
        "ax_residual",
        "ay_residual",
        "sigma",
        "physical_dt",
        "electrical_tolerance",
        "tdgl_max_normalized_timestep",
        "thermal_max_substep",
        "oscillation",
        "status",
    )

    def select(row):
        return {
            key: row.get(key)
            for key in fields
        }

    return {
        "physical_step": attempt["physical_step"],
        "attempt": attempt["attempt"],
        "converged": attempt["converged"],
        "restart_reason": attempt["restart_reason"],
        "start": select(history[0]),
        "end": select(history[-1]),
    }


# ============================================================================
# Benchmark run
# ============================================================================


def run_single(
    config,
    dt,
    sigma,
    ordering,
):
    simulation = build_simulation(
        config["simulation_config"]
    )

    if config["initial_conditions"].get(
        "uniform_psi_one",
        False,
    ):
        simulation.fields.psi[:] = 1.0 + 0.0j

    thermal_model = ThermalModel(
        bath_temperature=config["thermal"]["bath_temperature"],
        thermal_relaxation_rate=config["thermal"]["thermal_relaxation_rate"],
    )

    voltage_left = config["electrical"]["voltage_left"]
    voltage_right = config["electrical"]["voltage_right"]

    external_heat = None
    if config["external_heat"]["enabled"]:
        external_heat = make_gaussian_heat_source(
            simulation,
            amplitude=config["external_heat"]["amplitude"],
            radius_cells=config["external_heat"]["radius_cells"],
        )

    physical_steps = int(config["physical_steps"])
    convergence_config = config["convergence"]
    max_iterations = int(convergence_config["max_iterations"])
    tolerance = float(convergence_config["tolerance"])
    check_interval = int(convergence_config["check_interval"])
    prediction_window = int(convergence_config["prediction_window"])
    reasonable_iterations = int(convergence_config["reasonable_iterations"])

    relaxation_config = config["adaptive_relaxation"]
    adaptive_relaxation = bool(relaxation_config["enabled"])
    minimum_sigma = float(relaxation_config["minimum_sigma"])
    reduction_factor = float(relaxation_config["reduction_factor"])
    required_reversals = int(relaxation_config["required_reversals"])

    controller_config = config["adaptive_controller"]
    bounds = controller_config["bounds"]

    simulation_thermal_max_substep = get_nested_attribute(
        simulation,
        "config.thermal.max_substep",
        1e-13,
    )
    simulation_electrical_tolerance = get_nested_attribute(
        simulation,
        "config.electrical.solver.tolerance",
        1e-10,
    )
    tdgl_max_normalized_timestep = config["initial_conditions"]["initial_tdgl_max_normalized_timestep"]
    print(tdgl_max_normalized_timestep)

    parameters = NumericalParameters(
        physical_dt=float(dt),
        sigma=float(sigma),
        electrical_tolerance=float(simulation_electrical_tolerance),
        tdgl_max_normalized_timestep=tdgl_max_normalized_timestep,
        thermal_max_substep=float(simulation_thermal_max_substep),
    )

    all_history = []
    attempt_histories = []
    adaptation_events = []

    total_iterations = 0
    failed_steps = 0
    total_restarts = 0
    checkpoint_count = 0
    adaptation_recommendation_count = 0
    oscillation_count = 0
    oscillation_restart_count = 0
    successful_adaptations = 0
    prediction_errors = []
    final_result = None
    physical_step_summaries = []

    kymograph_temperature = []
    kymograph_psi = []
    kymograph_times = []

    mesh = simulation.mesh
    center_y = mesh.ny // 2

    tracker = 0
    if hasattr(mesh, "dx"):
        kymograph_x = (
            np.arange(mesh.nx, dtype=float)
            - (mesh.nx - 1) / 2.0
        ) * float(mesh.dx) * 1e6
        kymograph_x_label = "Position across film (µm)"
    else:
        kymograph_x = np.arange(mesh.nx, dtype=float)
        kymograph_x_label = "Position across film (cell index)"

    physical_time = 0.0
    start_time = time.perf_counter()

    # These track the failed attempt immediately preceding the current retry.
    # If its adaptation did not materially improve convergence, the next
    # escalation is a reduction of the coupled physical timestep.
    previous_failed_result = None
    previous_adaptation_parameters = []

    adaptation_effectiveness_config = controller_config.get(
        "adaptation_effectiveness",
        {},
    )
    minimum_residual_improvement = float(
        adaptation_effectiveness_config.get(
            "minimum_residual_improvement",
            0.10,
        )
    )
    minimum_iteration_improvement = float(
        adaptation_effectiveness_config.get(
            "minimum_iteration_improvement",
            0.10,
        )
    )

    timing_history = {
        "tdgl": [],
        "electrical": [],
        "thermal": []
        }

    # Each outer physical step represents one requested interval of the
    # original dt. If convergence forces dt downward, that interval is
    # completed using additional accepted substeps. Thus reducing dt does
    # not silently shorten the requested physical simulation time.
    target_step_duration = float(dt)

    interval = max(1, physical_steps // 100)

    for physical_step in range(physical_steps):
        time_start = time.time()


        accepted_state = copy.deepcopy(
            copy_trial_state(simulation)
        )
        outer_step_start_events = len(adaptation_events)

        remaining_step_time = target_step_duration
        timestep_attempt = 0
        timestep_converged = False
        substep_index = 0
        substep_start_time = physical_time


        while remaining_step_time > max(
            1e-30,
            target_step_duration * 1e-15,
        ):
            # The current adaptive dt may not overshoot the requested outer
            # physical-step endpoint.
            parameters.physical_dt = min(
                parameters.physical_dt,
                remaining_step_time,
            )

            timestep_attempt += 1

            if timestep_attempt > int(
                controller_config["retry"]["max_restarts_per_timestep"]
            ):
                failed_steps += 1
                physical_step_summaries.append({
                    "physical_step": physical_step + 1,
                    "attempts": timestep_attempt - 1,
                    "substeps_completed": substep_index,
                    "converged": False,
                    "final_dt": parameters.physical_dt,
                    "physical_time": physical_time,
                    "reason": "maximum_restarts_exceeded",
                })
                timestep_converged = False
                break

            apply_trial_state(
                simulation,
                accepted_state,
            )

            attempt_start_events = len(adaptation_events)
            result = solve_coupled_attempt(
                simulation=simulation,
                accepted_state=accepted_state,
                parameters=parameters,
                ordering=ordering,
                thermal_model=thermal_model,
                voltage_left=voltage_left,
                voltage_right=voltage_right,
                external_heat=external_heat,
                tolerance=tolerance,
                max_iterations=max_iterations,
                check_interval=check_interval,
                prediction_window=prediction_window,
                reasonable_iterations=reasonable_iterations,
                adaptive_relaxation=adaptive_relaxation,
                minimum_sigma=minimum_sigma,
                reduction_factor=reduction_factor,
                required_reversals=required_reversals,
                adaptive_controller_config=controller_config,
            )
            final_result = result
            total_iterations += result.iterations
            checkpoint_count += sum(result.checkpoint_history)
            adaptation_recommendation_count += sum(result.adaptation_history)
            oscillation_count += sum(result.oscillation_history)
            prediction_errors.extend(result.prediction_errors)
            for component, elapsed in result.timing_history.items():
                timing_history[component].extend(elapsed)
            
            if result.converged:
                attempt_history = build_attempt_history(
                    result=result,
                    physical_step=physical_step + 1,
                    attempt=timestep_attempt,
                    parameters=parameters,
                    restart_reason="converged",
                )

                attempt_histories.append({
                    "physical_step": physical_step + 1,
                    "attempt": timestep_attempt,
                    "substep": substep_index + 1,
                    "converged": True,
                    "restart_reason": "converged",
                    "history": attempt_history,
                })
                all_history.extend(attempt_history)

                timestep_converged = True
                accepted_state = copy.deepcopy(
                    copy_trial_state(simulation)
                )

                # Only an accepted physical state advances physical time and
                # enters the kymographs.
                accepted_dt = parameters.physical_dt
                physical_time += accepted_dt
                remaining_step_time = max(
                    0.0,
                    remaining_step_time - accepted_dt,
                )
                substep_index += 1

                kymograph_temperature.append(
                    simulation.fields.temperature[center_y, :].copy()
                )
                kymograph_psi.append(
                    np.abs(
                        simulation.fields.psi[center_y, :]
                    ).copy()
                )
                kymograph_times.append(physical_time)

                # A successful substep completes this retry sequence. If the
                # outer requested interval still has time remaining, continue
                # with another accepted substep from this newly accepted state.
                previous_failed_result = None
                previous_adaptation_parameters = []

                if remaining_step_time <= max(
                    1e-30,
                    target_step_duration * 1e-15,
                ):
                    physical_step_summaries.append({
                        "physical_step": physical_step + 1,
                        "attempts": timestep_attempt,
                        "restarts": timestep_attempt - 1,
                        "substeps_completed": substep_index,
                        "converged": True,
                        "final_dt": parameters.physical_dt,
                        "final_sigma": parameters.sigma,
                        "physical_time": physical_time,
                        "step_duration": physical_time - substep_start_time,
                    })

                    for event in adaptation_events[outer_step_start_events:]:
                        if event.successful_after_restart is None:
                            event.successful_after_restart = True
                            successful_adaptations += 1
                    if physical_step % interval == 0:
                        percentage = (physical_step / physical_steps)*100
                        print("Completion percentage:", percentage, "%")
                        time_end =time.time()
                        time_per_percent = round(time_end-time_start, 3)
                        print(time_per_percent, "seconds")
        

                # Reset the attempt counter for a new substep. This prevents
                # successful subdivision from consuming the restart budget.
                timestep_attempt = 0
                continue

            # Failed coupling attempt.
            timestep_converged = False

            if not bool(controller_config["enabled"]):
                failed_steps += 1
                attempt_history = build_attempt_history(
                    result=result,
                    physical_step=physical_step + 1,
                    attempt=timestep_attempt,
                    parameters=parameters,
                    restart_reason="failed_adaptive_controller_disabled",
                )
                attempt_histories.append({
                    "physical_step": physical_step + 1,
                    "attempt": timestep_attempt,
                    "substep": substep_index + 1,
                    "converged": False,
                    "restart_reason": "failed_adaptive_controller_disabled",
                    "history": attempt_history,
                })
                all_history.extend(attempt_history)
                physical_step_summaries.append({
                    "physical_step": physical_step + 1,
                    "attempts": timestep_attempt,
                    "substeps_completed": substep_index,
                    "converged": False,
                    "final_dt": parameters.physical_dt,
                    "physical_time": physical_time,
                    "reason": "adaptive_controller_disabled",
                })
                break

            final_status = (
                result.status_history[-1]
                if result.status_history
                else ConvergenceStatus.START.value
            )

            oscillating = any(result.oscillation_history)
            actions = []
            problematic = []

            # If the previous adaptation was component-specific and did not
            # help enough, escalate immediately to coupled physical dt.
            adaptation_was_ineffective = False
            if (
                previous_failed_result is not None
                and previous_adaptation_parameters
            ):
                adaptation_was_ineffective = not adaptation_was_effective(
                    previous_failed_result,
                    result,
                    minimum_residual_improvement=minimum_residual_improvement,
                    minimum_iteration_improvement=minimum_iteration_improvement,
                )

            if adaptation_was_ineffective:
                actions = choose_physical_dt_fallback(
                    parameters,
                    controller_config,
                    bounds,
                )
            else:
                if (
                    oscillating
                    and controller_config["behavior"][
                        "reduce_sigma_only_for_oscillation"
                    ]
                ):
                    old_sigma = parameters.sigma
                    new_sigma = max(
                        minimum_sigma,
                        old_sigma * reduction_factor,
                    )

                    if new_sigma < old_sigma:
                        actions.append(
                            (
                                "sigma",
                                old_sigma,
                                new_sigma,
                            )
                        )
                        oscillation_restart_count += 1
                else:
                    residual_histories = {
                        "temperature": result.temperature_residual,
                        "psi": result.psi_residual,
                        "voltage": result.voltage_residual,
                    }
                    problematic = identify_problematic_components(
                        residual_histories=residual_histories,
                        tolerance=tolerance,
                        controller_config=controller_config,
                    )
                    actions = choose_adaptation(
                        problematic_components=problematic,
                        parameters=parameters,
                        controller_config=controller_config,
                        bounds=bounds,
                    )

                # If the normal hierarchy found nothing it can change, use
                # the same physical-dt fallback before declaring failure.
                if not actions:
                    actions = choose_physical_dt_fallback(
                        parameters,
                        controller_config,
                        bounds,
                    )

            if not actions:
                failed_steps += 1
                attempt_history = build_attempt_history(
                    result=result,
                    physical_step=physical_step + 1,
                    attempt=timestep_attempt,
                    parameters=parameters,
                    restart_reason="failed_no_further_adaptation_possible",
                )
                attempt_histories.append({
                    "physical_step": physical_step + 1,
                    "attempt": timestep_attempt,
                    "substep": substep_index + 1,
                    "converged": False,
                    "restart_reason": "failed_no_further_adaptation_possible",
                    "history": attempt_history,
                })
                all_history.extend(attempt_history)
                physical_step_summaries.append({
                    "physical_step": physical_step + 1,
                    "attempts": timestep_attempt,
                    "substeps_completed": substep_index,
                    "converged": False,
                    "final_dt": parameters.physical_dt,
                    "physical_time": physical_time,
                    "reason": "no_further_adaptation_possible",
                })
                break

            if any(parameter == "sigma" for parameter, _, _ in actions):
                restart_reason = "oscillation"
            elif any(parameter == "physical_dt" for parameter, _, _ in actions):
                restart_reason = (
                    "physical timestep reduction"
                    if not adaptation_was_ineffective
                    else "ineffective adaptation so physical timestep reduction"
                )

            else:
                parameters_changed = [
                    parameter
                    for parameter, _, _ in actions
                ]
                restart_reason = (
                    "numerical adaptation: "
                    + ", ".join(parameters_changed)
                )
            print(restart_reason)
            attempt_history = build_attempt_history(
                result=result,
                physical_step=physical_step + 1,
                attempt=timestep_attempt,
                parameters=parameters,
                restart_reason=restart_reason,
            )
            attempt_histories.append({
                "physical_step": physical_step + 1,
                "attempt": timestep_attempt,
                "substep": substep_index + 1,
                "converged": False,
                "restart_reason": restart_reason,
                "history": attempt_history,
            })
            all_history.extend(attempt_history)

            for parameter, old_value, new_value in actions:
                setattr(
                    parameters,
                    parameter,
                    new_value,
                )

                event = AdaptationEvent(
                    physical_step=physical_step + 1,
                    attempt=timestep_attempt,
                    iteration=result.iterations,
                    reason=(
                        "oscillation"
                        if parameter == "sigma"
                        else (
                            "component_adaptation"
                            if parameter != "physical_dt"
                            else (
                                "physical_timestep_fallback"
                                if adaptation_was_ineffective
                                else "physical_timestep_retry"
                            )
                        )
                    ),
                    parameter=parameter,
                    old_value=old_value,
                    new_value=new_value,
                    residual=result.final_residual,
                    status=final_status,
                    predicted_iterations=(
                        result.estimated_remaining_history[-1]
                        if result.estimated_remaining_history
                        else None
                    ),
                    prediction_confidence=(
                        result.confidence_history[-1]
                        if result.confidence_history
                        else 0.0
                    ),
                    problematic_components=(
                        problematic
                        if not oscillating
                        else []
                    ),
                )
                adaptation_events.append(event)

            total_restarts += 1
            print(total_restarts)

            if total_restarts > int(
                controller_config["retry"]["max_total_restarts"]
            ):
                failed_steps += 1
                physical_step_summaries.append({
                    "physical_step": physical_step + 1,
                    "attempts": timestep_attempt,
                    "substeps_completed": substep_index,
                    "converged": False,
                    "final_dt": parameters.physical_dt,
                    "physical_time": physical_time,
                    "reason": "maximum_total_restarts_exceeded",
                })
                break

            # The next attempt always begins from the last accepted state.
            apply_trial_state(
                simulation,
                accepted_state,
            )

            previous_failed_result = result
            previous_adaptation_parameters = [
                parameter
                for parameter, _, _ in actions
            ]

        # Stop processing later outer steps after an unrecoverable failure.
        if not timestep_converged and remaining_step_time > max(
            1e-30,
            target_step_duration * 1e-15,
        ):
            continue

    runtime = time.perf_counter() - start_time
    fields = simulation.fields

    valid_prediction_errors = [
        value
        for value in prediction_errors
        if np.isfinite(value)
    ]

    status_counts = {}
    for row in all_history:
        status = row["status"]
        status_counts[status] = status_counts.get(status, 0) + 1


    adaptation_reason_counts = {}
    adaptation_parameter_counts = {}
    for event in adaptation_events:
        adaptation_reason_counts[event.reason] = (
            adaptation_reason_counts.get(event.reason, 0) + 1
        )
        adaptation_parameter_counts[event.parameter] = (
            adaptation_parameter_counts.get(event.parameter, 0) + 1
        )

    unsuccessful_adaptations = sum(
        1
        for event in adaptation_events
        if event.successful_after_restart is False
    )

    mean_prediction_error = (
        float(np.mean(valid_prediction_errors))
        if valid_prediction_errors
        else None
    )

    mean_prediction_confidence = (
        float(np.mean([
            row["prediction_confidence"]
            for row in all_history
        ]))
        if all_history
        else None
    )


    return {
        "dt": dt,
        "initial_sigma": sigma,
        "final_sigma": parameters.sigma,
        "final_physical_dt": parameters.physical_dt,
        "final_electrical_tolerance": parameters.electrical_tolerance,
        "final_tdgl_max_normalized_timestep": parameters.tdgl_max_normalized_timestep,
        "final_thermal_max_substep": parameters.thermal_max_substep,
        "ordering": ordering,
        "converged": failed_steps == 0,
        "physical_steps": physical_steps,
        "failed_steps": failed_steps,
        "total_iterations": total_iterations,
        "average_iterations": total_iterations / physical_steps,
        "final_residual": (
            final_result.final_residual
            if final_result is not None
            else float("inf")
        ),
        "runtime_seconds": runtime,
        "runtime_per_physical_step": runtime / physical_steps,
        "final_physical_time": physical_time,
        "final_temperature_mean": float(np.mean(fields.temperature)),
        "final_temperature_max": float(np.max(fields.temperature)),
        "final_psi_amplitude_mean": float(np.mean(np.abs(fields.psi))),
        "final_psi_amplitude_min": float(np.min(np.abs(fields.psi))),
        "final_current_mean": float(np.mean(np.sqrt(
            fields.current_density_x ** 2
            + fields.current_density_y ** 2
        ))),
        "final_heat_mean": float(np.mean(fields.heat_source)),
        "checkpoint_count": checkpoint_count,
        "adaptation_recommendation_count": adaptation_recommendation_count,
        "adaptation_event_count": len(adaptation_events),
        "restart_count": total_restarts,
        "oscillation_count": oscillation_count,
        "oscillation_restart_count": oscillation_restart_count,
        "successful_adaptation_count": successful_adaptations,
        "unsuccessful_adaptation_count": unsuccessful_adaptations,
        "mean_prediction_error": mean_prediction_error,
        "mean_prediction_confidence": mean_prediction_confidence,
        "status_counts": status_counts,
        "adaptation_reason_counts": adaptation_reason_counts,
        "adaptation_parameter_counts": adaptation_parameter_counts,
        "physical_step_summaries": physical_step_summaries,
        "attempt_summaries": [
            build_attempt_summary(attempt)
            for attempt in attempt_histories
        ],
        "adaptation_events": [
            {
                "physical_step": event.physical_step,
                "attempt": event.attempt,
                "iteration": event.iteration,
                "reason": event.reason,
                "parameter": event.parameter,
                "old_value": event.old_value,
                "new_value": event.new_value,
                "residual": event.residual,
                "status": event.status,
                "predicted_iterations": event.predicted_iterations,
                "prediction_confidence": event.prediction_confidence,
                "problematic_components": event.problematic_components,
                "successful_after_restart": event.successful_after_restart,
            }
            for event in adaptation_events
        ],
        # Detailed histories remain available to CSV/diagnostic plotting but
        # are explicitly excluded from summary.json by save_results().
        "history": all_history,
        "attempt_histories": attempt_histories,
        "_kymograph_temperature": kymograph_temperature,
        "_kymograph_psi": kymograph_psi,
        "_kymograph_times": kymograph_times,
        "_kymograph_x": kymograph_x,
        "_kymograph_x_label": kymograph_x_label,
        "timing_history": timing_history
    }


# ============================================================================
# Plotting
# ============================================================================


def plot_attempt_residual(
    attempt,
    output_path,
    tolerance,
):
    """
    Plot the convergence behavior of exactly one coupling attempt.

    Each physical timestep attempt gets its own graph. This prevents
    restarted attempts from being visually stitched together and makes
    the numerical adaptation decision visible at the point where it occurs.
    """

    history = attempt["history"]

    if not history:
        return

    iterations = [
        row["iteration"]
        for row in history
    ]

    total = [
        row["residual"]
        for row in history
    ]

    temperature = [
        row["temperature_residual"]
        for row in history
    ]

    psi = [
        row["psi_residual"]
        for row in history
    ]

    voltage = [
        row["voltage_residual"]
        for row in history
    ]

    physical_step = attempt[
        "physical_step"
    ]

    attempt_number = attempt[
        "attempt"
    ]

    reason = attempt[
        "restart_reason"
    ]

    plt.figure(
        figsize=(9, 6)
    )

    plt.semilogy(
        iterations,
        total,
        marker="o",
        label="Total",
    )

    plt.semilogy(
        iterations,
        temperature,
        marker="o",
        label="Thermal",
    )

    plt.semilogy(
        iterations,
        psi,
        marker="o",
        label="TDGL",
    )

    plt.semilogy(
        iterations,
        voltage,
        marker="o",
        label="Electrical",
    )

    plt.axhline(
        tolerance,
        linestyle="--",
        label="Tolerance",
    )

    plt.xlabel(
        "Coupling iteration"
    )

    plt.ylabel(
        "Relative residual"
    )

    if reason == "converged":

        title = (
            f"Step {physical_step} — "
            f"Attempt {attempt_number} — "
            f"Converged"
        )

    else:

        title = (
            f"Step {physical_step} — "
            f"Attempt {attempt_number} — "
            f"Restart: {reason}"
        )

    plt.title(
        title
    )

    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    plt.savefig(
        output_path,
        dpi=150,
    )

    plt.close()


# ============================================================================
# Legacy / secondary diagnostic plots
# ============================================================================


def plot_field_residuals(
    result,
    output_path,
):
    history = result["history"]

    if not history:
        return

    iterations = [
        row["iteration"]
        for row in history
    ]

    fields = (
        "temperature",
        "psi",
        "voltage",
        "ax",
        "ay",
    )

    plt.figure()

    for field in fields:

        values = [
            row[
                f"{field}_residual"
            ]
            for row in history
        ]

        plt.semilogy(
            iterations,
            values,
            marker="o",
            label=field,
        )

    plt.xlabel(
        "Coupling iteration"
    )

    plt.ylabel(
        "Relative residual"
    )

    plt.title(
        "Individual field residuals"
    )

    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    plt.savefig(
        output_path,
        dpi=150,
    )

    plt.close()


def plot_sigma_history(
    result,
    output_path,
):
    history = result["history"]

    if not history:
        return

    iterations = [
        row["iteration"]
        for row in history
    ]

    sigma = [
        row["sigma"]
        for row in history
    ]

    plt.figure()

    plt.plot(
        iterations,
        sigma,
        marker="o",
    )

    plt.xlabel(
        "Coupling iteration"
    )

    plt.ylabel(
        "Relaxation sigma"
    )

    plt.title(
        "Relaxation history"
    )

    plt.grid(True)
    plt.tight_layout()

    plt.savefig(
        output_path,
        dpi=150,
    )

    plt.close()


def plot_convergence_prediction(
    result,
    output_path,
):
    history = result["history"]

    valid_rows = [
        row
        for row in history
        if (
            row[
                "estimated_iterations_remaining"
            ]
            is not None
            and row[
                "actual_iterations_remaining"
            ]
            is not None
        )
    ]

    if not valid_rows:
        return

    iterations = [
        row["iteration"]
        for row in valid_rows
    ]

    predicted = [
        row[
            "estimated_iterations_remaining"
        ]
        for row in valid_rows
    ]

    actual = [
        row[
            "actual_iterations_remaining"
        ]
        for row in valid_rows
    ]

    plt.figure()

    plt.plot(
        iterations,
        predicted,
        marker="o",
        label="predicted",
    )

    plt.plot(
        iterations,
        actual,
        marker="x",
        label="actual",
    )

    plt.xlabel(
        "Coupling iteration"
    )

    plt.ylabel(
        "Iterations remaining"
    )

    plt.title(
        "Convergence prediction"
    )

    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    plt.savefig(
        output_path,
        dpi=150,
    )

    plt.close()


def plot_prediction_confidence(
    result,
    output_path,
):
    history = result["history"]

    if not history:
        return

    iterations = [
        row["iteration"]
        for row in history
    ]

    confidence = [
        row[
            "prediction_confidence"
        ]
        for row in history
    ]

    plt.figure()

    plt.plot(
        iterations,
        confidence,
        marker="o",
    )

    plt.xlabel(
        "Coupling iteration"
    )

    plt.ylabel(
        "Prediction confidence"
    )

    plt.ylim(
        0.0,
        1.05,
    )

    plt.title(
        "Convergence prediction confidence"
    )

    plt.grid(True)
    plt.tight_layout()

    plt.savefig(
        output_path,
        dpi=150,
    )

    plt.close()


def plot_kymograph(
    data,
    times,
    x_positions,
    x_label,
    output_path,
    title,
    colorbar_label,
):
    """Plot a centerline field against physical simulation time."""
    if not data or not times:
        return

    values = np.asarray(data, dtype=float)
    times = np.asarray(times, dtype=float)
    x_positions = np.asarray(x_positions, dtype=float)

    if values.ndim != 2 or values.shape[0] != len(times):
        return

    if values.shape[1] != len(x_positions):
        return

    plt.figure(figsize=(9, 6))

    # pcolormesh uses the actual time coordinates, so adaptive/subdivided
    # timesteps are represented at their true physical times rather than being
    # incorrectly treated as uniformly spaced rows.
    mesh = plt.pcolormesh(
        x_positions,
        times,
        values,
        shading="auto",
    )

    plt.xlabel(x_label)
    plt.ylabel("Physical simulation time (s)")
    plt.title(title)
    plt.colorbar(mesh, label=colorbar_label)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()


def plot_coupling_timings(timing_history, output_path):
    tdgl = np.asarray(timing_history["tdgl"])
    electrical = np.asarray(timing_history["electrical"])
    thermal = np.asarray(timing_history["thermal"])

    total = tdgl + electrical + thermal

    # Avoid division by zero
    percentages = {
        "tdgl": np.divide(
            tdgl,
            total,
            out=np.zeros_like(tdgl),
            where=total != 0,
        ) * 100,

        "electrical": np.divide(
            electrical,
            total,
            out=np.zeros_like(electrical),
            where=total != 0,
        ) * 100,

        "thermal": np.divide(
            thermal,
            total,
            out=np.zeros_like(thermal),
            where=total != 0,
        ) * 100,
    }

    iterations = range(1, len(tdgl) + 1)

    plt.figure(figsize=(10, 6))

    plt.plot(
        iterations,
        percentages["tdgl"],
        label="TDGL",
    )

    plt.plot(
        iterations,
        percentages["electrical"],
        label="Electrical",
    )

    plt.plot(
        iterations,
        percentages["thermal"],
        label="Thermal",
    )

    plt.xlabel("Coupling iteration")
    plt.ylabel("Runtime fraction (%)")
    plt.title("Coupled solver subsystem runtime distribution")
    plt.ylim(0, 100)
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close
# ============================================================================
# Configuration / output
# ============================================================================


def load_config(
    path,
):
    with open(
        path,
        "r",
        encoding="utf-8",
    ) as f:

        return json.load(f)


def select_orderings(
    config,
):
    mode = config.get(
        "order_mode",
        "one",
    )

    if mode == "all":
        return list(
            DEFAULT_ORDERINGS
        )

    if mode == "one":

        ordering = config[
            "single_ordering"
        ]

        if ordering not in DEFAULT_ORDERINGS:

            raise ValueError(
                f"Unknown ordering: "
                f"{ordering}"
            )

        return [ordering]

    raise ValueError(
        "order_mode must be "
        "'one' or 'all'"
    )


def save_results(
    results,
    output_directory,
    save_history,
    make_plots,
    tolerance,
):
    output_directory.mkdir(
        parents=True,
        exist_ok=False,
    )

    excluded_from_json = {
        "history",
        "attempt_histories",
        "_kymograph_temperature",
        "_kymograph_psi",
        "_kymograph_times",
        "_kymograph_x",
        "_kymograph_x_label",
    }

    compact_results = []
    for result in results:
        compact = {
            key: value
            for key, value in result.items()
            if key not in excluded_from_json
        }
        compact_results.append(compact)

    with open(
        output_directory / "summary.json",
        "w",
        encoding="utf-8",
    ) as f:
        # IMPORTANT: write the compact collection, not the full last result.
        json.dump(
            compact_results,
            f,
            indent=2,
        )

    if compact_results:
        fieldnames = list(compact_results[0].keys())
        with open(
            output_directory / "summary.csv",
            "w",
            newline="",
            encoding="utf-8",
        ) as f:
            writer = csv.DictWriter(
                f,
                fieldnames=fieldnames,
                extrasaction="ignore",
            )
            writer.writeheader()
            writer.writerows(compact_results)

    if save_history:
        history_fieldnames = [
            "dt",
            "initial_sigma",
            "ordering",
            "physical_step",
            "attempt",
            "restart_reason",
            "iteration",
            "residual",
            "temperature_residual",
            "psi_residual",
            "voltage_residual",
            "ax_residual",
            "ay_residual",
            "sigma",
            "physical_dt",
            "electrical_tolerance",
            "tdgl_max_normalized_timestep",
            "thermal_max_substep",
            "oscillation",
            "residual_ratio",
            "log_slope",
            "estimated_iterations_remaining",
            "actual_iterations_remaining",
            "prediction_error",
            "prediction_confidence",
            "status",
            "checkpoint",
            "adaptation_recommended",
        ]

        with open(
            output_directory / "residual_history.csv",
            "w",
            newline="",
            encoding="utf-8",
        ) as f:
            writer = csv.DictWriter(
                f,
                fieldnames=history_fieldnames,
                extrasaction="ignore",
            )
            writer.writeheader()

            for result in results:
                history = result["history"]

                save_interval = max(
                    1,
                    math.ceil(len(history) / 20),
                )

                for index, row in enumerate(history):
                    if index % save_interval != 0:
                        continue
                    writer.writerow({
                        "dt": result["dt"],
                        "initial_sigma": result["initial_sigma"],
                        "ordering": result["ordering"],
                        **row,
                    })

        event_fieldnames = [
            "dt",
            "initial_sigma",
            "ordering",
            "physical_step",
            "attempt",
            "iteration",
            "reason",
            "parameter",
            "old_value",
            "new_value",
            "residual",
            "status",
            "predicted_iterations",
            "prediction_confidence",
            "problematic_components",
            "successful_after_restart",
        ]

        with open(
            output_directory / "adaptation_events.csv",
            "w",
            newline="",
            encoding="utf-8",
        ) as f:
            writer = csv.DictWriter(
                f,
                fieldnames=event_fieldnames,
            )
            writer.writeheader()

            for result in results:
                for event in result["adaptation_events"]:
                    writer.writerow({
                        "dt": result["dt"],
                        "initial_sigma": result["initial_sigma"],
                        "ordering": result["ordering"],
                        **event,
                    })

    if make_plots:
        for index, result in enumerate(results, start=1):
            prefix = f"run_{index:03d}"

            attempts = result["attempt_histories"]

            save_interval = max(
                1,
                math.ceil(len(attempts) / 20),
            )

            for attempt_index, attempt in enumerate(attempts):
                if attempt_index % save_interval != 0:
                    continue

                step = attempt["physical_step"]
                attempt_number = attempt["attempt"]
                reason = attempt["restart_reason"]

                safe_reason = (
                    reason.lower()
                    .replace(" ", "_")
                    .replace(":", "")
                    .replace(",", "")
                )

                filename = (
                    f"{prefix}"
                    f"_step_{step:03d}"
                    f"_attempt_{attempt_number:02d}"
                    f"_{safe_reason}.png"
                )

                plot_attempt_residual(
                    attempt=attempt,
                    output_path=output_directory / filename,
                    tolerance=tolerance,
                )

            plot_kymograph(
                data=result["_kymograph_temperature"],
                times=result["_kymograph_times"],
                x_positions=result["_kymograph_x"],
                x_label=result["_kymograph_x_label"],
                output_path=(
                    output_directory
                    / f"{prefix}_temperature_kymograph.png"
                ),
                title=(
                    f"{result['ordering']} — "
                    "Temperature centerline kymograph"
                ),
                colorbar_label="Temperature",
            )

            plot_kymograph(
                data=result["_kymograph_psi"],
                times=result["_kymograph_times"],
                x_positions=result["_kymograph_x"],
                x_label=result["_kymograph_x_label"],
                output_path=(
                    output_directory
                    / f"{prefix}_psi_kymograph.png"
                ),
                title=(
                    f"{result['ordering']} — "
                    "|ψ| centerline kymograph"
                ),
                colorbar_label="|ψ|",
            )
            plot_coupling_timings(result['timing_history'], output_directory / "coupling_timing.png")


# ============================================================================
# Main
# ============================================================================


def main():
    parser = argparse.ArgumentParser(
        description=(
            "SHS adaptive coupled convergence benchmark"
        )
    )

    parser.add_argument(
        "--config",
        default=str(
            CONFIG_PATH
        ),
    )

    parser.add_argument(
        "--iterations",
        type=int,
        default=None,
    )

    parser.add_argument(
        "--steps",
        type=int,
        default=None,
    )

    parser.add_argument(
        "--orders",
        choices=(
            "one",
            "all",
        ),
        default=None,
    )

    parser.add_argument(
        "--ordering",
        default=None,
    )

    parser.add_argument(
        "--dt",
        type=float,
        nargs="+",
        default=None,
    )

    parser.add_argument(
        "--sigma",
        type=float,
        nargs="+",
        default=None,
    )

    args = parser.parse_args()

    config = load_config(
        args.config
    )

    if args.iterations is not None:
        config[
            "convergence"
        ][
            "max_iterations"
        ] = args.iterations

    if args.steps is not None:
        config[
            "physical_steps"
        ] = args.steps

    if args.dt is not None:
        config[
            "dt_values"
        ] = args.dt

    if args.sigma is not None:
        config[
            "sigma_values"
        ] = args.sigma

    if args.orders is not None:
        config[
            "order_mode"
        ] = args.orders

    if args.ordering is not None:
        config[
            "single_ordering"
        ] = args.ordering

    orderings = select_orderings(
        config
    )

    results = []

    for dt in config[
        "dt_values"
    ]:

        for sigma in config[
            "sigma_values"
        ]:

            for ordering in orderings:

                print()
                print(
                    "=" * 70
                )

                print(
                    f"dt = {dt:.3e}"
                )

                print(
                    f"sigma = {sigma:.3f}"
                )

                print(
                    f"ordering = {ordering}"
                )

                print(
                    "max iterations = "
                    f"{config['convergence']['max_iterations']}"
                )

                print(
                    "tolerance = "
                    f"{config['convergence']['tolerance']:.3e}"
                )

                print(
                    "adaptive controller = "
                    f"{config['adaptive_controller']['enabled']}"
                )

                print(
                    "=" * 70
                )

                result = run_single(
                    config=config,
                    dt=dt,
                    sigma=sigma,
                    ordering=ordering,
                )


                results.append(
                    result
                )

                
                status = (
                    "CONVERGED"
                    if result[
                        "converged"
                    ]
                    else "FAILED"
                )

                print(
                    f"status = {status}"
                )

                print(
                    "iterations = "
                    f"{result['total_iterations']}"
                )

                print(
                    "final residual = "
                    f"{result['final_residual']:.6e}"
                )

                print(
                    "final sigma = "
                    f"{result['final_sigma']:.3f}"
                )

                print(
                    "final physical dt = "
                    f"{result['final_physical_dt']:.3e}"
                )

                print(
                    "restarts = "
                    f"{result['restart_count']}"
                )

                print(
                    "checkpoints = "
                    f"{result['checkpoint_count']}"
                )

                print(
                    "adaptation events = "
                    f"{result['adaptation_event_count']}"
                )

                print(
                    "oscillation detections = "
                    f"{result['oscillation_count']}"
                )

                print(
                    "oscillation restarts = "
                    f"{result['oscillation_restart_count']}"
                )

                print(
                    "successful adaptations = "
                    f"{result['successful_adaptation_count']}"
                )

                print(
                    "unsuccessful adaptations = "
                    f"{result['unsuccessful_adaptation_count']}"
                )

                if (
                    result[
                        "mean_prediction_error"
                    ]
                    is not None
                ):

                    print(
                        "mean prediction error = "
                        f"{result['mean_prediction_error']:.3%}"
                    )

                if (
                    result[
                        "mean_prediction_confidence"
                    ]
                    is not None
                ):

                    print(
                        "mean prediction confidence = "
                        f"{result['mean_prediction_confidence']:.3f}"
                    )

                print(
                    "adaptation parameters = "
                    f"{result['adaptation_parameter_counts']}"
                )

                print(
                    "runtime = "
                    f"{result['runtime_seconds']:.3f} s"
                )


    base_directory = Path(
        config[
            "output"
        ][
            "directory"
        ]
    )
    output_directory = base_directory
    counter = 1
    while output_directory.exists():
        output_directory = base_directory.parent / f"{base_directory.name}_{counter}"
        counter+=1


    save_results(
        results=results,
        output_directory=output_directory,
        save_history=config[
            "output"
        ][
            "save_history"
        ],
        make_plots=config[
            "output"
        ][
            "make_plots"
        ],
        tolerance=config[
            "convergence"
        ][
            "tolerance"
        ],
    )

    print()
    print(
        "Benchmark complete."
    )

    print(
        "Results written to: "
        f"{output_directory}"
    )


if __name__ == "__main__":
    main()