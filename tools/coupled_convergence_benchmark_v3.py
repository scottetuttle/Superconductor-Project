"""
Coupled convergence benchmark for SHS.

This benchmark measures convergence of the coupled TDGL/electrical/thermal
system at a single physical timestep.

Important distinction:

    Physical time:
        t_n -> t_(n+1)

    Coupling iterations:
        k = 0, 1, 2, ...

The coupling iterations do NOT advance physical time repeatedly.

At every coupling iteration, TDGL starts from the accepted state at t_n,
while the trial temperature/vector potential/etc. are used to evaluate the
prospective state at t_(n+1).

Relaxation is applied between the PREVIOUS TRIAL state and the newly
CALCULATED state:

    trial_(k+1) =
        trial_k + sigma * (calculated_k - trial_k)

This is important. Relaxing against the original accepted state every
iteration would prevent the Picard iteration from behaving correctly.
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from shs.config.builder import build_simulation
from shs.solvers.tdgl_solver import tdgl_step
from shs.solvers.electrical_solver import electrical_step
from shs.solvers.thermal_solver import thermal_step
from shs.physics.thermal import ThermalModel

# Adjust this import if your project uses a different TDGL module path.
from shs.tdgl import TDGLModel
from shs.tdgl.parameters import TDGLParameters


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

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


@dataclass
class TrialState:
    temperature: np.ndarray
    psi: np.ndarray
    voltage: np.ndarray
    vector_potential_x: np.ndarray
    vector_potential_y: np.ndarray


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
    final_sigma: float


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def relative_residual(old, new, floor=1e-14):
    """
    Relative L2 difference between two fields.

        ||new-old|| / max(||old||, floor)
    """
    old = np.asarray(old)
    new = np.asarray(new)

    numerator = np.linalg.norm((new - old).ravel())
    denominator = max(np.linalg.norm(old.ravel()), floor)

    return float(numerator / denominator)


def copy_trial_state(simulation) -> TrialState:
    fields = simulation.fields

    return TrialState(
        temperature=fields.temperature.copy(),
        psi=fields.psi.copy(),
        voltage=fields.voltage.copy(),
        vector_potential_x=fields.vector_potential_x.copy(),
        vector_potential_y=fields.vector_potential_y.copy(),
    )


def apply_trial_state(simulation, trial: TrialState):
    fields = simulation.fields

    fields.temperature = trial.temperature.copy()
    fields.psi = trial.psi.copy()
    fields.voltage = trial.voltage.copy()
    fields.vector_potential_x = trial.vector_potential_x.copy()
    fields.vector_potential_y = trial.vector_potential_y.copy()


def calculate_residuals(old: TrialState, new: TrialState):
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
    old_trial: TrialState,
    calculated: TrialState,
    sigma: float,
) -> TrialState:
    """
    Under-relax the NEW trial state from the PREVIOUS trial state.

        x_(k+1) = x_k + sigma * (G(x_k) - x_k)
    """

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


def state_update_directions(old: TrialState, new: TrialState):
    """
    Return the direction vectors of the state changes.

    Used only for heuristic oscillation detection.
    """

    return {
        "temperature": (
            new.temperature - old.temperature
        ).ravel(),

        "psi": (
            new.psi - old.psi
        ).view(np.float64).ravel(),

        "voltage": (
            new.voltage - old.voltage
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
    """
    Detect a simple direction reversal.

    This is deliberately heuristic. It is not treated as proof of
    mathematical instability.
    """

    oscillating = False

    for name in previous_delta:
        old = previous_delta[name]
        new = current_delta[name]

        old_norm = np.linalg.norm(old)
        new_norm = np.linalg.norm(new)

        if old_norm < 1e-14 or new_norm < 1e-14:
            continue

        dot = np.vdot(old, new).real

        if dot < 0.0:
            oscillating = True
            break

    return oscillating


# ---------------------------------------------------------------------------
# Individual physics operations
# ---------------------------------------------------------------------------

def run_tdgl(
    simulation,
    dt,
    tdgl_model,
):
    tdgl_step(
        simulation,
        dt,
        tdgl_model,
    )


def run_electrical(
    simulation,
    electrical_voltage_left,
    electrical_voltage_right,
):
    fields = simulation.fields

    electrical_step(
        fields,
        simulation.mesh,
        simulation.material_map,
        simulation.contact_map,
        voltage_left=electrical_voltage_left,
        voltage_right=electrical_voltage_right,
        superconducting_fraction=np.abs(fields.psi) ** 2,
        superconducting_current_x=(
            fields.supercurrent_density_x
        ),
        superconducting_current_y=(
            fields.supercurrent_density_y
        ),
    )


def run_thermal(
    simulation,
    dt,
    thermal_model,
    external_heat=None,
):
    if external_heat is not None:
        original_heat = simulation.fields.heat_source.copy()

        simulation.fields.heat_source = (
            original_heat + external_heat
        )

        thermal_step(
            simulation,
            dt,
            thermal_model,
        )

        simulation.fields.heat_source = original_heat

    else:
        thermal_step(
            simulation,
            dt,
            thermal_model,
        )


# ---------------------------------------------------------------------------
# Coupling order
# ---------------------------------------------------------------------------

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

    for operation in operations:

        if operation == "tdgl":
            run_tdgl(
                simulation,
                dt,
                tdgl_model,
            )

        elif operation == "electrical":
            run_electrical(
                simulation,
                voltage_left,
                voltage_right,
            )

        elif operation == "thermal":
            run_thermal(
                simulation,
                dt,
                thermal_model,
                external_heat=external_heat,
            )

        else:
            raise ValueError(
                f"Unknown coupling operation: {operation}"
            )


# ---------------------------------------------------------------------------
# One physical timestep
# ---------------------------------------------------------------------------

def solve_coupled_timestep(
    simulation,
    accepted_state: TrialState,
    dt,
    ordering,
    tdgl_model,
    thermal_model,
    voltage_left,
    voltage_right,
    external_heat,
    initial_sigma,
    tolerance,
    max_iterations,
    adaptive_relaxation,
    minimum_sigma,
    reduction_factor,
    required_reversals,
):
    """
    Solve one physical timestep using a fixed-point coupling iteration.

    accepted_state is the physical state at t_n.

    trial_state is the current estimate of the state at t_(n+1).
    """

    # Initial guess for t_(n+1) is the accepted state at t_n.
    trial_state = copy.deepcopy(accepted_state)

    residual_history = []
    temperature_history = []
    psi_history = []
    voltage_history = []
    ax_history = []
    ay_history = []

    sigma_history = []
    oscillation_history = []

    sigma = initial_sigma

    previous_delta = None
    consecutive_reversals = 0

    converged = False

    for iteration in range(1, max_iterations + 1):

        # ---------------------------------------------------------------
        # IMPORTANT:
        #
        # TDGL must start from psi_n every iteration.
        # Temperature/vector-potential/etc. come from the current trial.
        # ---------------------------------------------------------------

        # Start from the accepted physical state.
        apply_trial_state(
            simulation,
            accepted_state,
        )

        # Replace the coupled trial variables with the current estimate.
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

        # ψ must start from ψ_n, not the previous trial ψ.
        simulation.fields.psi = (
            accepted_state.psi.copy()
        )

        # ---------------------------------------------------------------
        # Calculate G(trial)
        # ---------------------------------------------------------------

        execute_ordering(
            simulation=simulation,
            ordering=ordering,
            dt=dt,
            tdgl_model=tdgl_model,
            thermal_model=thermal_model,
            voltage_left=voltage_left,
            voltage_right=voltage_right,
            external_heat=external_heat,
        )

        calculated = copy_trial_state(simulation)

        # ---------------------------------------------------------------
        # Calculate convergence BEFORE relaxation.
        # ---------------------------------------------------------------

        residuals = calculate_residuals(
            trial_state,
            calculated,
        )

        residual_history.append(
            residuals["total"]
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

        sigma_history.append(sigma)

        # ---------------------------------------------------------------
        # Oscillation detection
        # ---------------------------------------------------------------

        current_delta = state_update_directions(
            trial_state,
            calculated,
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

        if (
            adaptive_relaxation
            and consecutive_reversals
            >= required_reversals
        ):
            new_sigma = max(
                minimum_sigma,
                sigma * reduction_factor,
            )

            if new_sigma < sigma:
                sigma = new_sigma

            consecutive_reversals = 0

        # ---------------------------------------------------------------
        # Convergence
        # ---------------------------------------------------------------
        print("Total", residuals["total"])
        print("psi", residuals["psi"])
        print("Temp", residuals["temperature"])
        print("Voltage", residuals["voltage"])

        if residuals["total"] < tolerance:
            converged = True

            trial_state = calculated

            break

        # ---------------------------------------------------------------
        # Correct Picard relaxation:
        #
        # new trial = old trial + sigma * (calculated-old trial)
        # ---------------------------------------------------------------

        trial_state = relax_trial_state(
            trial_state,
            calculated,
            sigma,
        )

        previous_delta = current_delta

    # Put final trial state into simulation.
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
        iterations=len(residual_history),
        final_residual=final_residual,
        residual_history=residual_history,
        temperature_residual=temperature_history,
        psi_residual=psi_history,
        voltage_residual=voltage_history,
        ax_residual=ax_history,
        ay_residual=ay_history,
        sigma_history=sigma_history,
        oscillation_history=oscillation_history,
        final_sigma=sigma,
    )


# ---------------------------------------------------------------------------
# External heat
# ---------------------------------------------------------------------------

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

    cx = (mesh.nx - 1) / 2.0
    cy = (mesh.ny - 1) / 2.0

    r2 = (
        (x - cx) ** 2
        + (y - cy) ** 2
    )

    source = amplitude * np.exp(
        -r2 / (2.0 * radius_cells ** 2)
    )

    return source


# ---------------------------------------------------------------------------
# One benchmark run
# ---------------------------------------------------------------------------

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
        bath_temperature=config["thermal"][
            "bath_temperature"
        ],
        thermal_relaxation_rate=config["thermal"][
            "thermal_relaxation_rate"
        ],
    )

    tdgl_model = TDGLModel(
        TDGLParameters()
    )

    voltage_left = config["electrical"][
        "voltage_left"
    ]

    voltage_right = config["electrical"][
        "voltage_right"
    ]

    external_heat = None

    if config["external_heat"]["enabled"]:
        external_heat = make_gaussian_heat_source(
            simulation,
            amplitude=config["external_heat"][
                "amplitude"
            ],
            radius_cells=config["external_heat"][
                "radius_cells"
            ],
        )

    physical_steps = int(
        config["physical_steps"]
    )

    max_iterations = int(
        config["max_iterations"]
    )

    tolerance = float(
        config["tolerance"]
    )

    adaptive_config = config[
        "adaptive_relaxation"
    ]

    adaptive_enabled = bool(
        adaptive_config["enabled"]
    )

    minimum_sigma = float(
        adaptive_config["minimum_sigma"]
    )

    reduction_factor = float(
        adaptive_config["reduction_factor"]
    )

    required_reversals = int(
        adaptive_config["required_reversals"]
    )

    all_residual_history = []

    total_iterations = 0
    failed_steps = 0

    final_result = None

    start_time = time.perf_counter()

    for physical_step in range(
        physical_steps
    ):
        accepted_state = copy_trial_state(
            simulation
        )

        result = solve_coupled_timestep(
            simulation=simulation,
            accepted_state=accepted_state,
            dt=dt,
            ordering=ordering,
            tdgl_model=tdgl_model,
            thermal_model=thermal_model,
            voltage_left=voltage_left,
            voltage_right=voltage_right,
            external_heat=external_heat,
            initial_sigma=sigma,
            tolerance=tolerance,
            max_iterations=max_iterations,
            adaptive_relaxation=adaptive_enabled,
            minimum_sigma=minimum_sigma,
            reduction_factor=reduction_factor,
            required_reversals=required_reversals,
        )

        total_iterations += result.iterations

        if not result.converged:
            failed_steps += 1

        final_result = result

        for iteration, residual in enumerate(
            result.residual_history,
            start=1,
        ):
            all_residual_history.append(
                {
                    "physical_step": physical_step + 1,
                    "iteration": iteration,
                    "residual": residual,
                    "temperature_residual":
                        result.temperature_residual[
                            iteration - 1
                        ],
                    "psi_residual":
                        result.psi_residual[
                            iteration - 1
                        ],
                    "voltage_residual":
                        result.voltage_residual[
                            iteration - 1
                        ],
                    "ax_residual":
                        result.ax_residual[
                            iteration - 1
                        ],
                    "ay_residual":
                        result.ay_residual[
                            iteration - 1
                        ],
                    "sigma":
                        result.sigma_history[
                            iteration - 1
                        ],
                    "oscillation":
                        result.oscillation_history[
                            iteration - 1
                        ],
                }
            )

    runtime = (
        time.perf_counter()
        - start_time
    )

    fields = simulation.fields

    return {
        "dt": dt,
        "initial_sigma": sigma,
        "final_sigma": (
            final_result.final_sigma
            if final_result is not None
            else sigma
        ),
        "ordering": ordering,
        "converged": (
            failed_steps == 0
        ),
        "physical_steps": physical_steps,
        "failed_steps": failed_steps,
        "total_iterations": total_iterations,
        "average_iterations": (
            total_iterations / physical_steps
        ),
        "final_residual": (
            final_result.final_residual
            if final_result is not None
            else float("inf")
        ),
        "runtime_seconds": runtime,
        "runtime_per_physical_step": (
            runtime / physical_steps
        ),
        "final_temperature_mean": float(
            np.mean(fields.temperature)
        ),
        "final_temperature_max": float(
            np.max(fields.temperature)
        ),
        "final_psi_amplitude_mean": float(
            np.mean(np.abs(fields.psi))
        ),
        "final_psi_amplitude_min": float(
            np.min(np.abs(fields.psi))
        ),
        "final_current_mean": float(
            np.mean(
                np.sqrt(
                    fields.current_density_x ** 2
                    + fields.current_density_y ** 2
                )
            )
        ),
        "final_heat_mean": float(
            np.mean(fields.heat_source)
        ),
        "history": all_residual_history,
    }


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_residual_history(
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

    residuals = [
        row["residual"]
        for row in history
    ]

    plt.figure()

    plt.semilogy(
        iterations,
        residuals,
        marker="o",
    )

    plt.xlabel(
        "Coupling iteration"
    )

    plt.ylabel(
        "Total residual"
    )

    plt.title(
        "Coupled convergence"
    )

    plt.grid(True)

    plt.tight_layout()

    plt.savefig(
        output_path,
        dpi=150,
    )

    plt.close()


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
            row[f"{field}_residual"]
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


# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------

def load_config(path):
    with open(
        path,
        "r",
        encoding="utf-8",
    ) as f:
        return json.load(f)


def select_orderings(config):
    mode = config.get(
        "order_mode",
        "one",
    )

    if mode == "all":
        return list(DEFAULT_ORDERINGS)

    if mode == "one":
        ordering = config[
            "single_ordering"
        ]

        if ordering not in DEFAULT_ORDERINGS:
            raise ValueError(
                f"Unknown ordering: {ordering}"
            )

        return [ordering]

    raise ValueError(
        "order_mode must be 'one' or 'all'"
    )


def save_results(
    results,
    output_directory,
    save_history,
    make_plots,
):
    output_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    compact_results = []

    for result in results:
        compact = {
            key: value
            for key, value in result.items()
            if key != "history"
        }

        compact_results.append(
            compact
        )

    with open(
        output_directory / "summary.json",
        "w",
        encoding="utf-8",
    ) as f:
        json.dump(
            compact_results,
            f,
            indent=2,
        )

    if compact_results:
        fieldnames = list(
            compact_results[0].keys()
        )

        with open(
            output_directory / "summary.csv",
            "w",
            newline="",
            encoding="utf-8",
        ) as f:
            writer = csv.DictWriter(
                f,
                fieldnames=fieldnames,
            )

            writer.writeheader()

            writer.writerows(
                compact_results
            )

    if save_history:
        with open(
            output_directory
            / "residual_history.csv",
            "w",
            newline="",
            encoding="utf-8",
        ) as f:
            fieldnames = [
                "dt",
                "initial_sigma",
                "ordering",
                "physical_step",
                "iteration",
                "residual",
                "temperature_residual",
                "psi_residual",
                "voltage_residual",
                "ax_residual",
                "ay_residual",
                "sigma",
                "oscillation",
            ]

            writer = csv.DictWriter(
                f,
                fieldnames=fieldnames,
            )

            writer.writeheader()

            for result in results:
                for row in result["history"]:
                    writer.writerow(
                        {
                            "dt": result["dt"],
                            "initial_sigma":
                                result[
                                    "initial_sigma"
                                ],
                            "ordering":
                                result[
                                    "ordering"
                                ],
                            **row,
                        }
                    )

    if make_plots:
        for index, result in enumerate(
            results,
            start=1,
        ):
            prefix = (
                f"run_{index:03d}"
            )

            plot_residual_history(
                result,
                output_directory
                / f"{prefix}_total_residual.png",
            )

            plot_field_residuals(
                result,
                output_directory
                / f"{prefix}_field_residuals.png",
            )

            plot_sigma_history(
                result,
                output_directory
                / f"{prefix}_sigma.png",
            )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "SHS coupled convergence benchmark"
        )
    )

    parser.add_argument(
        "--config",
        default=str(CONFIG_PATH),
    )

    parser.add_argument(
        "--iterations",
        type=int,
        default=None,
        help=(
            "Override maximum coupling iterations"
        ),
    )

    parser.add_argument(
        "--steps",
        type=int,
        default=None,
        help=(
            "Override physical timestep count"
        ),
    )

    parser.add_argument(
        "--orders",
        choices=("one", "all"),
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
        config["max_iterations"] = (
            args.iterations
        )

    if args.steps is not None:
        config["physical_steps"] = (
            args.steps
        )

    if args.dt is not None:
        config["dt_values"] = args.dt

    if args.sigma is not None:
        config["sigma_values"] = args.sigma

    if args.orders is not None:
        config["order_mode"] = (
            args.orders
        )

    if args.ordering is not None:
        config["single_ordering"] = (
            args.ordering
        )

    orderings = select_orderings(
        config
    )

    results = []

    for dt in config["dt_values"]:

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
                    f"max iterations = "
                    f"{config['max_iterations']}"
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
                    if result["converged"]
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
                    "runtime = "
                    f"{result['runtime_seconds']:.3f} s"
                )

    output_directory = Path(
        config["output"]["directory"]
    )

    save_results(
        results=results,
        output_directory=output_directory,
        save_history=config["output"][
            "save_history"
        ],
        make_plots=config["output"][
            "make_plots"
        ],
    )

    print()
    print(
        "Benchmark complete."
    )

    print(
        f"Results written to: "
        f"{output_directory}"
    )


if __name__ == "__main__":
    main()