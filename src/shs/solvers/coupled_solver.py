"""
Coupled superconducting electrothermal solver.

Coordinates the SHS physics solvers.

For fully configured simulations, each physical timestep uses
an inner fixed-point iteration:

    accepted state
          |
          v
      trial state
          |
          v
        TDGL
          |
          v
      Electrical
          |
          v
      Joule heating
          |
          v
       Thermal
          |
          v
       residual
          |
          +---- not converged ----> relaxation
          |                              |
          |                              v
          +------------------------ new trial state
          |
          v
      accepted state

The individual physics equations remain implemented in their
respective physics and solver modules.

A compatibility path is retained for simulations that were
constructed without a SimulationConfig. This is required for
existing unit tests and programmatic callers that construct
Simulation objects directly.
"""

from copy import deepcopy
from dataclasses import dataclass, field

import numpy as np

from shs.config.simulation_state import Simulation

from shs.physics.thermal import ThermalModel

from shs.tdgl.model import TDGLModel
from shs.tdgl.parameters import TDGLParameters

from shs.solvers.tdgl_solver import tdgl_step
from shs.solvers.electrical_solver import electrical_step
from shs.solvers.thermal_solver import thermal_step

from shs.numerics.convergence import (
    ConvergenceController,
    ConvergenceMonitor,
    ConvergenceStatus,
)


@dataclass
class CouplingStepResult:
    """
    Diagnostics from one physical coupled timestep.
    """

    converged: bool

    iterations: int

    final_residual: float

    status: ConvergenceStatus

    estimated_iterations_remaining: float | None = None

    prediction_confidence: float = 0.0

    adaptation_recommended: bool = False

    checkpoint_count: int = 0


@dataclass
class CoupledSolverResult:
    """
    Stores coupled simulation results and convergence diagnostics.
    """

    simulation: Simulation

    steps: int

    converged: bool = True

    failed_steps: int = 0

    total_coupling_iterations: int = 0

    coupling_steps: list = field(
        default_factory=list
    )


def build_thermal_model(
    simulation: Simulation,
):
    """
    Construct a ThermalModel from SimulationConfig.

    A simulation without configuration must provide its
    ThermalModel explicitly.
    """

    if simulation.config is None:
        raise ValueError(
            "A ThermalModel must be supplied when "
            "simulation.config is None."
        )

    config = simulation.config

    return ThermalModel(
        bath_temperature=(
            config.thermal.bath_temperature
        ),
        thermal_relaxation_rate=(
            config.thermal.thermal_relaxation_rate
        ),
        max_substep=(
            config.thermal.max_substep
        ),
    )


def build_tdgl_model(
    simulation: Simulation,
):
    """
    Construct a TDGLModel from SimulationConfig.

    When no configuration exists, use the existing TDGL parameter
    defaults. This preserves compatibility with programmatically
    constructed simulations.
    """

    if simulation.config is None:
        return TDGLModel(
            parameters=TDGLParameters()
        )

    config = simulation.config

    parameters = TDGLParameters(
        u=config.tdgl.u,
        gamma=config.tdgl.gamma,
        kappa=config.tdgl.kappa,
        max_normalized_timestep=(
            config.tdgl.max_normalized_timestep
        ),
    )

    return TDGLModel(
        parameters=parameters
    )


def build_convergence_controller(
    simulation: Simulation,
):
    """
    Construct the coupled convergence controller.

    Fully configured simulations use the convergence settings from
    SimulationConfig when available.

    Missing coupling configuration falls back to the defaults used
    by the convergence infrastructure.
    """

    tolerance = 1e-8
    prediction_window = 10
    reasonable_iterations = 50
    max_iterations = 1000
    check_interval = 10

    if simulation.config is not None:

        coupling_config = getattr(
            simulation.config,
            "coupling",
            None,
        )

        if coupling_config is not None:

            tolerance = getattr(
                coupling_config,
                "tolerance",
                tolerance,
            )

            prediction_window = getattr(
                coupling_config,
                "prediction_window",
                prediction_window,
            )

            reasonable_iterations = getattr(
                coupling_config,
                "reasonable_iterations",
                reasonable_iterations,
            )

            max_iterations = getattr(
                coupling_config,
                "max_iterations",
                max_iterations,
            )

            check_interval = getattr(
                coupling_config,
                "check_interval",
                check_interval,
            )

    monitor = ConvergenceMonitor(
        tolerance=tolerance,
        prediction_window=prediction_window,
        reasonable_iterations=reasonable_iterations,
    )

    return ConvergenceController(
        monitor=monitor,
        max_iterations=max_iterations,
        check_interval=check_interval,
    )


def _relative_field_residual(
    old_value,
    new_value,
):
    """
    Calculate a normalized relative residual between two fields.

    For complex-valued fields, the complex difference is used
    directly so that phase changes are included.
    """

    old_array = np.asarray(
        old_value
    )

    new_array = np.asarray(
        new_value
    )

    difference = (
        new_array - old_array
    )

    numerator = np.linalg.norm(
        difference.ravel()
    )

    denominator = max(
        np.linalg.norm(
            old_array.ravel()
        ),
        np.linalg.norm(
            new_array.ravel()
        ),
        1.0,
    )

    return float(
        numerator / denominator
    )


def _coupling_residual(
    previous_fields,
    calculated_fields,
):
    """
    Calculate the coupled-state residual.

    Primary evolving fields are monitored rather than every derived
    diagnostic quantity.

    The maximum normalized field residual is used as the coupled
    residual.
    """

    residuals = []

    #
    # Temperature.
    #

    if hasattr(
        previous_fields,
        "temperature",
    ):

        residuals.append(
            _relative_field_residual(
                previous_fields.temperature,
                calculated_fields.temperature,
            )
        )

    #
    # Complex superconducting order parameter.
    #

    if hasattr(
        previous_fields,
        "psi",
    ):

        residuals.append(
            _relative_field_residual(
                previous_fields.psi,
                calculated_fields.psi,
            )
        )

    elif hasattr(
        previous_fields,
        "order_parameter",
    ):

        residuals.append(
            _relative_field_residual(
                previous_fields.order_parameter,
                calculated_fields.order_parameter,
            )
        )

    #
    # Electric potential / voltage.
    #

    if hasattr(
        previous_fields,
        "voltage",
    ):

        residuals.append(
            _relative_field_residual(
                previous_fields.voltage,
                calculated_fields.voltage,
            )
        )

    #
    # Vector potential.
    #

    for name in (
        "vector_potential_x",
        "vector_potential_y",
    ):

        if hasattr(
            previous_fields,
            name,
        ):

            residuals.append(
                _relative_field_residual(
                    getattr(
                        previous_fields,
                        name,
                    ),
                    getattr(
                        calculated_fields,
                        name,
                    ),
                )
            )

    if not residuals:
        raise RuntimeError(
            "No primary coupled fields are available "
            "for convergence monitoring."
        )

    return max(
        residuals
    )


def _relax_fields(
    previous_fields,
    calculated_fields,
    relaxation,
):
    """
    Apply under-relaxation to the calculated state.

    new = old + sigma * (calculated - old)
    """

    if not (
        0.0 < relaxation <= 1.0
    ):
        raise ValueError(
            "Coupling relaxation must be greater than "
            "0 and less than or equal to 1."
        )

    relaxed_fields = deepcopy(
        previous_fields
    )

    previous_values = vars(
        previous_fields
    )

    calculated_values = vars(
        calculated_fields
    )

    relaxed_values = vars(
        relaxed_fields
    )

    for name, old_value in (
        previous_values.items()
    ):

        if name not in calculated_values:
            continue

        new_value = calculated_values[
            name
        ]

        if (
            isinstance(
                old_value,
                np.ndarray,
            )
            and isinstance(
                new_value,
                np.ndarray,
            )
            and old_value.shape
            == new_value.shape
            and np.issubdtype(
                old_value.dtype,
                np.number,
            )
        ):

            relaxed_values[name] = (
                old_value
                + relaxation
                * (
                    new_value
                    - old_value
                )
            )

    return relaxed_fields


def _execute_physics_pass(
    simulation,
    dt,
    thermal_model,
    tdgl_model,
    voltage_left,
    voltage_right,
    include_tdgl=True,
):
    """
    Execute one sequential physics pass.

    The ordering is currently:

        TDGL
          |
          v
        Electrical
          |
          v
        Joule heating
          |
          v
        Thermal

    This function performs exactly one pass. It does not decide
    whether the resulting state has converged.
    """

    #
    # --------------------------------------------------
    # 1. TDGL
    # --------------------------------------------------
    #

    if include_tdgl:

        tdgl_step(
            simulation=simulation,
            dt=dt,
            tdgl_model=tdgl_model,
        )

    #
    # --------------------------------------------------
    # 2. Electrical transport
    # --------------------------------------------------
    #

    if simulation.config is not None:

        electrical_config = getattr(
            simulation.config,
            "electrical",
            None,
        )

    else:

        electrical_config = None

    if electrical_config is not None:

        solver_config = getattr(
            electrical_config,
            "solver",
            None,
        )

    else:

        solver_config = None

    if solver_config is not None:

        solver_tolerance = getattr(
            solver_config,
            "tolerance",
            1e-12,
        )

        solver_max_iterations = getattr(
            solver_config,
            "max_iterations",
            10000,
        )

        solver_omega = getattr(
            solver_config,
            "omega",
            1.7,
        )

    else:

        solver_tolerance = 1e-12
        solver_max_iterations = 10000
        solver_omega = 1.7

    #
    # The electrical solver calculates normal current,
    # total current, electric field, and Joule heating.
    #

    electrical_step(
        fields=simulation.fields,
        mesh=simulation.mesh,
        material_map=simulation.material_map,
        contact_map=simulation.contact_map,
        voltage_left=voltage_left,
        voltage_right=voltage_right,
        superconducting_current_x=(
            simulation.fields.supercurrent_density_x
        ),
        superconducting_current_y=(
            simulation.fields.supercurrent_density_y
        ),
        solver_tolerance=solver_tolerance,
        solver_max_iterations=solver_max_iterations,
        solver_omega=solver_omega,
    )

    #
    # --------------------------------------------------
    # 3. Thermal evolution
    # --------------------------------------------------
    #
    # electrical_step has already calculated Joule heating
    # and stored it in fields.heat_source.
    #

    thermal_step(
        simulation=simulation,
        dt=dt,
        thermal_model=thermal_model,
    )

    return simulation.fields


def _legacy_coupled_step(
    simulation,
    dt,
    thermal_model,
    voltage_left,
    voltage_right,
):
    """
    Execute one legacy electrothermal timestep.

    This path is used when Simulation.config is None.

    The existing electrothermal test constructs Simulation
    directly rather than through SimulationConfig. It therefore
    needs to continue using the original electrothermal sequence.

    TDGL is intentionally not called here. The legacy test uses
    dt=1e-6 seconds, which is not intended to be a TDGL timestep.
    """

    previous_fields = deepcopy(
        simulation.fields
    )

    _execute_physics_pass(
        simulation=simulation,
        dt=dt,
        thermal_model=thermal_model,
        tdgl_model=None,
        voltage_left=voltage_left,
        voltage_right=voltage_right,
        include_tdgl=False,
    )

    residual = _coupling_residual(
        previous_fields,
        simulation.fields,
    )

    return (
        simulation,
        CouplingStepResult(
            converged=True,
            iterations=1,
            final_residual=residual,
            status=ConvergenceStatus.CONVERGED,
            estimated_iterations_remaining=0.0,
            prediction_confidence=1.0,
            adaptation_recommended=False,
            checkpoint_count=0,
        ),
    )


def coupled_step(
    simulation: Simulation,
    dt: float = None,
    thermal_model: ThermalModel = None,
    tdgl_model: TDGLModel = None,
    voltage_left: float = None,
    voltage_right: float = None,
    relaxation: float = 1.0,
):
    """
    Perform one physical coupled timestep.

    Configured simulations use the self-consistent inner coupling
    iteration.

    Simulations without SimulationConfig use the legacy
    electrothermal compatibility path.
    """

    config = simulation.config

    #
    # --------------------------------------------------
    # Determine physical timestep.
    # --------------------------------------------------
    #

    if dt is None:

        if config is None:
            raise ValueError(
                "dt must be supplied when "
                "simulation.config is None."
            )

        dt = config.dt

    dt = float(
        dt
    )

    if dt <= 0.0:
        raise ValueError(
            "Coupled timestep dt must be positive."
        )

    #
    # --------------------------------------------------
    # Build models if necessary.
    # --------------------------------------------------
    #

    if thermal_model is None:

        thermal_model = build_thermal_model(
            simulation
        )

    if tdgl_model is None:

        tdgl_model = build_tdgl_model(
            simulation
        )

    #
    # --------------------------------------------------
    # Determine electrical boundary conditions.
    # --------------------------------------------------
    #

    if voltage_left is None:

        if config is not None:

            electrical_config = getattr(
                config,
                "electrical",
                None,
            )

            if electrical_config is not None:

                voltage_left = getattr(
                    electrical_config,
                    "voltage_left",
                    1e-3,
                )

            else:

                voltage_left = 1e-3

        else:

            #
            # Legacy default.
            #
            # This matches the 1 mV scale used by the NbN
            # hotspot configuration.
            #

            voltage_left = 1e-3

    if voltage_right is None:

        if config is not None:

            electrical_config = getattr(
                config,
                "electrical",
                None,
            )

            if electrical_config is not None:

                voltage_right = getattr(
                    electrical_config,
                    "voltage_right",
                    0.0,
                )

            else:

                voltage_right = 0.0

        else:

            voltage_right = 0.0

    #
    # --------------------------------------------------
    # Legacy compatibility path.
    # --------------------------------------------------
    #

    if config is None:

        return _legacy_coupled_step(
            simulation=simulation,
            dt=dt,
            thermal_model=thermal_model,
            voltage_left=voltage_left,
            voltage_right=voltage_right,
        )

    #
    # --------------------------------------------------
    # Fully configured self-consistent path.
    # --------------------------------------------------
    #

    accepted_fields = deepcopy(
        simulation.fields
    )

    trial_fields = deepcopy(
        accepted_fields
    )

    controller = build_convergence_controller(
        simulation
    )

    checkpoint_count = 0

    #
    # --------------------------------------------------
    # Inner fixed-point iteration.
    # --------------------------------------------------
    #

    for iteration in range(
        1,
        controller.max_iterations + 1,
    ):

        #
        # Start this iteration from the current trial state.
        #

        simulation.fields = deepcopy(
            trial_fields
        )

        previous_fields = deepcopy(
            simulation.fields
        )

        #
        # Execute the complete coupled physics sequence.
        #

        calculated_fields = (
            _execute_physics_pass(
                simulation=simulation,
                dt=dt,
                thermal_model=thermal_model,
                tdgl_model=tdgl_model,
                voltage_left=voltage_left,
                voltage_right=voltage_right,
                include_tdgl=True,
            )
        )

        calculated_fields = deepcopy(
            calculated_fields
        )

        #
        # Calculate residual BEFORE relaxation.
        #
        # This is important. Relaxation must not artificially
        # make the residual appear smaller.
        #

        residual = _coupling_residual(
            previous_fields,
            calculated_fields,
        )

        assessment = controller.record(
            iteration=iteration,
            residual=residual,
        )

        #
        # Convergence is always checked immediately.
        #

        if (
            assessment.status
            == ConvergenceStatus.CONVERGED
        ):

            simulation.fields = (
                calculated_fields
            )

            return (
                simulation,
                CouplingStepResult(
                    converged=True,
                    iterations=iteration,
                    final_residual=residual,
                    status=assessment.status,
                    estimated_iterations_remaining=(
                        assessment
                        .estimated_iterations_remaining
                    ),
                    prediction_confidence=(
                        assessment
                        .prediction_confidence
                    ),
                    adaptation_recommended=(
                        assessment
                        .adaptation_recommended
                    ),
                    checkpoint_count=checkpoint_count,
                ),
            )

        #
        # Strategic convergence checkpoint.
        #
        # The controller still receives every residual, but
        # should_check() identifies the iterations at which
        # higher-level adaptation will eventually occur.
        #

        if controller.should_check(
            iteration
        ):

            checkpoint_count += 1

        #
        # Generate the next trial state.
        #

        trial_fields = _relax_fields(
            previous_fields=previous_fields,
            calculated_fields=calculated_fields,
            relaxation=relaxation,
        )

    #
    # --------------------------------------------------
    # Coupling failed to converge.
    # --------------------------------------------------
    #

    #
    # Never leave the simulation at a partially converged state.
    #

    simulation.fields = deepcopy(
        accepted_fields
    )

    final_assessment = controller.decision(
        controller.max_iterations
    )

    return (
        simulation,
        CouplingStepResult(
            converged=False,
            iterations=controller.max_iterations,
            final_residual=(
                controller.monitor.residual
            ),
            status=final_assessment.status,
            estimated_iterations_remaining=(
                final_assessment
                .estimated_iterations_remaining
            ),
            prediction_confidence=(
                final_assessment
                .prediction_confidence
            ),
            adaptation_recommended=True,
            checkpoint_count=checkpoint_count,
        ),
    )


def run_coupled_simulation(
    simulation: Simulation,
    thermal_model: ThermalModel = None,
    tdgl_model: TDGLModel = None,
    steps: int = None,
    dt: float = None,
    voltage_left: float = None,
    voltage_right: float = None,
    relaxation: float = 1.0,
):
    """
    Run multiple physical coupled timesteps.

    Configured simulations use self-consistent coupling.

    Simulations without SimulationConfig use the legacy
    electrothermal compatibility path.
    """

    config = simulation.config

    #
    # --------------------------------------------------
    # Determine timestep.
    # --------------------------------------------------
    #

    if dt is None:

        if config is None:
            raise ValueError(
                "dt must be supplied when "
                "simulation.config is None."
            )

        dt = config.dt

    dt = float(
        dt
    )

    if dt <= 0.0:
        raise ValueError(
            "Coupled timestep dt must be positive."
        )

    #
    # --------------------------------------------------
    # Build models.
    # --------------------------------------------------
    #

    if thermal_model is None:

        thermal_model = build_thermal_model(
            simulation
        )

    if tdgl_model is None:

        tdgl_model = build_tdgl_model(
            simulation
        )

    #
    # --------------------------------------------------
    # Determine number of physical steps.
    # --------------------------------------------------
    #

    if steps is None:

        if config is None:
            raise ValueError(
                "steps must be supplied when "
                "simulation.config is None."
            )

        steps = max(
            1,
            int(
                round(
                    config.duration / dt
                )
            ),
        )

    steps = int(
        steps
    )

    if steps < 1:
        raise ValueError(
            "Number of coupled steps must be positive."
        )

    #
    # --------------------------------------------------
    # Run physical timestep loop.
    # --------------------------------------------------
    #

    coupling_steps = []

    failed_steps = 0

    total_iterations = 0

    all_converged = True

    for _ in range(
        steps
    ):

        (
            simulation,
            step_result,
        ) = coupled_step(
            simulation=simulation,
            dt=dt,
            thermal_model=thermal_model,
            tdgl_model=tdgl_model,
            voltage_left=voltage_left,
            voltage_right=voltage_right,
            relaxation=relaxation,
        )

        coupling_steps.append(
            step_result
        )

        total_iterations += (
            step_result.iterations
        )

        if not step_result.converged:

            failed_steps += 1

            all_converged = False

            #
            # Do not advance beyond a failed physical timestep.
            #

            break

    return CoupledSolverResult(
        simulation=simulation,
        steps=len(
            coupling_steps
        ),
        converged=all_converged,
        failed_steps=failed_steps,
        total_coupling_iterations=(
            total_iterations
        ),
        coupling_steps=coupling_steps,
    )