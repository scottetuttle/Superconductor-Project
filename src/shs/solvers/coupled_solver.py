"""Transactional fixed-point coupling of one physical TDGL/electrothermal step.

Every candidate integrates from the same accepted T and psi. Trial endpoint
fields supply coupling inputs only. A failed candidate never mutates the
accepted state. Adaptive retries remain in the benchmark tools.
"""
from copy import copy, deepcopy
from dataclasses import dataclass, field
import numpy as np
from shs.physics.thermal import ThermalModel
from shs.tdgl.model import TDGLModel
from shs.tdgl.parameters import TDGLParameters
from shs.solvers.tdgl_solver import tdgl_step
from shs.solvers.thermal_solver import thermal_step
from shs.solvers.electrical_solver import electrical_simulation_step
from shs.solvers.magnetic_solver import magnetic_screening_step
from shs.numerics.convergence import ConvergenceController, ConvergenceMonitor, ConvergenceStatus
from shs.utils.defaults import default_section
from shs.optics.moving_laser import apply_laser_heat_source


_COUPLING_DEFAULTS = default_section("coupling")


@dataclass
class CouplingStepResult:
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
    simulation: object
    steps: int
    converged: bool = True
    failed_steps: int = 0
    total_coupling_iterations: int = 0
    coupling_steps: list = field(default_factory=list)
    elapsed_time: float = 0.0


def build_thermal_model(simulation):
    if simulation.config is None:
        raise ValueError('Supply a ThermalModel when config is None.')
    c = simulation.config.thermal
    return ThermalModel(c.bath_temperature, c.thermal_relaxation_rate,
                        c.max_substep, c.stability_safety_factor,
                        c.model, c.electron_heat_capacity_fraction,
                        c.electron_thermal_conductivity_fraction,
                        c.electron_phonon_coupling_W_m3_K,
                        c.phonon_escape_rate_per_s)


def build_tdgl_model(simulation):
    c = getattr(simulation.config, 'tdgl', None)
    return TDGLModel(TDGLParameters() if c is None else TDGLParameters(
        u=c.u, gamma=c.gamma, kappa=c.kappa,
        time_integrator=c.time_integrator,
        normalization=c.normalization,
        temperature_model=c.temperature_model,
        include_scalar_potential=c.include_scalar_potential,
        max_normalized_timestep=c.max_normalized_timestep,
        stability_safety_factor=c.stability_safety_factor))


def build_convergence_controller(simulation):
    c = getattr(simulation.config, 'coupling', None)
    monitor = ConvergenceMonitor(
        tolerance=getattr(c, 'tolerance', _COUPLING_DEFAULTS['tolerance']),
        prediction_window=getattr(c, 'prediction_window', _COUPLING_DEFAULTS['prediction_window']),
        reasonable_iterations=getattr(c, 'reasonable_iterations', _COUPLING_DEFAULTS['reasonable_iterations']),
        minimum_iterations=getattr(c, 'minimum_iterations', _COUPLING_DEFAULTS['minimum_iterations']),
        fast_ratio=getattr(c, 'fast_ratio', _COUPLING_DEFAULTS['fast_ratio']),
        healthy_ratio=getattr(c, 'healthy_ratio', _COUPLING_DEFAULTS['healthy_ratio']),
        stall_ratio=getattr(c, 'stall_ratio', _COUPLING_DEFAULTS['stall_ratio']),
        oscillation_window=getattr(c, 'oscillation_window', _COUPLING_DEFAULTS['oscillation_window']))
    return ConvergenceController(
        monitor,
        getattr(c, 'max_iterations', _COUPLING_DEFAULTS['max_iterations']),
        getattr(c, 'check_interval', _COUPLING_DEFAULTS['check_interval']),
    )


FIELD_SCALES = _COUPLING_DEFAULTS['field_scales']


def _relative_field_residual(old_value, new_value, scale=1.0):
    old, new = np.asarray(old_value), np.asarray(new_value)
    rms = lambda a: float(np.linalg.norm(a.ravel())/np.sqrt(a.size))
    return rms(new-old)/max(rms(old), rms(new), scale)


def _coupling_residual(previous_fields, calculated_fields, scales=None):
    scales = FIELD_SCALES if scales is None else scales
    residuals = [_relative_field_residual(getattr(previous_fields, name),
                 getattr(calculated_fields, name), scale) for name, scale in scales.items()]
    if (previous_fields.phonon_temperature is not None
            and calculated_fields.phonon_temperature is not None):
        residuals.append(_relative_field_residual(
            previous_fields.phonon_temperature, calculated_fields.phonon_temperature,
            scales.get('temperature', 1.0)))
    return max(residuals)


def _relax_fields(previous_fields, calculated_fields, relaxation):
    if not np.isfinite(relaxation) or not 0 < relaxation <= 1:
        raise ValueError('Coupling relaxation must lie in (0, 1].')
    # Derived quantities are recomputed on the next pass, never independently mixed.
    result = deepcopy(calculated_fields)
    for name in FIELD_SCALES:
        old = getattr(previous_fields, name)
        setattr(result, name, old + relaxation*(getattr(calculated_fields, name)-old))
    if (previous_fields.phonon_temperature is not None
            and calculated_fields.phonon_temperature is not None):
        result.phonon_temperature = (
            previous_fields.phonon_temperature + relaxation*(
                calculated_fields.phonon_temperature-previous_fields.phonon_temperature))
    if hasattr(result, "applied_vector_potential_x"):
        result.induced_vector_potential_x = (
            result.vector_potential_x - result.applied_vector_potential_x
        )
        result.induced_vector_potential_y = (
            result.vector_potential_y - result.applied_vector_potential_y
        )
    return result


def _execute_physics_pass(simulation, dt, thermal_model, tdgl_model,
                          voltage_left, voltage_right, include_tdgl=True,
                          accepted_fields=None, laser_time_s=0.0):
    apply_laser_heat_source(simulation, laser_time_s)
    if include_tdgl:
        tdgl_step(simulation, dt, tdgl_model,
                  initial_psi=None if accepted_fields is None else accepted_fields.psi)
    electrical_tolerance = None
    coupling = getattr(simulation.config, 'coupling', None)
    if coupling is not None:
        electrical = simulation.config.electrical
        left = electrical.voltage_left if voltage_left is None else voltage_left
        right = electrical.voltage_right if voltage_right is None else voltage_right
        voltage_scale = max(coupling.field_scales['voltage'], abs(left), abs(right))
        # Resolve each inner electrical solve more accurately than the outer
        # voltage criterion, without relaxing the user's requested tolerance.
        electrical_tolerance = min(electrical.solver.tolerance,
                                   coupling.electrical_tolerance_fraction
                                   * coupling.tolerance * voltage_scale)
    electrical_simulation_step(simulation, voltage_left, voltage_right,
                               include_superconductivity=include_tdgl,
                               solver_tolerance=electrical_tolerance)
    magnetic_screening_step(simulation)
    thermal_step(simulation, dt, thermal_model,
                 initial_temperature=None if accepted_fields is None else accepted_fields.temperature,
                 initial_phonon_temperature=None if accepted_fields is None else accepted_fields.phonon_temperature)
    return simulation.fields


def coupled_step(simulation, dt=None, thermal_model=None, tdgl_model=None,
                 voltage_left=None, voltage_right=None, relaxation=1.0,
                 time_s=0.0):
    config = simulation.config
    dt = getattr(config, 'dt', None) if dt is None else dt
    if dt is None or not np.isfinite(dt) or dt <= 0:
        raise ValueError('Supply a positive finite physical timestep.')
    if not np.isfinite(relaxation) or not 0 < relaxation <= 1:
        raise ValueError('Coupling relaxation must lie in (0, 1].')
    thermal_model = build_thermal_model(simulation) if thermal_model is None else thermal_model
    tdgl_model = build_tdgl_model(simulation) if tdgl_model is None else tdgl_model
    accepted = deepcopy(simulation.fields)
    trial = deepcopy(accepted)
    work = copy(simulation)
    controller = build_convergence_controller(simulation)
    c = getattr(config, 'coupling', None)
    scales = getattr(c, 'field_scales', FIELD_SCALES)
    legacy = config is None
    magnetic = getattr(config, "electromagnetic", None)
    screening_enabled = bool(getattr(magnetic, "include_self_field", False))
    screening_tolerance = getattr(magnetic, "screening_tolerance", np.inf)
    budget = 1 if legacy else controller.max_iterations
    if screening_enabled:
        budget = min(budget, magnetic.screening_max_iterations)
    checkpoints = 0
    for iteration in range(1, budget+1):
        work.fields = deepcopy(trial)
        calculated = _execute_physics_pass(work, dt, thermal_model, tdgl_model,
                         voltage_left, voltage_right, not legacy, accepted,
                         time_s + 0.5 * dt)
        residual = _coupling_residual(trial, calculated, scales)
        assessment = controller.record(iteration, residual)
        if controller.should_check(iteration):
            checkpoints += 1
        screening_converged = (
            not screening_enabled
            or calculated.magnetic_screening_residual <= screening_tolerance
        )
        if legacy or (
            assessment.status == ConvergenceStatus.CONVERGED and screening_converged
        ):
            # Commit atomically and preserve references to the Fields container.
            simulation.fields.__dict__.update(deepcopy(vars(calculated)))
            return simulation, CouplingStepResult(True, iteration, residual,
                ConvergenceStatus.CONVERGED, 0.0, assessment.prediction_confidence,
                False, checkpoints)
        trial = _relax_fields(trial, calculated, relaxation)
    assessment = controller.decision(budget)
    return simulation, CouplingStepResult(False, budget, residual, assessment.status,
        assessment.estimated_iterations_remaining, assessment.prediction_confidence,
        True, checkpoints)


def run_coupled_simulation(simulation, thermal_model=None, tdgl_model=None,
                           steps=None, dt=None, voltage_left=None,
                           voltage_right=None, relaxation=1.0):
    config = simulation.config
    dt = getattr(config, 'dt', None) if dt is None else dt
    if dt is None or not np.isfinite(dt) or dt <= 0:
        raise ValueError('Supply a positive finite physical timestep.')
    duration = None
    if steps is None:
        if config is None:
            raise ValueError('Supply steps when config is None.')
        duration = config.duration
        steps = int(np.ceil(duration/dt))
    if not isinstance(steps, (int, np.integer)) or steps < 1:
        raise ValueError('steps must be a positive integer.')
    thermal_model = build_thermal_model(simulation) if thermal_model is None else thermal_model
    tdgl_model = build_tdgl_model(simulation) if tdgl_model is None else tdgl_model
    results, elapsed = [], 0.0
    for _ in range(steps):
        actual_dt = dt if duration is None else min(dt, duration-elapsed)
        if actual_dt <= 0:
            break
        _, result = coupled_step(simulation, actual_dt, thermal_model, tdgl_model,
                                  voltage_left, voltage_right, relaxation,
                                  time_s=elapsed)
        results.append(result)
        if not result.converged:
            break
        elapsed += actual_dt
    failures = sum(not r.converged for r in results)
    return CoupledSolverResult(simulation, len(results), not failures, failures,
                                sum(r.iterations for r in results), results, elapsed)
