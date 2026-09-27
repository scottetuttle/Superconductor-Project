import numpy as np

from shs.geometry import (
    load_geometry,
    create_mesh,
)

from shs.physics import Fields

from shs.mapping import (
    build_region_map,
    build_material_map,
    build_contact_map,
)
from shs.mapping.contact_map import ContactMap

from shs.materials.database import get_material

from shs.config.simulation_state import Simulation

from shs.tdgl import (
    TDGLModel,
    TDGLParameters,
    TDGLBoundarySet,
    TDGLBoundaryCondition,
    TDGLBoundarySide,
    TDGLBoundaryType,
)

from shs.solvers import tdgl_step
from shs.solvers.tdgl_solver import generalized_tdgl_update, tdgl_scales


def create_nbn_simulation(
    temperature=3.0,
):
    """
    Construct the standard NbN TDGL test simulation.
    """

    geometry = load_geometry(
        "configs/geometry/NbN_film.json"
    )

    geometry.contacts = [c for c in geometry.contacts if c.contact_type == 'current']
    geometry.film.nx = geometry.film.ny = 8

    mesh = create_mesh(
        geometry
    )

    material = get_material(
        "NbN"
    )

    region_map = build_region_map(
        geometry,
        mesh
    )

    material_map = build_material_map(
        region_map,
        material
    )

    fields = Fields.create(
        mesh,
        initial_temperature=temperature
    )

    tdgl_boundaries = TDGLBoundarySet()

    for side in TDGLBoundarySide:

        tdgl_boundaries.add(
            TDGLBoundaryCondition(
                side=side,
                type=TDGLBoundaryType.INSULATING,
        )
    )

    left = np.zeros((mesh.ny, mesh.nx), dtype=bool)
    right = np.zeros((mesh.ny, mesh.nx), dtype=bool)
    left[:, 0] = True
    right[:, -1] = True

    simulation = Simulation(
    config=None,
    geometry=geometry,
    mesh=mesh,
    region_map=region_map,
    material_map=material_map,
    boundaries=None,
    tdgl_boundaries=tdgl_boundaries,
    fields=fields,
    contact_map=ContactMap(
        {"left_current": left, "right_current": right},
        {"left_current": "current", "right_current": "current"},
    ),
)


    return simulation


def test_scalar_potential_rotates_phase_without_changing_equilibrium_amplitude():
    simulation = create_nbn_simulation(temperature=3.0)
    model = TDGLModel(TDGLParameters(include_scalar_potential=True, normalization="legacy_gl", temperature_model="one_minus_t_over_tc"))
    reduced_temperature = 3.0 / simulation.material_map.materials[0].Tc
    amplitude = model.equilibrium_amplitude(reduced_temperature)
    simulation.fields.psi.fill(amplitude + 0j)
    simulation.fields.voltage.fill(1e-3)
    scales = tdgl_scales(simulation, model)
    normalized_dt = 1e-3
    tdgl_step(simulation, normalized_dt * scales.time_scale, model)
    expected_phase = -float(
        scales.scalar_potential_to_dimensionless(1e-3)
    ) * normalized_dt
    assert np.allclose(np.abs(simulation.fields.psi), amplitude, atol=1e-12)
    assert np.allclose(np.angle(simulation.fields.psi), expected_phase, atol=1e-12)


def test_constant_voltage_gauge_offset_only_changes_global_phase():
    first = create_nbn_simulation(temperature=3.0)
    second = create_nbn_simulation(temperature=3.0)
    model = TDGLModel(TDGLParameters(include_scalar_potential=True, normalization="legacy_gl", temperature_model="one_minus_t_over_tc"))
    first.fields.voltage.fill(0.0)
    second.fields.voltage.fill(2e-4)
    scales = tdgl_scales(first, model)
    dt = 1e-3 * scales.time_scale
    tdgl_step(first, dt, model)
    tdgl_step(second, dt, model)
    phase = -float(scales.scalar_potential_to_dimensionless(2e-4)) * 1e-3
    assert np.allclose(second.fields.psi, first.fields.psi * np.exp(1j * phase))
    assert np.allclose(
        second.fields.supercurrent_density_x,
        first.fields.supercurrent_density_x,
    )


def test_generalized_update_reduces_to_phase_rotated_euler_at_zero_gamma():
    psi = np.array([[0.6 + 0.2j]])
    rhs = np.array([[0.1 - 0.3j]])
    potential = np.array([[0.4]])
    dt, u = 1e-3, 5.79
    result = generalized_tdgl_update(psi, rhs, potential, dt, u, 0.0)
    expected = np.exp(-1j * potential * dt) * (psi + dt * rhs / u)
    assert np.allclose(result, expected)


def test_nonzero_gamma_generalized_update_is_finite_for_small_step():
    psi = np.full((3, 4), 0.7 + 0.1j)
    rhs = np.full((3, 4), -0.2 + 0.05j)
    result = generalized_tdgl_update(psi, rhs, np.zeros((3, 4)), 1e-3, 5.79, 2.0)
    assert np.all(np.isfinite(result))
    assert np.all(np.abs(result) < np.abs(psi))


def test_tdgl_low_temperature_stability():

    simulation = create_nbn_simulation(
        temperature=3.0
    )

    model = TDGLModel(
        TDGLParameters(normalization="legacy_gl", temperature_model="one_minus_t_over_tc")
    )

    simulation.fields.psi[:] = (
        1.0 + 0.0j
    )

    initial = np.abs(
        simulation.fields.psi.copy()
    )

    tdgl_step(
        simulation,
        dt=0.001 * model.characteristic_time(simulation.material_map.materials[0].Tc),
        tdgl_model=model,
    )

    final = np.abs(
        simulation.fields.psi
    )

    assert np.all(
        np.isfinite(final)
    )

    assert np.all(
        final >= 0.0
    )

    assert np.max(
        np.abs(final - initial)
    ) < 1e-2


def test_tdgl_uniform_state_matches_analytical_solution():

    temperature = 3.0

    simulation = create_nbn_simulation(
        temperature=temperature
    )

    model = TDGLModel(
        TDGLParameters(normalization="legacy_gl", temperature_model="one_minus_t_over_tc")
    )

    simulation.fields.psi[:] = (
        1.0 + 0.0j
    )

    initial_amplitude = np.mean(
        np.abs(
            simulation.fields.psi
        )
    )

    dt = 0.001
    num_steps = 1000

    for _ in range(num_steps):

        tdgl_step(
            simulation,
            dt=dt * model.characteristic_time(simulation.material_map.materials[0].Tc),
            tdgl_model=model,
        )

    final_amplitude = np.mean(
        np.abs(
            simulation.fields.psi
        )
    )

    #
    # Total dimensionless simulation time.
    #

    total_time = (
        num_steps *
        dt
    )

    #
    # Reduced temperature.
    #

    reduced_temperature = (
        temperature /
        simulation.material_map.Tc[0, 0]
    )

    a = (
        1.0 -
        reduced_temperature
    )

    #
    # Analytical solution for uniform TDGL.
    #
    # u dpsi/dt = a psi - psi^3
    #

    expected_amplitude = np.sqrt(
        a /
        (
            1.0 +
            (
                a /
                initial_amplitude**2
                -
                1.0
            ) *
            np.exp(
                -2.0 *
                a *
                total_time /
                model.parameters.u
            )
        )
    )

    print(
        "Initial |psi|:",
        initial_amplitude
    )

    print(
        "Final |psi|:",
        final_amplitude
    )

    print(
        "Analytical |psi|:",
        expected_amplitude
    )

    print(
        "Absolute error:",
        abs(
            final_amplitude -
            expected_amplitude
        )
    )

    #
    # The numerical solution should closely
    # reproduce the analytical transient.
    #

    assert np.isclose(
        final_amplitude,
        expected_amplitude,
        atol=5e-4,
    )

    #
    # The system must also be relaxing downward
    # from its initial state.
    #

    assert (
        final_amplitude <
        initial_amplitude
    )

def test_tdgl_normal_state_above_tc():

    simulation = create_nbn_simulation(
        temperature=1.1 *
        get_material("NbN").Tc
    )

    model = TDGLModel(
        TDGLParameters(normalization="legacy_gl", temperature_model="one_minus_t_over_tc")
    )

    simulation.fields.psi[:] = (
        1.0 + 0.0j
    )

    initial_amplitude = np.mean(
        np.abs(
            simulation.fields.psi
        )
    )

    for _ in range(2000):

        tdgl_step(
            simulation,
            dt=0.001 * model.characteristic_time(simulation.material_map.materials[0].Tc),
            tdgl_model=model,
        )

    final_amplitude = np.mean(
        np.abs(
            simulation.fields.psi
        )
    )

    #
    # Above Tc, the superconducting state is
    # unstable and the order parameter must decay.
    #

    assert final_amplitude < initial_amplitude

def test_tdgl_equilibrium_above_tc():

    model = TDGLModel(
        TDGLParameters(normalization="legacy_gl", temperature_model="one_minus_t_over_tc")
    )

    equilibrium = model.equilibrium_amplitude(
        1.1
    )

    assert np.isclose(
        equilibrium,
        0.0
    )

def test_tdgl_uniform_state_approaches_equilibrium():

    temperature = 3.0

    simulation = create_nbn_simulation(
        temperature=temperature
    )

    model = TDGLModel(
        TDGLParameters(normalization="legacy_gl", temperature_model="one_minus_t_over_tc")
    )

    simulation.fields.psi[:] = (
        1.0 + 0.0j
    )

    dt = 0.001
    dt = 0.01
    num_steps = 10000

    for _ in range(num_steps):

        tdgl_step(
            simulation,
            dt=dt * model.characteristic_time(simulation.material_map.materials[0].Tc),
            tdgl_model=model,
        )

    final_amplitude = np.mean(
        np.abs(
            simulation.fields.psi
        )
    )

    reduced_temperature = (
        temperature /
        simulation.material_map.Tc[0, 0]
    )

    expected_amplitude = np.sqrt(
        max(
            1.0 -
            reduced_temperature,
            0.0,
        )
    )

    print(
        "Final |psi|:",
        final_amplitude
    )

    print(
        "Equilibrium |psi|:",
        expected_amplitude
    )

    assert np.isclose(
        final_amplitude,
        expected_amplitude,
        atol=5e-4,
    )

def test_tdgl_solver_enforces_insulating_boundaries():

    simulation = create_nbn_simulation(
        temperature=3.0
    )

    model = TDGLModel(
        TDGLParameters(normalization="legacy_gl", temperature_model="one_minus_t_over_tc")
    )

    simulation.fields.psi[:] = (
        1.0 + 0.0j
    )

    simulation.fields.psi[:, 1] = (
        0.8 + 0.2j
    )

    tdgl_step(
        simulation,
        dt=0.001 * model.characteristic_time(simulation.material_map.materials[0].Tc),
        tdgl_model=model,
    )

    psi = simulation.fields.psi

    Ax = simulation.fields.vector_potential_x
    Ay = simulation.fields.vector_potential_y

    Dx_left = (
        (psi[:, 1] - psi[:, 0])
        / (simulation.mesh.dx / simulation.material_map.materials[0].coherence_length)
        - 1j * Ax[:, 0] * psi[:, 0]
    )

    assert np.allclose(
        Dx_left,
        0.0,
        atol=1e-10,
    )


def test_tdgl_solver_enforces_normal_contact_boundaries():
    simulation = create_nbn_simulation(temperature=3.0)
    for side in (TDGLBoundarySide.LEFT, TDGLBoundarySide.RIGHT):
        simulation.tdgl_boundaries.add(
            TDGLBoundaryCondition(side=side, type=TDGLBoundaryType.NORMAL_CONTACT)
        )

    model = TDGLModel(TDGLParameters(normalization="legacy_gl", temperature_model="one_minus_t_over_tc"))
    simulation.fields.psi.fill(0.8 + 0.1j)
    tdgl_step(
        simulation,
        dt=0.001 * model.characteristic_time(
            simulation.material_map.materials[0].Tc
        ),
        tdgl_model=model,
    )

    assert np.all(simulation.fields.psi[:, 0] == 0.0)
    assert np.all(simulation.fields.psi[:, -1] == 0.0)
    assert np.any(np.abs(simulation.fields.psi[:, 1:-1]) > 0.0)


def test_tdgl_solver_uses_partial_contact_mask_and_insulates_uncovered_edge():
    simulation = create_nbn_simulation(temperature=3.0)
    partial = np.zeros_like(simulation.fields.psi, dtype=bool)
    partial[2:6, 0] = True
    simulation.contact_map.contact_masks["left_current"] = partial
    simulation.tdgl_boundaries.add(TDGLBoundaryCondition(
        side=TDGLBoundarySide.LEFT,
        type=TDGLBoundaryType.NORMAL_CONTACT_MASK,
    ))
    model = TDGLModel(TDGLParameters(
        normalization="legacy_gl",
        temperature_model="one_minus_t_over_tc",
    ))
    simulation.fields.psi.fill(0.8 + 0.1j)

    tdgl_step(
        simulation,
        dt=0.001 * model.characteristic_time(
            simulation.material_map.materials[0].Tc
        ),
        tdgl_model=model,
    )

    psi = simulation.fields.psi
    assert np.all(psi[2:6, 0] == 0.0)
    assert np.allclose(psi[[0, 1, 6, 7], 0], psi[[0, 1, 6, 7], 1])
