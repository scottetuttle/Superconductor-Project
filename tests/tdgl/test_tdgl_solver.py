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

from shs.materials.database import get_material

from shs.config.simulation_state import Simulation

from shs.tdgl import (
    TDGLModel,
    TDGLParameters,
)

from shs.solvers import tdgl_step


def create_nbn_simulation(
    temperature=3.0,
):
    """
    Construct the standard NbN TDGL test simulation.
    """

    geometry = load_geometry(
        "configs/geometry/NbN_film.json"
    )

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

    simulation = Simulation(
        config=None,
        geometry=geometry,
        mesh=mesh,
        region_map=region_map,
        material_map=material_map,
        boundaries=None,
        fields=fields,
        contact_map=build_contact_map(
            geometry,
            mesh
        ),
    )

    return simulation


def test_tdgl_low_temperature_stability():

    simulation = create_nbn_simulation(
        temperature=3.0
    )

    model = TDGLModel(
        TDGLParameters()
    )

    simulation.fields.psi[:] = (
        1.0 + 0.0j
    )

    initial = np.abs(
        simulation.fields.psi.copy()
    )

    tdgl_step(
        simulation,
        dt=0.001,
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
        TDGLParameters()
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
            dt=dt,
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
        TDGLParameters()
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
            dt=0.001,
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
        TDGLParameters()
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
        TDGLParameters()
    )

    simulation.fields.psi[:] = (
        1.0 + 0.0j
    )

    dt = 0.001
    num_steps = 100000

    for _ in range(num_steps):

        tdgl_step(
            simulation,
            dt=dt,
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