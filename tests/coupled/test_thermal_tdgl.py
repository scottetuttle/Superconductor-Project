import numpy as np

from shs.config.builder import build_simulation

from shs.tdgl import (
    TDGLModel,
    TDGLParameters,
)

from shs.physics.thermal import ThermalModel

from shs.solvers.coupled import thermal_tdgl_step

from shs.physics.electrical import ElectricalModel

def test_tdgl_suppression_generates_joule_heating():

    simulation = build_simulation(
        "configs/simulations/nbn_hotspot_test.json"
    )

    # Suppressed superconductivity
    simulation.fields.psi[:] = (
        0.5 + 0.0j
    )

    # Apply an electric field
    simulation.fields.electric_field_x[:] = 1.0

    simulation.fields.electric_field_y[:] = 0.0

    tdgl_model = TDGLModel(
        TDGLParameters()
    )

    thermal_model = ThermalModel(
        bath_temperature=3.0,
        thermal_relaxation_rate=0.0,
    )

    electrical_model = ElectricalModel()

    thermal_tdgl_step(
        simulation,
        dt=1e-13,
        tdgl_model=tdgl_model,
        thermal_model=thermal_model,
        electrical_model=electrical_model,
    )

    assert np.all(simulation.fields.joule_heat_source >= 0.0)
    assert np.any(simulation.fields.joule_heat_source > 0.0)
def test_supercurrent_is_not_dissipative():
    # Compare electrical solves with identical bias/conductivity and different
    # divergence-free supercurrents; only the normal current contributes heat.
    from shs.solvers.electrical_solver import electrical_step
    simulation = build_simulation("configs/simulations/nbn_hotspot_test.json")
    fields = simulation.fields
    args = dict(fields=fields, mesh=simulation.mesh,
                material_map=simulation.material_map, contact_map=simulation.contact_map,
                solver_tolerance=1e-9, solver_max_iterations=20000)
    electrical_step(**args)
    # Closed, divergence-free link loop wholly inside the film.
    jx, jy = np.zeros_like(fields.voltage), np.zeros_like(fields.voltage)
    jx[40, 40], jx[41, 40] = 10, -10
    jy[40, 41], jy[40, 40] = 10, -10
    electrical_step(**args, superconducting_current_x=jx, superconducting_current_y=jy)
    expected_heat = (
        fields.normal_current_density_x * fields.electric_field_x
        + fields.normal_current_density_y * fields.electric_field_y
    )
    assert np.allclose(fields.joule_heat_source, expected_heat)
    assert np.isclose(
        fields.current_density_x[40, 40]
        - fields.normal_current_density_x[40, 40],
        10,
    )


def test_electrothermal_feedback():

    simulation = build_simulation(
        "configs/simulations/nbn_hotspot_test.json"
    )

    simulation.boundaries = None

    # Start with a superconducting state.
    simulation.fields.psi[:] = (
        1.0 + 0.0j
    )

    # Apply an electric field so that the normal
    # component produces dissipative heating.
    simulation.fields.electric_field_x[:] = 1.0
    simulation.fields.electric_field_y[:] = 0.0

    initial_temperature = np.mean(
        simulation.fields.temperature
    )

    initial_amplitude = np.mean(
        np.abs(simulation.fields.psi)
    )

    tdgl_model = TDGLModel(
        TDGLParameters()
    )

    thermal_model = ThermalModel(
        bath_temperature=3.0,
        thermal_relaxation_rate=0.0,
    )

    electrical_model = ElectricalModel()

    # Advance the coupled system.
    for _ in range(100):

        thermal_tdgl_step(
            simulation,
            dt=1e-13,
            tdgl_model=tdgl_model,
            thermal_model=thermal_model,
            electrical_model=electrical_model,
        )

    final_temperature = np.mean(
        simulation.fields.temperature
    )

    final_amplitude = np.mean(
        np.abs(simulation.fields.psi)
    )

    # Electrical dissipation must heat the system.
    assert final_temperature > initial_temperature

    # The resulting thermal increase must suppress
    # the superconducting order parameter.
    assert final_amplitude < initial_amplitude
