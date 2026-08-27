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

    assert np.all(
        simulation.fields.heat_source > 0.0
    )
def test_supercurrent_is_not_dissipative():

    simulation = build_simulation(
        "configs/simulations/nbn_hotspot_test.json"
    )

    simulation.fields.psi[:] = (
        1.0 + 0.0j
    )

    simulation.fields.supercurrent_density_x[:] = 10.0
    simulation.fields.supercurrent_density_y[:] = 0.0

    simulation.fields.electric_field_x[:] = 0.0
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
    print(
        "mean |psi|^2:",
        np.mean(np.abs(simulation.fields.psi) ** 2)
    )

    print(
        "mean normal current:",
        np.mean(
            np.sqrt(
                simulation.fields.current_density_x**2 +
                simulation.fields.current_density_y**2
            )
        )
    )

    print(
        "mean heat:",
        np.mean(simulation.fields.heat_source)
    )
    assert np.allclose(
        simulation.fields.heat_source,
        0.0,
    )

def test_electrothermal_feedback():

    simulation = build_simulation(
        "configs/simulations/nbn_hotspot_test.json"
    )

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