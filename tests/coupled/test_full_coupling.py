import numpy as np

from shs.config.builder import build_simulation

from shs.solvers.electrical_solver import electrical_step
from shs.tdgl import (
    TDGLModel,
    TDGLParameters,
)

from shs.physics.thermal import ThermalModel
from shs.physics.electrical import ElectricalModel

from shs.solvers.coupled import (
    thermal_tdgl_step,
)
import matplotlib.pyplot as plt




def test_heating_suppresses_superconductivity():

    simulation = build_simulation(
        "configs/simulations/nbn_hotspot_test.json"
    )

    simulation.fields.psi[:] = 1.0 + 0.0j

    tdgl_model = TDGLModel(
        TDGLParameters()
    )

    thermal_model = ThermalModel(
        bath_temperature=3.0,
        thermal_relaxation_rate=0.0,
    )

    electrical_model = ElectricalModel()

    initial_amplitude = np.mean(
        np.abs(simulation.fields.psi)
    )

    for _ in range(100):

        thermal_tdgl_step(
            simulation,
            dt=1e-13,
            tdgl_model=tdgl_model,
            thermal_model=thermal_model,
            electrical_model=electrical_model,
        )

        fields = simulation.fields

        assert np.all(
            np.isfinite(fields.psi)
        )

        assert np.all(
            np.isfinite(fields.temperature)
        )

        assert np.all(
            np.isfinite(fields.voltage)
        )

        assert np.all(
            np.isfinite(fields.heat_source)
        )

    final_amplitude = np.mean(
        np.abs(simulation.fields.psi)
    )

    left = simulation.contact_map.contact_masks[
        "left_current"
    ]

    right = simulation.contact_map.contact_masks[
        "right_current"
    ]

    left_voltage = np.mean(
        simulation.fields.voltage[left]
    )

    right_voltage = np.mean(
        simulation.fields.voltage[right]
    )

    voltage_drop = (
        left_voltage - right_voltage
    )

    # The superconducting order parameter was suppressed.
    assert final_amplitude < initial_amplitude

    # The electrical solution retained the imposed voltage.
    assert np.isclose(
        voltage_drop,
        1e-3,
        rtol=1e-6,
        atol=1e-9,
    )

    # Dissipative heating is present.
    assert np.any(
        simulation.fields.heat_source > 0.0
    )


#    assert final_temperature > initial_temperature

#    assert final_amplitude < initial_amplitude


def test_suppressed_superconductivity_increases_normal_current():

    simulation = build_simulation(
        "configs/simulations/nbn_hotspot_test.json"
    )

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

    # Strongly superconducting state

    simulation.fields.psi[:] = (
        1.0 + 0.0j
    )

    thermal_tdgl_step(
        simulation,
        dt=1e-13,
        tdgl_model=tdgl_model,
        thermal_model=thermal_model,
        electrical_model=electrical_model,
    )

    superconducting_current = (
        simulation.fields.current_density_x.copy()
    )

    # Suppressed superconducting state

    simulation.fields.psi[:] = (
        0.5 + 0.0j
    )

    thermal_tdgl_step(
        simulation,
        dt=1e-13,
        tdgl_model=tdgl_model,
        thermal_model=thermal_model,
        electrical_model=electrical_model,
    )

    suppressed_current = (
        simulation.fields.current_density_x.copy()
    )

    assert np.mean(
        np.abs(suppressed_current)
    ) > np.mean(
        np.abs(superconducting_current)
    )


def test_full_thermal_electrical_tdgl_feedback():

    simulation = build_simulation(
        "configs/simulations/nbn_hotspot_test.json"
    )

    simulation.fields.psi[:] = (
        1.0 + 0.0j
    )

    center_y = simulation.mesh.ny // 2
    center_x = simulation.mesh.nx // 2

    simulation.fields.electric_field_x[
        center_y - 2:center_y + 3,
        center_x - 2:center_x + 3,
    ] = 1.0

    initial_temperature = (
        simulation.fields.temperature.copy()
    )

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

    # Heating occurred where the electric field was applied.
    heated_region = (
        simulation.fields.temperature[
            center_y - 2:center_y + 3,
            center_x - 2:center_x + 3,
        ]
    )

    initial_region = (
        initial_temperature[
            center_y - 2:center_y + 3,
            center_x - 2:center_x + 3,
        ]
    )

    assert np.all(
        heated_region >= initial_region
    )

    # Dissipative heating exists.
    assert np.any(
        simulation.fields.heat_source > 0.0
    )

    # The coupled solver produced a finite TDGL state.
    assert np.all(
        np.isfinite(simulation.fields.psi)
    )

    # Temperature remains finite.
    assert np.all(
        np.isfinite(simulation.fields.temperature)
    )

def test_electrical_conductivity_is_suppressed_by_superconductivity():

    simulation = build_simulation(
        "configs/simulations/nbn_hotspot_test.json"
    )

    simulation.fields.psi[:] = 0.5 + 0j

    simulation.fields.electric_field_x[:] = 1.0
    simulation.fields.electric_field_y[:] = 0.0

    electrical_model = ElectricalModel()

    Jx, Jy = electrical_model.normal_current(
        simulation.fields.electric_field_x,
        simulation.fields.electric_field_y,
        simulation.material_map.normal_resistivity,
        np.abs(simulation.fields.psi) ** 2,
    )

    expected_fraction = 1.0 - 0.25

    expected = (
        expected_fraction
        *
        simulation.material_map.electrical_conductivity
    )

    assert np.allclose(
        Jx,
        expected,
    )

    assert np.allclose(
        Jy,
        0.0,
    )

def test_electrical_solver_separates_supercurrent_and_normal_current():

    simulation = build_simulation(
        "configs/simulations/nbn_hotspot_test.json"
    )

    simulation.fields.psi[:] = 0.5 + 0j

    simulation.fields.supercurrent_density_x[:] = 10.0
    simulation.fields.supercurrent_density_y[:] = 0.0

    electrical_step(
        fields=simulation.fields,
        mesh=simulation.mesh,
        material_map=simulation.material_map,
        contact_map=simulation.contact_map,
        superconducting_current_x=(
            simulation.fields.supercurrent_density_x
        ),
        superconducting_current_y=(
            simulation.fields.supercurrent_density_y
        ),
        superconducting_fraction=(
            np.abs(simulation.fields.psi) ** 2
        ),
    )

    assert np.allclose(
        simulation.fields.current_density_x,
        (
            simulation.fields.normal_current_density_x
            +
            simulation.fields.supercurrent_density_x
        ),
    )

def test_electrical_solver_joule_heating_uses_normal_current():

    simulation = build_simulation(
        "configs/simulations/nbn_hotspot_test.json"
    )

    simulation.fields.psi[:] = 0.5 + 0j

    simulation.fields.supercurrent_density_x[:] = 100.0
    simulation.fields.supercurrent_density_y[:] = 0.0

    electrical_step(
        fields=simulation.fields,
        mesh=simulation.mesh,
        material_map=simulation.material_map,
        contact_map=simulation.contact_map,
        superconducting_current_x=(
            simulation.fields.supercurrent_density_x
        ),
        superconducting_current_y=(
            simulation.fields.supercurrent_density_y
        ),
        superconducting_fraction=(
            np.abs(simulation.fields.psi) ** 2
        ),
    )

    expected_heat = (
        simulation.fields.normal_current_density_x
        *
        simulation.fields.electric_field_x
        +
        simulation.fields.normal_current_density_y
        *
        simulation.fields.electric_field_y
    )

    assert np.allclose(
        simulation.fields.heat_source,
        expected_heat,
    )