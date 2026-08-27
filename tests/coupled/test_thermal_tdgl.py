import numpy as np

from shs.config.builder import build_simulation

from shs.tdgl import (
    TDGLModel,
    TDGLParameters,
)

from shs.physics.thermal import ThermalModel

from shs.solvers.coupled import thermal_tdgl_step


def test_thermal_tdgl_step_updates_both_fields():

    simulation = build_simulation(
        "configs/simulations/nbn_hotspot_test.json"
    )

    simulation.fields.psi[:] = (
        1.0 + 0.0j
    )

    initial_temperature = (
        simulation.fields.temperature.copy()
    )

    initial_psi = (
        simulation.fields.psi.copy()
    )

    tdgl_model = TDGLModel(
        TDGLParameters()
    )

    thermal_model = ThermalModel(
        bath_temperature=3.0,
        thermal_relaxation_rate=0.0,
    )

    thermal_tdgl_step(
        simulation,
        dt=0.001,
        tdgl_model=tdgl_model,
        thermal_model=thermal_model,
    )

    assert np.all(
        np.isfinite(
            simulation.fields.psi
        )
    )

    assert np.all(
        np.isfinite(
            simulation.fields.temperature
        )
    )

    assert not np.allclose(
        simulation.fields.psi,
        initial_psi,
    )

    assert np.allclose(
        simulation.fields.temperature,
        initial_temperature,
    )