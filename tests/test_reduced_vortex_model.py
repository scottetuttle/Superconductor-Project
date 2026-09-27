import numpy as np

from shs.physics.reduced_vortex import (
    ReducedPinningSite, ReducedVortexMaterial, ReducedVortexState,
    reduced_vortex_forces, reduced_vortex_step,
)
from tools.reduced_vortex_diagnostics import load_config, protocol_value


def _material():
    return ReducedVortexMaterial(8.6, 40e-9, 90e-9, 15e-9)


def test_thermal_force_points_toward_hotspot():
    state = ReducedVortexState([[1e-6, 0]], [1], [0, 0], thermal_delta_K=2.0)
    force, components = reduced_vortex_forces(
        state, _material(), thermal_sigma_m=0.5e-6,
    )
    assert components["thermal"][0, 0] < 0
    assert force[0, 0] < 0
    assert force[0, 1] == 0


def test_pinning_force_points_toward_defect():
    state = ReducedVortexState([[1e-7, 0]], [1], [0, 0])
    site = ReducedPinningSite(0, 0, 1e-12, 1e-7)
    _, components = reduced_vortex_forces(
        state, _material(), thermal_sigma_m=0.5e-6, pinning_sites=[site],
    )
    assert components["pinning"][0, 0] < 0


def test_overdamped_step_moves_vortex_toward_thermal_center():
    state = ReducedVortexState([[1e-6, 0]], [1], [0, 0], thermal_delta_K=2.0)
    before = state.positions_m.copy()
    reduced_vortex_step(
        state, _material(), 1e-10, commanded_center_m=[0, 0],
        commanded_delta_K=2.0, thermal_sigma_m=0.5e-6,
        thermal_response_time_s=1e-9,
    )
    assert state.positions_m[0, 0] < before[0, 0]


def test_manipulation_configuration_and_protocol_are_valid():
    config = load_config("tools/config/reduced_optical_tweezers.json")
    center, delta = protocol_value(config["laser"]["protocol"], 7.03e-6)
    assert 8e-6 < center[0] < 22e-6
    assert np.isclose(delta, 3.0)
