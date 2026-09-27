import numpy as np

from shs.config.builder import build_simulation
from shs.physics.josephson import junction_observables, weak_link_suppression


def test_configured_weak_link_is_localized_and_measurable():
    simulation = build_simulation("configs/simulations/nbn_tdgl_weak_link.json")
    suppression = weak_link_suppression(simulation)
    assert suppression.max() == 0.18
    assert np.any(suppression > 0)
    assert np.all(suppression[:, 0] == 0)
    values = junction_observables(simulation)
    assert values["junction_voltage_V"] == 0
    assert values["junction_phase_difference_rad"] == 0
    assert np.isfinite(values["junction_mean_amplitude"])

