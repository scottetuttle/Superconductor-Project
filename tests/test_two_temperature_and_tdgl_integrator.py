"""Conservation and temporal-order checks for the optional baseline models."""
from copy import deepcopy

import numpy as np
import pytest

from shs.config.builder import build_simulation
from shs.geometry.mesh import create_mesh
from shs.mapping import build_region_map, build_material_map
from shs.mapping.contact_map import ContactMap
from shs.physics.fields import Fields
from shs.physics.thermal import ThermalModel
from shs.solvers.thermal_solver import thermal_step
from shs.solvers.tdgl_solver import tdgl_step, tdgl_scales
from shs.solvers.coupled_solver import build_tdgl_model
from shs.solvers.coupled_solver import coupled_step


@pytest.fixture
def small_film():
    sim = build_simulation("configs/simulations/nbn_tdgl_transport.json")
    sim.geometry.film.nx = sim.geometry.film.ny = 9
    sim.mesh = create_mesh(sim.geometry)
    sim.region_map = build_region_map(sim.geometry, sim.mesh)
    sim.material_map = build_material_map(sim.region_map, sim.material_map.materials[0])
    left = np.zeros((9, 9), dtype=bool)
    right = np.zeros((9, 9), dtype=bool)
    left[:, 0], right[:, -1] = True, True
    sim.contact_map = ContactMap(
        {"left_current": left, "right_current": right},
        {"left_current": "current", "right_current": "current"})
    sim.fields = Fields.create(sim.mesh, 4.0)
    sim.boundaries = None
    return sim


def test_electron_phonon_exchange_conserves_energy_and_reduces_difference(small_film):
    sim = small_film
    sim.material_map.thermal_conductivity.fill(0)
    sim.fields.temperature.fill(6)
    sim.fields.phonon_temperature.fill(4)
    model = ThermalModel(4, model="two_temperature", max_substep=1e-10,
                         electron_heat_capacity_fraction=.3,
                         electron_phonon_coupling_W_m3_K=1e15)
    ce = .3*sim.material_map.heat_capacity
    cp = .7*sim.material_map.heat_capacity
    before = ce*sim.fields.temperature + cp*sim.fields.phonon_temperature
    thermal_step(sim, 2e-10, model)
    after = ce*sim.fields.temperature + cp*sim.fields.phonon_temperature
    assert np.allclose(after, before, rtol=1e-14)
    assert np.all(sim.fields.temperature > sim.fields.phonon_temperature)
    assert np.max(sim.fields.temperature-sim.fields.phonon_temperature) < 2


def test_optical_energy_enters_electrons_then_flows_to_phonons(small_film):
    sim = small_film
    sim.material_map.thermal_conductivity.fill(0)
    sim.fields.laser_heat_source = np.full_like(sim.fields.temperature, 1e15)
    sim.fields.external_heat_source = np.zeros_like(sim.fields.temperature)
    model = ThermalModel(4, model="two_temperature", max_substep=1e-11,
                         electron_phonon_coupling_W_m3_K=1e15)
    thermal_step(sim, 1e-11, model)
    assert np.all(sim.fields.temperature > 4)
    assert np.all(sim.fields.phonon_temperature > 4)
    total = (.5*sim.material_map.heat_capacity*(sim.fields.temperature-4)
             + .5*sim.material_map.heat_capacity*(sim.fields.phonon_temperature-4))
    assert np.allclose(total, 1e4, rtol=1e-13)


def test_heun_reduces_tdgl_temporal_error_for_uniform_relaxation(small_film):
    sim = small_film
    sim.fields.psi.fill(.5+0j)
    sim.fields.temperature.fill(4)
    sim.config.tdgl.temperature_model = "one_minus_t_over_tc"
    model = build_tdgl_model(sim)
    scale = tdgl_scales(sim, model)
    dt = .002*scale.time_scale

    def evolve(method, count):
        work = deepcopy(sim)
        local_model = build_tdgl_model(work)
        local_model.parameters.time_integrator = method
        for _ in range(count):
            tdgl_step(work, dt*4/count, local_model)
        return work.fields.psi[4, 4]

    reference = evolve("heun", 256)
    euler = evolve("euler", 4)
    heun = evolve("heun", 4)
    assert abs(heun-reference) < abs(euler-reference)/10


def test_two_temperature_coupled_step_is_transactional(small_film):
    sim = small_film
    sim.config.thermal.model = "two_temperature"
    sim.config.thermal.electron_phonon_coupling_W_m3_K = 1e15
    sim.config.thermal.phonon_escape_rate_per_s = 1e10
    sim.config.thermal.thermal_relaxation_rate = 0
    sim.config.electrical.drive_mode = "voltage"
    sim.config.electrical.voltage_left = 0
    sim.config.electrical.voltage_right = 0
    sim.fields.temperature.fill(5)
    sim.fields.phonon_temperature.fill(4)
    initial = sim.fields.phonon_temperature.copy()
    _, result = coupled_step(sim, dt=1e-14)
    assert result.converged
    assert np.all(sim.fields.phonon_temperature > initial)
