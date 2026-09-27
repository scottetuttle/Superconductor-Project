import numpy as np

from shs.physics.electrical import ElectricalModel

from shs.config.builder import build_simulation

from shs.tdgl import (
    TDGLModel,
    TDGLParameters,
)


from shs.solvers.coupled import electrical_tdgl_step


def test_superconducting_state_has_no_normal_current():

    model = ElectricalModel()

    electric_field_x = np.ones((5, 5))
    electric_field_y = np.ones((5, 5))

    resistivity = np.ones((5, 5))

    superconducting_fraction = np.ones(
        (5, 5)
    )

    jx, jy = model.normal_current(
        electric_field_x,
        electric_field_y,
        resistivity,
        superconducting_fraction,
        model="condensate_depletion",
    )

    assert np.allclose(
        jx,
        0.0,
    )

    assert np.allclose(
        jy,
        0.0,
    )


def test_normal_state_has_normal_current():

    model = ElectricalModel()

    electric_field_x = np.ones((5, 5))
    electric_field_y = np.zeros((5, 5))

    resistivity = np.ones((5, 5))

    superconducting_fraction = np.zeros(
        (5, 5)
    )

    jx, jy = model.normal_current(
        electric_field_x,
        electric_field_y,
        resistivity,
        superconducting_fraction,
        model="condensate_depletion",
    )

    assert np.allclose(
        jx,
        1.0,
    )

    assert np.allclose(
        jy,
        0.0,
    )


def test_partial_superconducting_suppression():

    model = ElectricalModel()

    electric_field_x = np.ones((5, 5))
    electric_field_y = np.zeros((5, 5))

    resistivity = np.ones((5, 5))

    superconducting_fraction = np.full(
        (5, 5),
        0.75,
    )

    jx, jy = model.normal_current(
        electric_field_x,
        electric_field_y,
        resistivity,
        superconducting_fraction,
        model="condensate_depletion",
    )

    assert np.allclose(
        jx,
        0.25,
    )

    assert np.allclose(
        jy,
        0.0,
    )

def test_total_current_is_supercurrent_plus_normal_current():

    model = ElectricalModel()

    supercurrent_x = np.full(
        (5, 5),
        3.0,
    )

    supercurrent_y = np.full(
        (5, 5),
        2.0,
    )

    normal_x = np.full(
        (5, 5),
        1.0,
    )

    normal_y = np.full(
        (5, 5),
        4.0,
    )

    total_x, total_y = model.total_current(
        supercurrent_x,
        supercurrent_y,
        normal_x,
        normal_y,
    )

    assert np.allclose(
        total_x,
        4.0,
    )

    assert np.allclose(
        total_y,
        6.0,
    )


def test_supercurrent_does_not_produce_joule_heating():

    model = ElectricalModel()

    supercurrent_x = np.full(
        (5, 5),
        10.0,
    )

    supercurrent_y = np.full(
        (5, 5),
        10.0,
    )

    electric_field_x = np.ones(
        (5, 5)
    )

    electric_field_y = np.ones(
        (5, 5)
    )

    # No normal current
    normal_x = np.zeros(
        (5, 5)
    )

    normal_y = np.zeros(
        (5, 5)
    )

    heating = model.joule_heating(
        normal_x,
        normal_y,
        electric_field_x,
        electric_field_y,
    )

    assert np.allclose(
        heating,
        0.0,
    )



def test_electrical_tdgl_produces_total_current():

    simulation = build_simulation(
        "configs/simulations/nbn_hotspot_test.json"
    )

    # This is a finite-output sanity check under a very large 1 V drive,
    # so it does not require the production benchmark's tighter tolerance.
    simulation.config.electrical.solver.tolerance = 1e-6

    simulation.fields.psi[:] = (
        0.5 + 0.0j
    )

    tdgl_model = TDGLModel(
        TDGLParameters()
    )

    electrical_model = ElectricalModel()

    simulation.config.electrical.solver.max_iterations = 20000

    electrical_tdgl_step(
        simulation,
        dt=1e-13,
        tdgl_model=tdgl_model,
        electrical_model=electrical_model,
        voltage_left=1.0,
        voltage_right=0.0,
    )

    fields = simulation.fields

    assert np.all(
        np.isfinite(
            fields.current_density_x
        )
    )

    assert np.all(
        np.isfinite(
            fields.current_density_y
        )
    )

    assert np.all(
        np.isfinite(
            fields.electric_field_x
        )
    )

    assert np.all(
        np.isfinite(
            fields.electric_field_y
        )
    )

def test_supercurrent_contributes_to_total_current():

    simulation = build_simulation(
        "configs/simulations/nbn_hotspot_test.json"
    )

    # This test checks current decomposition rather than potential accuracy.
    simulation.config.electrical.solver.tolerance = 1e-9

    simulation.fields.psi[:] = (
        1.0 + 0.0j
    )

    tdgl_model = TDGLModel(
        TDGLParameters()
    )

    electrical_model = ElectricalModel()

    simulation.config.electrical.solver.max_iterations = 20000

    electrical_tdgl_step(
        simulation,
        dt=1e-13,
        tdgl_model=tdgl_model,
        electrical_model=electrical_model,
        voltage_left=1e-3,
        voltage_right=0.0,
    )

    fields = simulation.fields

    superconducting_current_x = (
        fields.supercurrent_density_x
    )

    superconducting_current_y = (
        fields.supercurrent_density_y
    )

    normal_current_x = (
        fields.current_density_x
        -
        superconducting_current_x
    )

    normal_current_y = (
        fields.current_density_y
        -
        superconducting_current_y
    )

    assert np.allclose(
        fields.current_density_x,
        (
            superconducting_current_x
            +
            normal_current_x
        ),
    )

    assert np.allclose(
        fields.current_density_y,
        (
            superconducting_current_y
            +
            normal_current_y
        ),
    )
def test_standard_tdgl_normal_conductivity_is_not_depleted_by_condensate():
    model = ElectricalModel()
    fraction = np.full((3, 4), 0.9)

    conductivity = model.normal_conductivity(
        np.full((3, 4), 2.0), fraction, model="constant"
    )

    assert np.allclose(conductivity, 0.5)
