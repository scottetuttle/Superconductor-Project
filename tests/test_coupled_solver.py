from shs.solvers import run_coupled_simulation

from shs.physics import Fields

from shs.geometry.database import load_geometry
from shs.geometry.mesh import create_mesh

from shs.mapping import (
    build_region_map,
    build_material_map,
    build_contact_map,
)

from shs.materials.database import get_material

from shs.physics.thermal import ThermalModel


from shs.config.simulation_state import Simulation



def test_coupled_electrothermal():

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


    contact_map = build_contact_map(
        geometry,
        mesh
    )


    fields = Fields.create(
        mesh,
        initial_temperature=3.0
    )


    simulation = Simulation(
        config=None,
        geometry=geometry,
        mesh=mesh,
        region_map=region_map,
        material_map=material_map,
        boundaries=None,
        fields=fields,
        contact_map=contact_map,
    )


    thermal_model = ThermalModel(
        bath_temperature=3.0,
        thermal_relaxation_rate=1.0,
        max_substep=1e-10,
)


    result = run_coupled_simulation(
        simulation,
        thermal_model,
        steps=5,
        dt=1e-12,
    )


    assert result.steps == 5


    assert (
        result.simulation.fields.temperature.shape
        ==
        (mesh.ny, mesh.nx)
    )


    assert (
        result.simulation.fields.heat_source.max()
        >
        0
    )
