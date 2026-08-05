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



def test_tdgl_low_temperature_stability():


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
        contact_map=build_contact_map(
            geometry,
            mesh
        ),
    )


    model = TDGLModel(
        TDGLParameters()
    )


    #
    # Start superconducting state
    #

    simulation.fields.psi[:] = 1.0 + 0.0j


    initial = np.abs(
        simulation.fields.psi.copy()
    )



    tdgl_step(
        simulation,
        dt=0.01,
        tdgl_model=model,
    )



    final = np.abs(
        simulation.fields.psi
    )


    assert np.allclose(
        final,
        initial,
        atol=1e-3
    )