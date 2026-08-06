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



def test_local_temperature_suppresses_order_parameter():


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



    #
    # Start uniform superconducting state
    #

    simulation.fields.psi[:] = (
        1.0 + 0.0j
    )



    #
    # Create thermal hotspot
    #

    center_x = mesh.nx // 2

    center_y = mesh.ny // 2


    radius = 10


    Y, X = np.ogrid[
        :mesh.ny,
        :mesh.nx
    ]


    hotspot = (

        (
            X-center_x
        )**2

        +

        (
            Y-center_y
        )**2

    ) < radius**2



    simulation.fields.temperature[
        hotspot
    ] = (
        0.95 *
        material.Tc
    )



    model = TDGLModel(
        TDGLParameters()
    )



    #
    # Advance TDGL
    #

    for _ in range(500):

        tdgl_step(
            simulation,
            dt=0.001,
            tdgl_model=model,
        )



    psi_abs = np.abs(
        simulation.fields.psi
    )


    hotspot_value = np.mean(
        psi_abs[hotspot]
    )


    outside_value = np.mean(
        psi_abs[~hotspot]
    )



    print(
        "Hotspot |psi|:",
        hotspot_value
    )

    print(
        "Outside |psi|:",
        outside_value
    )



    #
    # Superconductivity should be
    # suppressed in heated region
    #

    assert hotspot_value < outside_value