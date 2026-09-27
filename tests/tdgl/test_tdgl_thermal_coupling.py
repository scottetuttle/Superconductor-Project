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

from shs.tdgl.operators import covariant_laplacian

from shs.tdgl import (
    TDGLModel,
    TDGLParameters,
    TDGLBoundarySet,
    TDGLBoundaryCondition,
    TDGLBoundarySide,
    TDGLBoundaryType,
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
    tdgl_boundaries = TDGLBoundarySet()

    for side in TDGLBoundarySide:

        tdgl_boundaries.add(
            TDGLBoundaryCondition(
                side=side,
                type=TDGLBoundaryType.INSULATING,
            )
        )

    simulation = Simulation(
        config=None,
        geometry=geometry,
        mesh=mesh,
        region_map=region_map,
        material_map=material_map,
        boundaries=None,
        tdgl_boundaries=tdgl_boundaries,
        fields=fields,
        contact_map=build_contact_map(
            geometry,
            mesh
    ),
)

    #
    # Start from a uniform superconducting state.
    #

    simulation.fields.psi[:] = (
        1.0 + 0.0j
    )

    #
    # Create a local thermal hotspot.
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
            X - center_x
        )**2
        +
        (
            Y - center_y
        )**2
    ) < radius**2

    hotspot_temperature = (
        0.95 *
        material.Tc
    )

    simulation.fields.temperature[
        hotspot
    ] = hotspot_temperature

    model = TDGLModel(
        TDGLParameters(normalization="legacy_gl", temperature_model="one_minus_t_over_tc")
    )

    #
    # Evolve TDGL.
    #

    for _ in range(1000):

        tdgl_step(
            simulation,
            dt=0.001 * model.characteristic_time(simulation.material_map.materials[0].Tc),
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

    #
    # Expected local equilibrium amplitude.
    #

    hotspot_expected = np.sqrt(
        1.0 -
        hotspot_temperature /
        material.Tc
    )

    outside_temperature = 3.0

    outside_expected = np.sqrt(
        max(
            1.0 -
            outside_temperature /
            material.Tc,
            0.0,
        )
    )

    print(
        "Hotspot |psi|:",
        hotspot_value
    )

    print(
        "Hotspot expected |psi|:",
        hotspot_expected
    )

    print(
        "Outside |psi|:",
        outside_value
    )

    print(
        "Outside expected |psi|:",
        outside_expected
    )


    kinetic = covariant_laplacian(
        simulation.fields.psi,
        simulation.fields.vector_potential_x,
        simulation.fields.vector_potential_y,
        simulation.mesh.dx / material.coherence_length,
        simulation.mesh.dy / material.coherence_length,
    )

    nonlinear = (
        np.abs(simulation.fields.psi)**2 *
        simulation.fields.psi
    )



    print(
        "Hotspot mean |kinetic|:",
        np.mean(np.abs(kinetic[hotspot]))
    )

    

    print(
        "Hotspot mean |nonlinear|:",
        np.mean(np.abs(nonlinear[hotspot]))
    )

    
    #
    # Hotspot must suppress superconductivity.
    #

#
# Hotspot must suppress superconductivity.
#

    assert hotspot_value < outside_value

#
# The hotspot should produce a substantial suppression
# of the order parameter.
#
    suppression_ratio = hotspot_value / outside_value
    print(
        "Hotspot / outside |psi| ratio:",
        suppression_ratio
    )