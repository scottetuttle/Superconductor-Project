from shs.geometry.database import load_geometry
from shs.geometry.mesh import create_mesh

from shs.physics import Fields
from shs.solvers import electrical_step

from shs.mapping import (
    build_material_map,
    build_contact_map,
    build_region_map,
)

from shs.materials.database import get_material


def test_electrical_transport():

    geometry = load_geometry(
        "configs/geometry/NbN_film.json"
    )

    mesh = create_mesh(
        geometry
    )

    fields = Fields.create(
        mesh,
        initial_temperature=3.0,
    )

    material = get_material(
        "NbN"
    )

    material_map = build_material_map(
        build_region_map(
            geometry,
            mesh,
        ),
        material,
    )

    contact_map = build_contact_map(
        geometry,
        mesh,
    )

    updated = electrical_step(
        fields,
        mesh,
        material_map,
        contact_map,
    )

    assert updated.voltage.shape == (
        mesh.ny,
        mesh.nx,
    )

    assert updated.current_density_x.shape == (
        mesh.ny,
        mesh.nx,
    )

    assert updated.current_density_y.shape == (
        mesh.ny,
        mesh.nx,
    )

    assert updated.electric_field_x.shape == (
        mesh.ny,
        mesh.nx,
    )

    assert updated.electric_field_y.shape == (
        mesh.ny,
        mesh.nx,
    )

    assert updated.heat_source.shape == (
        mesh.ny,
        mesh.nx,
    )

    #
    # A voltage gradient should generate Joule heating.
    #

    assert updated.heat_source.max() > 0.0