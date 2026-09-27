from shs.config.builder import build_simulation


def test_build_simulation():

    simulation = build_simulation(
        "configs/simulations/nbn_hotspot_test.json"
    )


    assert simulation.geometry is not None

    assert simulation.mesh is not None

    assert simulation.region_map is not None

    assert simulation.material_map is not None

    assert simulation.boundaries is not None

    assert simulation.fields is not None


def test_builder_material():

    simulation = build_simulation(
        "configs/simulations/nbn_hotspot_test.json"
    )

    material = (
        simulation.material_map.materials[0]
    )

    assert material.name == "NbN"



def test_builder_boundary():

    simulation = build_simulation(
        "configs/simulations/nbn_hotspot_test.json"
    )

    from shs.boundaries.boundary import BoundarySide

    left = simulation.boundaries.get(
        BoundarySide.LEFT
    )

    assert left.temperature == simulation.config.boundaries["left"]["temperature"]

def test_builder_fields():

    simulation = build_simulation(
        "configs/simulations/nbn_hotspot_test.json"
    )

    assert simulation.fields.temperature[0, 0] == simulation.config.temperature

    assert simulation.fields.temperature.shape == (
        simulation.mesh.ny,
        simulation.mesh.nx,
    )
