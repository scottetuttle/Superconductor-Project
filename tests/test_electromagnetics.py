from shs.physics import ElectromagneticModel


def test_default_electromagnetic_model():

    model = ElectromagneticModel()


    assert model.reference_voltage == 1.0

    assert model.ground_voltage == 0.0

    assert model.include_self_field is False