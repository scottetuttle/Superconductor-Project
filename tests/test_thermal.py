from shs.physics import ThermalModel


def test_thermal_diffusivity():

    alpha = ThermalModel.thermal_diffusivity(
        thermal_conductivity=10,
        heat_capacity=2e6
    )

    assert alpha == 5e-6