from shs.physics import ThermalModel

def test_thermal_diffusivity():

    thermal = ThermalModel(
        thermal_conductivity=10,
        heat_capacity=2e6,
        bath_temperature=3,
        thermal_relaxation_rate=0.0
    )

    alpha = thermal.thermal_diffusivity()

    assert alpha == 5e-6