from shs.materials import get_material

def test_load_nbn():

    material = get_material("NbN")

    assert material.name == "NbN"
    assert material.Tc > 0
    assert material.thickness > 0

nbn = get_material("NbN")

print(nbn)
