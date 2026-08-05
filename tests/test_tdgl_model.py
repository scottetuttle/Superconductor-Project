from shs.tdgl import (
    TDGLParameters,
    TDGLModel,
)


def test_tdgl_equilibrium_amplitude():

    model = TDGLModel(
        TDGLParameters()
    )


    assert model.equilibrium_amplitude(
        0.0
    ) == 1.0


    assert model.equilibrium_amplitude(
        1.0
    ) == 0.0


    assert model.equilibrium_amplitude(
        1.2
    ) == 0.0