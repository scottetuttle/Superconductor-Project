import pytest

from shs.config.builder import build_simulation


@pytest.fixture
def nbn_simulation():

    return build_simulation(
        "configs/simulations/nbn_hotspot_test.json"
    )