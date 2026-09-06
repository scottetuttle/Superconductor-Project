"""
Simulation configuration loader.
"""

import json
from dataclasses import dataclass
from pathlib import Path

@dataclass
class SimulationConfig:

    name: str

    geometry: str
    material: str

    temperature: float

    current: float

    duration: float
    dt: float

    boundaries: dict

    


def load_simulation(filepath):

    filepath = Path(filepath)

    with open(filepath, "r") as file:
        data = json.load(file)

    return SimulationConfig(
        name=data["name"],

        geometry=data["geometry"],
        material=data["material"],

        temperature=data["temperature"],

        current=data["current"]["value"],

        duration=data["time"]["duration"],
        dt=data["time"]["dt"],

        boundaries=data["boundaries"]
)