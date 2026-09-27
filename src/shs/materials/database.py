"""
Material database and configuration loader.

Responsible for:
- Loading material definitions from config files
- Creating Material objects
- Providing access to known materials
"""

import json
from pathlib import Path

from .materials import Material



MATERIAL_PATH = Path(
    "configs/materials"
)


def get_material(name: str):
    """
    Retrieve material by name.
    """

    filepath = MATERIAL_PATH / f"{name}.json"

    if not filepath.exists():
        raise FileNotFoundError(
            f"Material {name} not found"
        )

    return load_material(filepath)

def load_material(filepath: str | Path) -> Material:
    """
    Load a material definition from a JSON file.

    Parameters
    ----------
    filepath:
        Path to material configuration file.

    Returns
    -------
    Material
        Initialized material dataclass.
    """

    filepath = Path(filepath)

    with open(filepath, "r") as file:
        data = json.load(file)

    data.pop("_documentation", None)

    return Material(**data)

