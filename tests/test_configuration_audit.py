import json
from pathlib import Path
from unittest.mock import patch

from shs.config.builder import build_simulation
from shs.geometry.database import load_geometry
from shs.materials.database import load_material
from shs.utils.output import reserve_output_directory


CONFIG_ROOTS = (Path("configs"), Path("tools/config"), Path("src/shs/config"))


def _leaf_paths(value, prefix=""):
    if isinstance(value, dict):
        for key, item in value.items():
            if key != "_documentation":
                yield from _leaf_paths(item, f"{prefix}.{key}" if prefix else key)
    elif isinstance(value, list) and value and isinstance(value[0], dict):
        keys = sorted(set().union(*(item.keys() for item in value if isinstance(item, dict))))
        for key in keys:
            example = next(item[key] for item in value if key in item)
            yield from _leaf_paths(example, f"{prefix}[].{key}")
    else:
        yield prefix


def test_every_json_configuration_is_valid_and_documents_every_field():
    for root in CONFIG_ROOTS:
        for path in root.rglob("*.json"):
            data = json.loads(path.read_text(encoding="utf-8"))
            documentation = data["_documentation"]
            assert documentation["purpose"]
            assert documentation["relies_on"]
            assert documentation["run_command"]
            assert documentation["authority"]
            missing = set(_leaf_paths(data)) - set(documentation["fields"])
            assert not missing, f"{path} has undocumented fields: {sorted(missing)}"


def test_all_reusable_physics_configs_load():
    for path in Path("configs/materials").glob("*.json"):
        load_material(path)
    for path in Path("configs/geometry").glob("*.json"):
        load_geometry(path)
    for path in Path("configs/simulations").glob("*.json"):
        build_simulation(path)


def test_tdgl_runner_does_not_duplicate_simulation_timestep_or_drive_controls():
    for path in Path("tools/config").glob("*.json"):
        data = json.loads(path.read_text(encoding="utf-8"))
        if "simulation_config" not in data or "run" not in data:
            continue
        if "sample_every" not in data["run"]:
            continue
        assert "dt" not in data["run"]
        assert "steps" not in data["run"]
        assert "transport" not in data


def test_output_directory_reservation_never_reuses_an_existing_run():
    created = set()

    def create_once(path, **_kwargs):
        key = str(path)
        if key in created:
            raise FileExistsError(key)
        created.add(key)

    base = Path("results") / "result"
    with patch.object(Path, "mkdir", autospec=True, side_effect=create_once):
        first = reserve_output_directory(base)
        second = reserve_output_directory(base)
    assert first == base
    assert second == Path("results") / "result_1"
