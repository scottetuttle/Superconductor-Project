"""Access packaged configuration defaults without importing config builders."""

import json
from copy import deepcopy
from pathlib import Path


_PATH = Path(__file__).parents[1] / "config" / "defaults.json"
with _PATH.open(encoding="utf-8") as _file:
    _DEFAULTS = json.load(_file)


def default_section(name):
    return deepcopy(_DEFAULTS[name])
