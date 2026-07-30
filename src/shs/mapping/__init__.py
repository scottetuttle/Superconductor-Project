"""
Mapping package.

Responsible for converting geometry into numerical property maps.
"""

from .region_map import RegionMap, build_region_map
from .material_map import MaterialMap, build_material_map

from .contact_map import (
    ContactMap,
    build_contact_map,
)