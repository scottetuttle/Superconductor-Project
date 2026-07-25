"""
Defines spatial regions inside a device geometry.
"""

from dataclasses import dataclass


@dataclass
class Region:
    """
    A named region of a device.

    Examples:
    - superconducting film
    - contact
    - hotspot
    """

    name: str
    region_type: str
    material: str | None = None