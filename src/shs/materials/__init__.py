"""
SHS Material System

Provides material definitions and database access.
"""

from .materials import Material
from .database import get_material

__all__ = [
    "Material",
    "get_material",
]