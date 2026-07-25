"""
Defines electrical contacts on device geometries.
"""

from dataclasses import dataclass


@dataclass
class Contact:
    """
    Electrical connection to a device.

    Parameters
    ----------
    name:
        Contact identifier.

    contact_type:
        Type of contact.

        Examples:
        - current
        - voltage

    location:
        Position description.
    """

    name: str
    contact_type: str
    location: str