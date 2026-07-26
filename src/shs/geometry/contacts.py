from dataclasses import dataclass


@dataclass
class Contact:
    """
    Represents an electrical contact on the device.
    """

    name: str

    contact_type: str

    location: str