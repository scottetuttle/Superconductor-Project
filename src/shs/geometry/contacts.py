from dataclasses import dataclass


@dataclass
class Contact:
    """
    Represents an electrical contact on the device.
    """

    name: str

    contact_type: str

    x: float
    y: float

    x_size: float
    y_size: float
