"""
Optical heating models.
"""

from dataclasses import dataclass
import numpy as np

from shs.geometry import Mesh


@dataclass
class GaussianHotspot:

    """
    Gaussian laser heating profile.

    Power density:

    Q(x,y)=A exp(-(r^2)/(2 sigma^2))
    """

    x0: float

    y0: float

    amplitude: float

    sigma: float


    def generate(
        self,
        mesh: Mesh
    ):
        """
        Generate heat source field.
        """

        X, Y = np.meshgrid(
            mesh.x,
            mesh.y,
            indexing="ij"
        )


        r2 = (
            (X-self.x0)**2
            +
            (Y-self.y0)**2
        )


        return (
            self.amplitude *
            np.exp(
                -r2 /
                (2*self.sigma**2)
            )
        )