"""Static material-pinning landscapes for normalized TDGL."""

import numpy as np


def pinning_suppression(simulation):
    """Return the dimensionless local reduction of the linear GL coefficient.

    Positive configured strengths lower the condensation energy in Gaussian
    defect regions. Overlapping sites add and are clipped at the configured
    maximum so a defect map cannot create an unintended unbounded coefficient.
    """
    shape = simulation.fields.psi.shape
    simulation_config = getattr(simulation, "config", None)
    config = getattr(simulation_config, "pinning", None)
    if config is None:
        return np.zeros(shape, dtype=float)
    if not config.enabled or not config.sites:
        return np.zeros(shape, dtype=float)

    cached = getattr(simulation, "_pinning_suppression_cache", None)
    signature = (
        config.maximum_suppression,
        tuple((site.x_m, site.y_m, site.sigma_m, site.strength)
              for site in config.sites),
        simulation.mesh.nx, simulation.mesh.ny,
        simulation.mesh.dx, simulation.mesh.dy,
    )
    if cached is not None and cached[0] == signature:
        return cached[1]

    x = np.arange(simulation.mesh.nx, dtype=float) * simulation.mesh.dx
    y = np.arange(simulation.mesh.ny, dtype=float) * simulation.mesh.dy
    xx, yy = np.meshgrid(x, y)
    suppression = np.zeros(shape, dtype=float)
    for site in config.sites:
        radius_squared = (xx - site.x_m) ** 2 + (yy - site.y_m) ** 2
        suppression += site.strength * np.exp(
            -radius_squared / (2.0 * site.sigma_m**2)
        )
    suppression = np.minimum(suppression, config.maximum_suppression)
    suppression *= np.asarray(simulation.region_map.active_mask, dtype=bool)
    simulation._pinning_suppression_cache = (signature, suppression)
    return suppression
