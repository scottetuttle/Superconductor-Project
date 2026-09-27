"""Moving Gaussian optical heating for thin superconducting films."""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class LaserPosition:
    x_m: float
    y_m: float


def laser_power_fraction(config, time_s):
    """Interpolate the dimensionless absorbed-power envelope."""
    if not np.isfinite(time_s) or time_s < 0:
        raise ValueError("Laser evaluation time must be finite and nonnegative.")
    points = config.waypoints
    if not points:
        raise ValueError("Laser path requires at least one waypoint.")
    if config.repeat_path and len(points) > 1:
        duration = points[-1].time_s - points[0].time_s
        if duration <= 0:
            raise ValueError("Repeating laser path requires positive duration.")
        time_s = points[0].time_s + ((time_s - points[0].time_s) % duration)
    if time_s <= points[0].time_s:
        return points[0].power_fraction
    if time_s >= points[-1].time_s:
        return points[-1].power_fraction
    times = np.array([point.time_s for point in points])
    upper = int(np.searchsorted(times, time_s, side="right"))
    first, second = points[upper - 1], points[upper]
    fraction = (time_s - first.time_s) / (second.time_s - first.time_s)
    return first.power_fraction + fraction * (
        second.power_fraction - first.power_fraction
    )


def laser_position(config, time_s):
    """Interpolate a configured piecewise-linear beam path in physical time."""
    if not np.isfinite(time_s) or time_s < 0:
        raise ValueError("Laser evaluation time must be finite and nonnegative.")
    points = config.waypoints
    if not points:
        raise ValueError("Laser path requires at least one waypoint.")
    if config.repeat_path and len(points) > 1:
        duration = points[-1].time_s - points[0].time_s
        if duration <= 0:
            raise ValueError("Repeating laser path requires positive duration.")
        time_s = points[0].time_s + ((time_s - points[0].time_s) % duration)
    if time_s <= points[0].time_s:
        return LaserPosition(points[0].x_m, points[0].y_m)
    if time_s >= points[-1].time_s:
        return LaserPosition(points[-1].x_m, points[-1].y_m)
    times = np.array([point.time_s for point in points])
    upper = int(np.searchsorted(times, time_s, side="right"))
    first, second = points[upper - 1], points[upper]
    fraction = (time_s - first.time_s) / (second.time_s - first.time_s)
    return LaserPosition(
        first.x_m + fraction * (second.x_m - first.x_m),
        first.y_m + fraction * (second.y_m - first.y_m),
    )


def gaussian_laser_heat_source(simulation, time_s):
    """Return volumetric optical heating in W/m^3 at ``time_s``.

    ``absorbed_power_W`` is the total power of the full Gaussian on an
    unbounded film. Power falling outside the active device or into a hole is
    not renormalized back into the film.
    """
    shape = simulation.fields.temperature.shape
    simulation_config = getattr(simulation, "config", None)
    config = getattr(simulation_config, "laser", None)
    if config is None:
        return np.zeros(shape), None
    if not config.enabled or config.absorbed_power_W == 0:
        return np.zeros(shape), None
    position = laser_position(config, time_s)
    coordinate_grid = getattr(simulation.mesh, "_laser_coordinate_grid", None)
    if coordinate_grid is None:
        coordinate_grid = np.meshgrid(simulation.mesh.x, simulation.mesh.y)
        simulation.mesh._laser_coordinate_grid = coordinate_grid
    xx, yy = coordinate_grid
    profile = np.exp(-((xx - position.x_m) ** 2 + (yy - position.y_m) ** 2)
                     / (2.0 * config.sigma_m**2))
    instantaneous_power = config.absorbed_power_W * laser_power_fraction(config, time_s)
    areal_power = instantaneous_power * profile / (2.0 * np.pi * config.sigma_m**2)
    heat = np.zeros(shape)
    active = simulation.region_map.active_mask
    heat[active] = areal_power[active] / simulation.material_map.thickness[active]
    return heat, position


def apply_laser_heat_source(simulation, time_s):
    config = getattr(getattr(simulation, "config", None), "laser", None)
    key = None if config is None else (
        float(time_s), config.enabled, config.absorbed_power_W, config.sigma_m,
        tuple((p.time_s, p.x_m, p.y_m, p.power_fraction) for p in config.waypoints),
    )
    cached = getattr(simulation, "_laser_heat_cache", None)
    if cached is not None and cached[0] == key:
        heat, position = cached[1], cached[2]
    else:
        heat, position = gaussian_laser_heat_source(simulation, time_s)
        simulation._laser_heat_cache = (key, heat, position)
    simulation.fields.laser_heat_source = heat
    simulation.fields.laser_position_x_m = None if position is None else position.x_m
    simulation.fields.laser_position_y_m = None if position is None else position.y_m
    return heat
