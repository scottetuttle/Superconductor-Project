# SHS Configuration Reference

Simulation JSON files are the authoritative place for experiment and solver controls. Packaged fallback values for programmatic callers are in `src/shs/config/defaults.json`. Fundamental SI constants are in `src/shs/utils/constants.py`. Values that are part of an equation or discretization, such as the factor two in a centered difference, remain in the implementation.

## Top-level simulation settings

| Setting | Meaning | Units |
| --- | --- | --- |
| `temperature` | Initial uniform device temperature. | K |
| `current.value` | Requested transport current when `electrical.drive_mode` is `current`. | A |
| `time.duration` | Requested physical simulated duration. | s |
| `time.dt` | Outer physical timestep. Smaller values improve temporal resolution and increase runtime. | s |

## Geometry holes

Geometry JSON may contain a `holes` list for one or more insulating cutouts in the superconducting film. Every hole requires a unique `name`, a `shape` of `circle` or `rectangle`, and center coordinates `x` and `y` in metres. A circle requires `radius`; a rectangle requires `width` and `height`, also in metres. Holes must remain strictly inside the outer film boundary. TDGL, magnetic-source, current-controlled electrical, and thermal operators honor these cutouts. The legacy voltage-driven electrical path retains its older outer-boundary quadrature and is currently qualified for symmetry and qualitative voltage checks rather than precision integrated-current comparisons.

## Thermal settings

| Setting | Meaning |
| --- | --- |
| `bath_temperature` | Temperature toward which bulk relaxation cools the film, in K. |
| `thermal_relaxation_rate` | Bulk cooling rate in s⁻¹. Zero disables bulk bath relaxation. |
| `max_substep` | Largest allowed internal explicit thermal step, in seconds. Smaller is more conservative and slower. |
| `stability_safety_factor` | Fraction of the calculated explicit stability limit allowed. Must be in `(0, 1]`; smaller is safer and slower. |

Thermal boundary entries select `fixed_temperature`, `insulating`, `convection`, or `fixed_heat_flux`. Boundary `temperature` is K, `heat_transfer` is W m⁻² K⁻¹, and positive `heat_flux` is inward W m⁻².

TDGL boundary entries select `insulating`, `normal_contact`, or `normal_contact_mask`. `normal_contact` sets `psi=0` on the complete configured side for compatibility with existing configurations. `normal_contact_mask` sets `psi=0` only where a geometry contact of type `current` intersects that side; uncovered nodes remain gauge-covariantly insulating. Masked TDGL contacts must lie on the outer mesh boundary.

## Moving laser settings

| Setting | Meaning |
| --- | --- |
| `laser.enabled` | Enables time-dependent Gaussian optical heating. |
| `laser.absorbed_power_W` | Full absorbed Gaussian power; power outside the active device or over a hole is not deposited. |
| `laser.sigma_m` | Gaussian standard deviation in metres. Resolve it by several mesh spacings. |
| `laser.repeat_path` | Repeats the path after its final waypoint time. |
| `laser.waypoints` | Ordered `time_s`, `x_m`, and `y_m` objects joined by linear motion. |
| `laser.waypoints[].power_fraction` | Optional fraction of `absorbed_power_W`, linearly interpolated between waypoints. Use it for physical capture ramps and gradual release. Defaults to `1`. |

## Static vortex pinning settings

| Setting | Meaning |
|---|---|
| `pinning.enabled` | Activates the permanent material-defect landscape. |
| `pinning.maximum_suppression` | Caps the summed dimensionless reduction of the linear GL coefficient. |
| `pinning.sites[].x_m`, `y_m` | Physical centre of a defect. |
| `pinning.sites[].sigma_m` | Gaussian defect width; it should span several mesh cells. |
| `pinning.sites[].strength` | Positive dimensionless reduction of the local linear GL coefficient. |

Positive strength reduces local condensation energy, so vortex attraction emerges from TDGL rather than an imposed vortex force. Overlapping sites add before clipping. Width and strength remain phenomenological until calibrated against a measured defect or pinning energy.

The current optical model heats the single effective film temperature. It does not yet resolve optical absorption or separate electron and phonon temperatures.

Electrical drive controls:

| Setting | Meaning |
|---|---|
| `electrical.drive_mode` | `voltage` applies the two configured terminal voltages. `current` applies balanced Neumann current fluxes and uses `current.value`. |
| `electrical.source_contact` | Geometry contact through which positive conventional current enters the film. |
| `electrical.sink_contact` | Geometry contact through which positive conventional current leaves the film. |
| `electrical.reference_voltage` | Mean scalar potential used to fix the arbitrary gauge of a pure-Neumann current solve. It does not change the electric field. |

Direct current drive currently requires `electrical.solver.backend` to be `sparse_direct`. Source and sink contacts must be distinct contiguous segments containing at least two nodes on one outer edge. The finite-volume solve integrates local film thickness, preserves every charge-continuity equation, and rejects incompatible net boundary flux.

Magnetic screening controls:

| Setting | Meaning |
|---|---|
| `electromagnetic.include_self_field` | Enables self-consistent induced vector potential within every coupled physical step. Screening is disabled by default while validation is incomplete. |
| `electromagnetic.screening_tolerance` | Required relative difference between the current-induced vector-potential estimate and the previous screening iterate. |
| `electromagnetic.screening_max_iterations` | Maximum coupled fixed-point iterations allowed when screening is active. |
| `electromagnetic.screening_step_size` | Damping applied to each induced-vector-potential update. Smaller values improve robustness and require more iterations. |
| `electromagnetic.screening_source_stride` | Samples every Nth sheet-current source in both directions. Keep at `1` for quantitative work; larger values are exploratory approximations. |

The implementation uses the thin-film integral `A_ind(r)=mu0/(4*pi) integral K(r')/|r-r'| dA'`, where `K=J*d`. The same-cell singularity on the collocated Cartesian grid is replaced by the analytic average for an equal-area circular cell. This approximation requires mesh-refinement validation.

## Electrical settings

| Setting | Meaning |
| --- | --- |
| `solver.backend` | Electrical linear solver. `sparse_direct` is the accurate fast choice for the current small structured meshes; `red_black_sor` remains available for comparisons. |
| `voltage_left`, `voltage_right` | Contact potentials in volts. Their difference drives the present transport solve. |
| `normal_conductivity_model` | `constant` uses the full normal-state conductivity required by the standard generalized-TDGL current equation. `condensate_depletion` retains the former `(1-|psi|²)sigma_n` interpolation for explicitly legacy comparisons. |
| `solver.tolerance` | Maximum diagonal-scaled electrical equation defect, in volts. Smaller is more accurate and may require more iterations. |
| `solver.max_iterations` | Hard iteration limit for one electrical potential solve. |
| `solver.omega` | Red-black SOR relaxation factor. It must be between 0 and 2; values above 1 may converge faster but can become less robust. |
| `solver.residual_check_interval` | Number of SOR sweeps between equation-defect checks. It has no effect with `sparse_direct`. |

## TDGL settings

| Setting | Meaning |
| --- | --- |
| `normalization` | `pytdgl` uses the mutually matched conductivity-based time, voltage, current, magnetic, and vector-potential scales. `legacy_gl` retains the former SHS scale for comparison only. |
| `temperature_model` | `tc_over_t_minus_one` uses `epsilon=clip(Tc/T-1,-1,1)`. `one_minus_t_over_tc` retains the former coefficient. |
| `u` | Dimensionless order-parameter relaxation parameter. |
| `gamma` | Kramer-Watts-Tobin amplitude/phase coupling parameter used by the generalized TDGL update. Zero recovers the simpler relaxation equation. |
| `kappa` | Ginzburg–Landau parameter. It is stored for electromagnetic coupling and is not yet active in the simplified evolution equation. |
| `max_normalized_timestep` | Largest allowed TDGL step in normalized time. Smaller is more conservative and slower. |
| `stability_safety_factor` | Fraction of the calculated explicit TDGL stability limit allowed. Must be in `(0, 1]`. |
| `include_scalar_potential` | Includes the gauge-covariant scalar-potential term using the solved electrical voltage. Keep it enabled for TDGL transport studies. |

TDGL boundary entries select `insulating` for zero gauge-covariant normal derivative or `normal_contact` for a superconducting-to-normal terminal where `psi = 0`. The latter currently applies to a complete mesh side.

`configs/simulations/nbn_tdgl_small_step.json` is the runnable conservative example. It advances ten 10 fs steps, enables scalar-potential coupling, uses normal contacts on the left and right, and uses the sparse electrical solve.

## Coupling settings

| Setting | Meaning |
| --- | --- |
| `tolerance` | Maximum normalized change across primary coupled fields required for acceptance. |
| `max_iterations` | Maximum fixed-point iterations within one physical timestep. |
| `minimum_iterations` | Minimum iterations before accepting convergence. |
| `check_interval` | Interval for convergence diagnostic checkpoints. |
| `prediction_window` | Recent residual count used to estimate convergence rate. |
| `reasonable_iterations` | Remaining-iteration estimate above which adaptation is recommended. |
| `fast_ratio` | Residual ratio at or below which convergence is classified as fast. |
| `healthy_ratio` | Residual ratio at or below which convergence is classified as healthy. |
| `stall_ratio` | Residual ratio at or above which convergence is classified as stalled. |
| `oscillation_window` | Recent residual count used to detect alternating convergence behavior. |
| `electrical_tolerance_fraction` | Requires the inner electrical defect to be this fraction of the outer voltage convergence target. |
| `field_scales` | Characteristic magnitudes used when a field is near zero. Temperature is K, voltage is V, vector potential is T m, and `psi` is dimensionless. |

## Shared numerical safeguards

| Setting | Meaning |
| --- | --- |
| `max_internal_substeps` | Refuses a single requested step if either explicit integrator would require more internal steps than this. This catches unit mistakes and impractical requests. |
| `contact_roundoff_ulps` | Floating-point tolerance, measured in representable-number spacings, used when mapping geometry endpoints to mesh cells. |

## Scientific constants

`src/shs/utils/constants.py` contains the reduced Planck constant, Boltzmann constant, elementary charge, vacuum permeability, and superconducting flux quantum in SI units. These values describe nature rather than a particular simulation and should not be placed in an experiment JSON.

Material-specific physical constants belong in `configs/materials/*.json`. Every material must explicitly supply `Tc`, coherence length, penetration depth, GL alpha and beta, TDGL `u`, thermal conductivity, volumetric heat capacity, normal resistivity, and thickness. There are no silent zero-valued material fallbacks.

## Reduced optical-manipulation model

Set `physics_backend` to `reduced_vortex` in a diagnostics configuration and run it through `tools/tdgl_diagnostics.py`. The reference file is `tools/config/reduced_optical_tweezers.json`.

The reduced model integrates an overdamped force balance using Bardeen-Stephen viscosity. It includes an attractive temperature-gradient force derived from the temperature dependence of vortex line energy, Gaussian pinning energies, bulk screened vortex-vortex interactions, Lorentz force from a uniform current density, and a first-order thermal lag behind the commanded laser.

`laser.protocol` defines time, position, and peak temperature rise. The peak must remain below `Tc`; moving normal regions require full TDGL. Pinning uses explicit `energy_per_length_J_per_m` and `sigma_m`. The run reports the largest transport lag, final destination error, numerical motion-limiter use, and whether the configured screening criteria were met.

This backend omits vortex-core deformation, nucleation and annihilation, stochastic activation, nonlocal Pearl interactions, and physical edge/image forces. Domain clipping is only a numerical boundary guard. Use the model to identify candidate parameter ranges, then qualify them with full TDGL and calibrated thermal inputs.

## Benchmark-only controls

`tools/config/coupled_convergence_benchmark.json` controls parameter sweeps, adaptive retry factors and bounds, external benchmark heating, output paths, and plotting. These settings affect benchmark experiments only and are intentionally separate from production simulation configuration.

Its `numerical_safeguards` section contains residual, parameter-comparison, and time-comparison floors used only by the experimental adaptive controller. `output.plot_dpi` controls saved-image resolution, while `output.progress_updates` controls approximate console progress frequency. `adaptive_controller.adaptation_effectiveness` sets the minimum residual and iteration improvements required to retain an attempted adaptation.

`electrical.inner_tolerance_policy` selects the initial accuracy of each algebraic electrical solve. `simulation` uses the simulation JSON value exactly. `coupling_scaled` uses `inner_tolerance_factor × coupling tolerance × voltage scale`, clamped by the adaptive-controller bounds. A factor of `0.01` makes the electrical solve roughly two orders of magnitude tighter than the outer coupling target without oversolving it by many additional orders.

`external_heat.amplitude` is peak volumetric heating in W/m³. `external_heat.radius_meters` is the Gaussian standard deviation in meters. Expressing the radius physically keeps the experiment unchanged when mesh resolution changes.

`adaptive_relaxation.method` selects `aitken` for dynamically computed, bounded fixed-point relaxation or `reduction` for the earlier oscillation-triggered multiplicative reduction. `minimum_sigma` and `maximum_sigma` bound Aitken updates; `reduction_factor` and `required_reversals` control the reduction method.
