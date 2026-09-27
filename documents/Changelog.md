# SHS Changelog

## 2026-09-16 - Configuration propagation audit

- Added end-to-end tests proving that simulation JSON controls geometry, material, initial temperature/current, thermal and TDGL model construction, laser deposition, and pinning landscapes.
- Added `configuration_provenance` and `effective_simulation_config` to diagnostic summaries.
- Added an immediate console warning when diagnostics step count and timestep replace a different simulation duration.
- Documented intentional precedence and configuration fields that are currently inactive.
- Full regression result: 261 passed.

## 2026-09-16 - Visual experiment configuration editor

- Added a local browser editor for rectangular devices with circular and rectangular holes.
- Added placement and editing of signed initial vortices and timed laser waypoints with power fractions.
- Added controls for device dimensions, mesh, material, temperature, current, magnetic field, laser size and power, timestep, output sampling, screening, scalar potential, and terminal boundary type.
- Added geometry, simulation, diagnostics, and reusable editor-project JSON exports plus structural validation warnings.

## 2026-09-16 - TDGL optical-protocol bridge

- Added optional per-waypoint `power_fraction` to the full TDGL laser configuration.
- Position and absorbed power now interpolate together, enabling capture ramps, constant-power transport, destination holds, and gradual release without applying an artificial vortex force.
- Preserved existing laser configurations with a default power fraction of one.
- Added interpolation, deposited-power, and validation tests. The full suite passes 254 tests.

## 2026-09-16 - Reduced optical vortex manipulation backend

- Added an optional overdamped point-vortex backend with Bardeen-Stephen drag, thermal attraction, explicit pinning energies, vortex interactions, Lorentz force, and thermal response lag.
- Added a selectable diagnostics path through `tools/tdgl_diagnostics.py`.
- Added a capture, 1 m/s transport, hold, and release configuration using illustrative niobium-scale parameters.
- Added trajectory, force, lag, destination, and numerical-limiter diagnostics.
- The reference screening run transports one vortex 14 micrometres, stays within 103 nm of the moving thermal centre, and remains at the destination after release without activating the motion limiter.

## 2026-09-16 - Static vortex pinning and capture diagnostics

- Added configurable Gaussian material defects that reduce the local linear TDGL coefficient and create condensation-energy pinning wells.
- Added additive overlap with a configurable suppression cap and cached the static landscape.
- Added nearest laser-vortex separation to each diagnostic frame and `laser_vortex_distance.png` to laser runs.
- Added summary fields that report separation while withholding a capture claim until sustained following and release criteria exist.
- Added tests for map geometry, clipping, inactive regions, caching, and local condensate suppression. The full suite passes 248 tests.

## 2026-09-16 - Moving laser heating and vortex trajectories

- Added absorbed Gaussian optical power following piecewise-linear physical waypoints.
- Coupled laser deposition into transactional electrothermal TDGL steps at timestep midpoints.
- Ensured light over holes or outside the device is not reassigned to active material.
- Cached optical coordinate grids and repeated coupling-iteration heat profiles.
- Added signed vortex identity tracking, CSV trajectory output, and a trajectory/laser-path plot.
- Added an editable moving-laser reference experiment and physical validation tests.

## 2026-09-15 - Fluxoid, gauge, Meissner, and screening qualification

- Added rectangular-contour physical London fluxoid and topological phase-winding diagnostics.
- Corrected a factor-of-four pyTDGL supercurrent conversion error exposed by the independent fluxoid comparison.
- Added discrete gauge-equivalence checks for supercurrent, fluxoid, induced vector potential, and total magnetic field.
- Added a Meissner-response test confirming that induced field opposes a positive applied field.
- Added kernel mesh-refinement checks and per-frame fluxoid diagnostics to the TDGL diagnostic tool.
- Three real screened 0.01 T, 10 nA gauge cases produced identical observables and a seeded-vortex London fluxoid of `0.999809 Phi0` versus topological `1 Phi0`.
- A zero-vortex screened case reduced the center field from `0.010000000 T` to `0.009999647 T`, consistent with weak screening for this film's large Pearl length.
- Validation: **229 tests passed in 99.50 seconds**.

## 2026-09-15 - Self-consistent magnetic-screening foundation

- Separated applied, induced, and total vector potentials in the evolving field state.
- Added an FFT thin-film Biot-Savart vector-potential operator using sheet current and an equal-area same-cell regularization.
- Added damped induced-field updates to the transactional coupled solver and required the magnetic residual to pass before accepting a screened step.
- Added screening tolerance, iteration limit, damping, and source-stride configuration.
- Kept screening disabled by default pending fluxoid, gauge, refinement, and external-reference validation.
- Added linearity, component direction, field composition, configuration, and coupled-convergence tests.
- A 100x100, 10 nA, 0.01 T screened smoke step converged in 21 iterations to a relative screening residual of `7.97e-4` with a requested tolerance of `1e-3`.
- Validation: **224 tests passed in 102.93 seconds**.



## 2026-09-15 - Direct-current diagnostics and finite-volume continuity

- Switched TDGL diagnostic `current` mode from repeated voltage searches to the production Neumann current-flux solver.
- Retained the former controller as explicit `voltage_search` mode for comparisons.
- Added reversed-current and bottom-to-top terminal validation alongside existing full and partial left-to-right cases.
- Replaced the boundary-inaccurate nodal continuity diagnostic with face-flux balance over physical finite-volume cells, including thickness and half-width boundary volumes.
- Validated two-step 10 nA coupled diagnostics at 0 T and 1 T with no voltage trials, relative current errors below `6e-13`, and relative continuity defects below `8e-16`.
- Validation: **221 tests passed in 98.71 seconds**.

## 2026-09-15 - Direct current-flux electrical drive

- Added `electrical.drive_mode`, source/sink contact selection, and a reference-voltage gauge setting.
- Added finite-volume Neumann current injection for contiguous partial or complete outer-edge contacts.
- Added an augmented sparse solve that constrains mean voltage without deleting a charge-continuity equation.
- Included local film thickness and edge control-volume weights in terminal-current integration.
- Retained voltage drive as the default for existing configurations and diagnostics.
- Added analytic full-edge and partial-contact current-conservation tests and a coupled 100x100 transport smoke check.
- Validation: **218 tests passed in 100.56 seconds**.

## 2026-09-15 - Physics-engine foundation and partial terminals

- Added `documents/PhysicsEngineImplementationContext.md` as a durable human- and AI-readable brief for the Josephson-junction, hotspot, electron-phonon, SQUID, magnetic-screening, circuit, and model-selection program.
- Added the explicit TDGL boundary type `normal_contact_mask`, which applies `psi=0` only where a geometry current-contact mask intersects the configured edge.
- Made uncovered portions of a masked-contact edge retain the gauge-covariant insulating condition.
- Preserved the existing `normal_contact` whole-side behavior so older configurations remain reproducible.
- Added validation for partial masks, invalid interior masks, solver-level partial terminals, and legacy full-side terminals.
- Validation: **216 tests passed in 99.80 seconds**.

## 2026-09-15 - Phase GIF visibility correction

- Expanded `order_parameter.gif` to show amplitude, absolute phase, and spatial phase change from the initial state.
- Removed the gauge-dependent uniform phase rotation from the phase-change view, masked phase where the order parameter is effectively zero, and dynamically scaled the remaining change so small transport-driven phase gradients are visible.
- Added the maximum spatial phase change to each frame title and documented why the fixed-range absolute-phase panel can appear stationary.
- Validation: **213 tests passed in 103.71 seconds**. A real 10-step, 10 nA diagnostic produced 11 distinct rendered phase-change panels.

## 2026-09-15 - TDGL normalization and normal transport

- Added an explicit `pytdgl` normalization with matched conductivity-based time, voltage, current-density, magnetic-field, and vector-potential scales.
- Made `epsilon=clip(Tc/T-1,-1,1)` the default temperature coefficient and retained the former equation under explicit legacy settings.
- Corrected the standard normal-current model to use full normal-state conductivity; retained condensate depletion only as a named legacy option.
- Added a near-`Tc`, 2 nm NbN reference material and updated the transport diagnostic configuration to satisfy its declared 2D TDGL regime checks.
- Added machine-readable temperature, thickness, and mesh-validity diagnostics and documented a deterministic multi-model architecture with a limited future role for machine learning.
- Validation: **212 tests passed in 99.33 seconds**.

## 2026-09-15 - Current and magnetic diagnostic expansion

- Added full-edge NbN transport geometry and a current-controlled diagnostic mode that retries voltage from the same accepted physical state.
- Added uniform perpendicular applied fields in symmetric and Landau gauges, numerical curl validation, and gauge transformation of the initial order parameter.
- Added FFT-accelerated thin-film Biot-Savart field diagnostics, a stored perpendicular magnetic field, cross-section current checks, magnetic flux/energy metrics, current-component maps, current streamlines, and magnetic-field plots.
- Added `documents/CurrentMagneticAudit_2026-09-15.md` with the equation audit, research basis, measured behavior, and remaining self-consistent-screening work.
- Validation: **207 tests passed in 59.63 seconds**. The supplied 0 T/1 T, 10 nA, 100-step comparison converged in both cases.

## 2026-09-15 - TDGL configuration diagnostics

- Added `tools/tdgl_diagnostics.py` for independent Cartesian configuration sweeps using the production coupled solver.
- Added signed gauge-invariant vortex detection, initial vortex and hotspot probes, diagnostic CSV/JSON histories, comparison plots, field snapshots, kymographs, and optional GIFs.
- Added a fully explained configuration at `tools/config/tdgl_diagnostics.json` and usage guidance in `documents/TDGLDiagnostics.md`.
- Corrected the 3 K small-step example's inherited fixed thermal boundaries from 9 K to 3 K after the diagnostic plots exposed the mismatch.
- Validation: **201 tests passed in 58.11 seconds**; the supplied two-case diagnostic sweep completed with both cases converged and a seeded +1 vortex tracked through all ten frames.

## 2026-09-15 - Small-step TDGL physics restart

- Added gauge-covariant scalar-potential evolution and activated the generalized Kramer-Watts-Tobin `gamma` update.
- Added full-side normal-contact TDGL boundaries and retained link-variable insulating boundaries.
- Added an accurate sparse-direct electrical backend and current-continuity, Joule-power, thermal-energy, temperature, and order-parameter diagnostics to the adaptive benchmark.
- Added `configs/simulations/nbn_tdgl_small_step.json`, a runnable ten-step, 10 fs TDGL example.
- Documented the remaining terminal, magnetic, calibration, thermal, and time-integration work in `TDGLPhysicsRoadmap_2026-09-15.md`.

**Validation:** **198 tests passed in 60.20 seconds**. The complete small-step example also converged for all ten steps with finite fields.

## 2026-09-14 — Adaptive benchmark audit

- Corrected full-state rollback, accepted-time-level integration, derived-field exchange, ordering consistency, and configuration propagation in the adaptive coupled benchmark.
- Added physical RMS residual scaling, coupling-scaled electrical accuracy, a meter-based centered heat source, safeguarded Aitken relaxation, solver-work accounting, and focused regression tests.
- Added `documents/AdaptiveBenchmarkAudit_2026-09-14.md` with verified behavior, remaining physics limitations, performance evidence, related-project research, and the recommended adaptive-controller architecture.
- Validation: **189 tests passed in 113.79 seconds**.

## 2026-09-14 — Configuration and scientific-constant audit

- Centralized fundamental SI values in `src/shs/utils/constants.py`.
- Added packaged JSON defaults for programmatic callers and removed duplicated solver defaults from production modules.
- Exposed thermal and TDGL stability factors, coupling classifications, inner-tolerance coordination, internal-substep limits, and contact roundoff control.
- Made all material properties explicit in material JSON instead of relying on hidden dataclass values.
- Added `documents/ConfigurationReference.md`, describing every simulation control, its units, and its effect.
- Preserved equation coefficients and indexing values beside their numerical methods because they define those methods rather than tune a run.
- Moved adaptive-benchmark fallback thresholds, progress frequency, and output resolution into the benchmark JSON reference.

**Validation:** `181 passed in 280.80 seconds` with the current 50 nm NbN geometry and configuration.

## 2026-09-14 — Solver corrections and validation

- Fixed coupling iterations to advance one physical timestep from the accepted state; failed steps and exceptions preserve that state.
- Applied SI conversions at TDGL interfaces and unified superconducting electrical coupling across solver entry points.
- Corrected electrical source signs, insulating edges, equation-residual checks, and coordination of inner and outer tolerances.
- Implemented conservative thermal diffusion, configured thermal boundaries, and mesh-dependent explicit timestep safeguards.
- Separated external and Joule heating; updated benchmark heating adapters to preserve external sources.
- Corrected TDGL boundary links, rectangular hotspot orientation, and contact-endpoint roundoff.
- Exposed coupling controls, declared dependencies, and updated usage and status documentation.
- Corrected stale configuration expectations, timestep units, and test assumptions; added 22 regression cases.

**Validation:** Original bounded baseline: 137 passed, 8 failed, with 14 tests excluded due to excessive timestep workloads. Final full suite: **181 passed in 53.01 seconds**, with no skips. The NbN smoke test completed two physical steps successfully.

Existing benchmark outputs were preserved. Full electromagnetic, coupled power-balance, and experimental validation remain outstanding.

See the [detailed correction report](SolverCorrections_2026-09-14.md) and [final test log](validation/2026-09-14/final_tests.log).
# 2026-09-15: Multiply connected TDGL validation

- Added circular and rectangular insulating holes, an active superconducting-domain mask, internal covariant Neumann boundaries, and named hole fluxoid contours.
- Added an NbN ring reference configuration and seven analytic/physical qualification tests.
- Confirmed the complete suite passes with 236 tests.
- Documented that normal-current and thermal hole handling remains the next required step before perforated electrothermal devices are quantitative.
# 2026-09-15: Perforated electrical and thermal transport

- Removed electrical and thermal face connections across inactive geometry.
- Reduced pure-Neumann current matrices to active nodes and rejected contacts that overlap holes.
- Prevented inactive cells from carrying current, producing Joule heat, diffusing heat, or coupling to bath relaxation.
- Added perforated-domain tests for symmetric flow, requested-current recovery, zero leakage, and thermal energy conservation.
- Added drawn polygonal devices with any number of holes as an explicit future geometry-system goal.
- Confirmed the complete suite passes with 240 tests.
# 2026-09-15: Human results dashboard

- Added `tools/build_results_dashboard.py` and generated `benchmark_results/index.html`.
- The dashboard indexes all diagnostic runs, metrics, validity warnings, plots, animations, and raw CSV/JSON files.
- Added a SQUID-readiness panel that separates validated ring capabilities from missing Josephson-junction requirements.
# 2026-09-15: Editable ring fluxoid exploration

- Added a screened NbN ring experiment with current bias, applied field, trapped winding, and an off-center hotspot.
- Generated plots, kymographs, diagnostics, and an order-parameter GIF for inspection through the results dashboard.
- Corrected diagnostic summaries and plot titles to distinguish self-consistent screening from diagnostic-only self-fields.
# 2026-09-16: Repository cleanup and long-duration TDGL plan

- Removed 1,970 superseded generated files and unused placeholders, reclaiming 229.15 MiB while preserving current user experiments and representative validation evidence.
- Removed two obsolete convergence benchmark variants and kept the maintained baseline and adaptive benchmark tools.
- Made generated benchmark-result directories ignored and retained `pyproject.toml` as the dependency authority.
- Added `CodebaseCleanup_2026-09-16.md` and `LongDurationTDGLImplementationPlan_2026-09-16.md`.
- Defined the next performance milestone as a validated sparse gauge-covariant IMEX TDGL backend with adaptive stepping, multirate coupling, and checkpoint/restart support.
- Validation: bytecode compilation succeeded and **259 tests passed**; the only warning was pytest being unable to write its ignored cache directory.
# 2026-09-16: Current-physics demonstration set

- Added four documented diagnostic configurations for vortex-antivortex relaxation, current-reversal vortex motion, fluxoid response to field bias, and laser capture/drag/release with a current comparison.
- Focused interpretation on validated qualitative behavior and clearly identified seeded vortices and uncalibrated absolute rates.
- Added `CurrentPhysicsDemonstrations_2026-09-16.md` with run commands, expected observables, and interpretation limits.
- End-to-end laser comparison validation completed both cases and reproduced the requested 5 nA terminal current with a relative current-continuity residual near `1.2e-14`.
# 2026-09-17: Human-facing physics and computation documentation

- Added a human documentation entry point with a capability-confidence table and result-quality checklist.
- Added subsystem-by-subsystem explanations of TDGL, normal current, heat, magnetic screening, vortices, holes, fluxoids, pinning, laser heating, and reduced vortex dynamics.
- Added the governing equations, normalization, discretizations, solver strategies, coupling order, rollback behavior, performance characteristics, and source map.
- Clearly separated implemented behavior, effective approximations, diagnostic models, and missing physics.
# 2026-09-17: Spatial TDGL Josephson weak link and researcher plots

- Added a configured vertical weak link that locally suppresses the TDGL linear coefficient.
- Added gauge-invariant junction phase, terminal voltage, weak-link amplitude, and bank-amplitude diagnostics.
- Added current-sweep tooling with I–V, differential resistance, phase response, condensate suppression, and transient plots.
- Added an NbN example and a five-point ±10 nA sweep; every case converged and showed the expected odd voltage/phase symmetry.
- Documented that the short sweep is nearly linear and does not yet establish a critical current or tunnel-junction behavior.
- Added explicit weak-link overlays, cross-device junction profiles, and integrated normal/supercurrent diagnostics. The +10 nA example carries about 9.7% supercurrent and 90.3% normal current after 200 fs.
- Extended transient plots with the normal/supercurrent partition; a 1 ps trial increased the supercurrent share to about 26.4%, confirming that phase-relaxation time dominates the short-run result.
- Added a stateful zero-current, positive-bias, current-inversion protocol for multiple current magnitudes with adaptive stabilization, per-step CSV data, reversal plots, and optional decimated GIF capture.
- Reused the constant current-drive sparse factorization and introduced a dedicated 100 nm × 50 nm junction mesh at `0.1xi` spacing. A ten-step benchmark improved from 25.4 s with repeated factorization to 14.2 s with reuse on the former mesh.
# 2026-09-18: Lumped SIS Josephson dynamics

- Added a selectable classical SIS/RCSJ physics backend with current and voltage bias, resistive and capacitive branches, DC/AC drive, Josephson voltage-phase evolution, and magnetic Fraunhofer modulation.
- Added optional Ambegaokar–Baratoff critical-current calculation from temperature, critical temperature, and tunnel resistance.
- Added component-resolved CSV output and phase, voltage, Josephson-current, and I–V plots.
- Validated the supplied `Ic = 2 µA` example: ±1.5 µA remains at zero average voltage, while ±2.5 and ±4 µA produce symmetric running-phase voltage states.
- Kept electron–phonon two-temperature dynamics as the next thermal addition; exported junction Joule power provides its future coupling input.
- Added finite-width junction interference from applied field and stationary electrode vortices, including effective critical-current suppression and phase offsets.
- Added microwave current-drive sweeps with Shapiro voltage normalization, locking errors, aperture-current textures, phase histories, and vortex/no-vortex comparisons.
- Validated 5 GHz integer locking at the expected 10.339 microvolt spacing for no-vortex, single-vortex, and opposite-vortex-pair scenarios.
- Expanded the study to a 3-by-4 sweep of vortex layouts and microwave amplitudes, including a zero-microwave control, with one normalized Shapiro matrix and separate CSV files for all twelve combinations.
- Added idealized electrode phase maps and optional phase/aperture-current GIFs at a selected microwave and DC bias.
- Added a 40 ns focused Shapiro run with 0.1 microamp bias spacing, voltage variance, phase-rate locking error, differential resistance, and a final-two-cycle phase-aware GIF.
- Added frequency-aware ODE stepping through `microwave_steps_per_cycle` so high-frequency drives cannot silently exceed the requested temporal resolution.
# 2026-09-18 — Configuration and runner audit

- Added complete inline documentation to every JSON configuration and packaged
  default, including purpose, dependencies, run command, authority, and every
  leaf field.
- Removed duplicate TDGL runner timestep/current controls; simulation configs
  now own physical time and sweeps own case-specific current values.
- Added atomic numbered output allocation to every result-producing runner so
  repeated runs preserve earlier results.
- Fixed reduced-vortex wrapper output reporting and removed two empty,
  unreferenced material files.
- Added configuration-loading, documentation-coverage, precedence, and output
  reservation regression tests.
# 2026-09-21 — Electron–phonon and TDGL time baseline

- Added an optional two-temperature film model with separate electron and
  phonon energy balances, conservative local exchange, phonon escape, and
  electron-only Joule/laser deposition.
- Added a gauge-compatible TDGL Heun predictor–corrector as the default while
  retaining Euler as a selectable comparison method.
- Added separate thermal-energy and temperature diagnostics, a runnable
  illustrative hotspot configuration, and conservation, temporal-order, and
  coupled-step regression tests.
- Documented the uncalibrated material parameters, whole-solver accuracy
  limits, and fluctuation–dissipation requirements before stochastic forcing.
