# Physics engine implementation context

## Purpose

This document is the durable implementation brief for humans and future AI agents working on SHS. The primary objective is physical accuracy for Josephson junctions, superconducting hotspots, electron-phonon nonequilibrium, SQUIDs, vortices, and related thin-film superconducting technologies. Runtime matters, but an inexpensive result must never be selected when its model is outside its declared validity range.

## Current trustworthy scope

SHS currently has a normalized generalized TDGL order-parameter solver, scalar-potential electrical continuity, conservative single-temperature heat diffusion, electrothermal fixed-point coupling, prescribed applied vector potential, and diagnostic-only thin-film Biot-Savart fields. The `pytdgl` normalization and constant normal-state conductivity model are the preferred path. The legacy normalization and condensate-depletion conductivity exist only for reproducibility.

The current model is suitable for controlled, small-step, two-dimensional dirty-film studies near `Tc`. It is not yet a quantitatively complete model for screened magnetic devices, real current terminals, low-temperature junctions, or nonequilibrium detector hotspots.

## Non-negotiable physical principles

1. Preserve gauge covariance in every spatial and temporal TDGL operation.
2. Measure charge continuity, terminal-current balance, energy balance, and fluxoid consistency; solver convergence alone is insufficient.
3. Keep normalization scales internally matched. Never mix `legacy_gl` and `pytdgl` scales.
4. Every physics profile must state its equations, assumptions, material inputs, validity tests, numerical methods, and required acceptance diagnostics.
5. Validate new physics against analytic limits, convergence studies, and an independent published implementation before claiming quantitative accuracy.
6. Model selection must be deterministic and validity constrained. Machine learning may later predict cost or flag anomalies, but it must not override physical validity rules.

## Target model profiles

- `normal_electrothermal`: normal current plus heat flow for inexpensive normal-state studies.
- `tdgl_applied_field`: generalized TDGL with a prescribed field when screening is demonstrably negligible.
- `tdgl_screened`: TDGL coupled self-consistently to induced vector potential for vortices and magnetic devices.
- `junction_rcsj`: lumped resistively and capacitively shunted junction dynamics.
- `junction_tdgl`: spatially resolved weak links or interface Josephson conditions.
- `tdgl_two_temperature`: TDGL coupled to electron and phonon temperatures.
- `tdgl_screened_two_temperature`: full thin-film magnetic and nonequilibrium electrothermal model.

## Ordered implementation program

### Stage 1: physically consistent terminals

Use geometry contact masks for partial-edge normal TDGL terminals. `normal_contact_mask` applies the covariant insulating condition to uncovered edge nodes and `psi=0` only on current-contact nodes; legacy `normal_contact` deliberately keeps its whole-side meaning. Add current-flux boundary terms to the scalar-potential equation, fix its pure-Neumann gauge with one reference potential, and verify equal signed source and sink currents. Current-controlled production simulations must no longer require repeated voltage trials.

Status: partial-edge TDGL terminal masks and direct scalar-potential current flux are implemented. Current mode uses a finite-volume pure-Neumann solve with a mean-voltage gauge constraint; full-edge and matching partial-edge contacts reproduce requested cross-section current. Reversed current and horizontal or vertical terminal orientations are tested. The TDGL diagnostics use this production current mode directly and retain voltage search only as an explicit legacy comparison. Coupled 0 T and 1 T checks reach 10 nA with finite-volume continuity defects near machine precision. Heterogeneous-thickness terminal studies and more complex multi-terminal rules remain future terminal work; magnetic screening is now the next primary stage.

### Stage 2: magnetic screening

Separate applied and induced vector potentials. Compute sheet current, solve the thin-film induced vector potential, and iterate current and field to a configured tolerance within each physical step. Add closed-contour flux and fluxoid diagnostics. Validate gauge invariance, the Meissner state, flux quantization, equilibrium vortices, and the negligible-screening limit.

Status: the first screening milestone and its initial qualification layer are implemented. Fields store applied, induced, and total vector potentials separately. A regularized FFT thin-film Biot-Savart convolution computes induced vector potential from sheet current, and the transactional coupled solver performs damped fixed-point updates until both its ordinary field residual and the relative induced-vector-potential residual pass. Physical and topological fluxoids, screened gauge equivalence, weak Meissner response, vortex quantization, and kernel refinement now have automated checks. These checks exposed and corrected a factor-of-four pyTDGL supercurrent conversion error. Screening remains disabled by default. Remaining work includes a strong-screening benchmark geometry, equilibrium-vortex refinement, hole/SQUID contours, heavy-ball acceleration, and quantitative trajectory comparison with pyTDGL.

Multiply connected geometry is now implemented for named circular and rectangular insulating holes. Active-domain gauge links prevent supercurrent leakage, electrical and thermal face conductances vanish across internal boundaries, pure-Neumann current solves use only active nodes, and named hole contours support fluxoid measurement. Tests cover geometry, isolation, gauge covariance, flux quantization, applied-field flux, requested-current recovery, symmetric electrical and thermal flow around a hole, and thermal energy conservation. The current-controlled finite-volume path is the quantitative transport reference; the legacy voltage-driven path still uses its older outer-boundary weighting. See `MultiplyConnectedValidation_2026-09-15.md`.

### Future geometry system

The geometry layer must eventually accept a drawn device outline with zero, one, or many cutouts instead of requiring a rectangular film plus analytic holes. The target pipeline should import or construct polygons, rasterize one authoritative active-domain mask at a declared resolution, preserve named boundaries and terminals, and detect disconnected or under-resolved features before solving. SVG or a simple interactive polygon editor is a suitable human-facing format, while the saved simulation input should contain deterministic physical coordinates and geometry metadata. All physics solvers should continue consuming the common active mask and face connectivity so the number of holes does not create solver-specific geometry logic.

### Stage 3: temporal accuracy and performance

Replace explicit TDGL diffusion with an IMEX or semi-implicit covariant solve and add an embedded temporal-error estimate. Add an implicit thermal method. Cache electrical matrix topology and factorizations. Retain rollback semantics. Demonstrate mesh, timestep, and solver-tolerance convergence before optimizing parallel execution.

### Stage 4: Josephson junctions

Support spatial weak links through local `Tc`, TDGL coefficients, conductivity, and thickness. Add gauge-invariant interface junction conditions for barriers too thin to resolve. Add an RCSJ profile and circuit coupling for inexpensive lumped simulations. Validate the DC and AC Josephson relations, current-phase curves, Shapiro steps, Fraunhofer patterns, and overdamped and hysteretic IV limits.

### Stage 5: hotspot and electron-phonon physics

Retain the single-temperature model as a fast profile. Add electron and phonon energy equations with temperature-dependent heat capacities and conductivities, electron-phonon transfer, phonon escape to the substrate, optical deposition profiles, and order-parameter relaxation energy. Add quasiparticle diffusion only when two-temperature physics is insufficient. Validate deposited-energy conservation, hotspot growth, recovery time, and retrapping behavior.

2026-09-21 checkpoint: an optional two-temperature branch now evolves separate
electron and phonon temperatures, conserves their combined energy under local
exchange, routes deposited power to electrons, and relaxes phonons to the bath.
Its heat-capacity/conductivity partitions and linear transfer coefficient are
illustrative pending material calibration and nonlinear transfer laws. TDGL now
defaults to a gauge-compatible Heun predictor–corrector for gamma=0 small-step
accuracy. See `documents/TwoTemperatureTDGLBaseline_2026-09-21.md` for model
boundaries, example config, and stochastic-noise requirements.

Status: moving Gaussian optical deposition is implemented as a single-temperature exploratory profile. Absorbed power follows piecewise-linear SI waypoints, respects holes and device boundaries, and is evaluated at timestep midpoints. Signed vortex trajectories and commanded laser paths are exported and plotted. Configurable Gaussian material defects now reduce the local linear GL coefficient, permitting explicit pinning-versus-laser competition. Nearest laser-vortex separation is recorded without automatically claiming capture. Two-temperature dynamics, absorption calibration, calibrated defect strengths, and stochastic activation remain required for quantitative optical capture thresholds.

2026-09-21 optical-search checkpoint: `tools/optical_vortex_search.py` now
sweeps applied field, absorbed power, and path speed with the optional
two-temperature thermal branch. It requires a continuously tracked vortex,
beam-position correlation, bounded lag, superconducting electron temperature,
and no-laser/stationary/reverse controls before claiming a candidate. A
confirmed case receives a complete image/GIF/CSV/JSON report. The first short
pilots have **not** shown capture; see
`documents/OpticalTweezerTwoTemperatureSearch_2026-09-21.md` for setup and
remaining model limitations. The older statement above describes the original
single-temperature profile; the new experiment is the two-temperature path.

### Stage 6: SQUID and circuit devices

Support multiply connected geometries, gauge-invariant junction phase differences, loop fluxoids, self and mutual inductance, multiple terminals, and external lumped circuits. Validate `Phi0`-periodic critical-current modulation, circulating currents, integer fluxoid states, and symmetric and asymmetric SQUID limits.

### Stage 7: fluctuations and automated model selection

Add thermodynamically consistent stochastic terms for switching distributions and vortex activation. Select the cheapest valid model with declared rules using temperature, geometry, screening estimates, requested observables, and live conservation diagnostics. Record every promotion or demotion and its reason.

## Immediate development checkpoint

The current checkpoint is Stage 1. Completion requires:

- Partial contacts affect only their covered TDGL boundary nodes.
- Uncovered portions of the same edge remain covariantly insulating.
- Electrical current injection uses signed boundary flux rather than a voltage-search wrapper.
- Arbitrary source and sink masks on outer edges are validated.
- Integrated source, sink, and cross-section currents agree within a configured tolerance.
- Tests cover full-edge and partial-edge contacts, reversed current, zero current, gauge offsets, and incompatible contact geometry.

After this checkpoint, magnetic screening is the highest-priority physics change.
