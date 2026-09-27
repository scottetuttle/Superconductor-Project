# Physics guide

This guide explains the physical meaning of each SHS subsystem. Equations and numerical implementation details are collected in [MathematicsAndSoftwareGuide.md](MathematicsAndSoftwareGuide.md).

## Superconducting order parameter and TDGL

The complex field `psi` describes the local superconducting condensate. Its magnitude indicates how strongly superconducting the material is, while its phase carries superflow and vortex winding.

- `|psi|` near its local equilibrium value means a developed condensate.
- `|psi|` approaching zero marks a normal or strongly suppressed region.
- A phase change of approximately `2π` around a closed loop marks one flux quantum of topological winding.
- A vortex has a suppressed core and signed phase winding. An antivortex has the opposite winding.

TDGL lets `psi` relax and move under temperature, current, vector potential, boundaries, pinning, and its own nonlinear condensation energy. Gauge-covariant links ensure that changing electromagnetic gauge does not change physical observables.

The implemented regime is a generalized, dissipative TDGL model for dirty thin films near `Tc`. It is useful for vortex motion, phase relaxation, current crowding, perforated films, and thermal suppression. It becomes less reliable far below `Tc`, in ballistic or very clean materials, at atomic-scale weak links, or when microscopic quasiparticle distributions matter. Thermal noise is absent, so stochastic switching and thermally activated vortex entry are not predicted.

## Normal electrical current

The electrical subsystem solves charge continuity. It combines superconducting current calculated from `psi` with dissipative normal current driven by the electric potential. In current-drive mode, specified source and sink contacts inject equal and opposite current. A reference-potential constraint removes the arbitrary voltage offset.

The preferred conductivity model uses the material's normal-state conductivity directly, as expected by the selected TDGL normalization. A legacy condensate-depletion option reduces normal conductivity as `|psi|²` grows; this is phenomenological and should not be treated as microscopic quasiparticle transport.

This subsystem handles two-terminal thin-film transport and flow around holes well enough for controlled studies. It does not include Hall conductivity, charge imbalance, frequency-dependent conductivity, capacitive junction current, external circuit impedance, or a kinetic quasiparticle equation. Partial contacts are supported, but complex multi-terminal circuitry remains future work.

## Heat and hotspots

Temperature evolves by heat diffusion, local heating, and relaxation to a bath. Joule heating comes from the dissipative normal current only. External and laser heat sources add energy independently. Holes do not conduct heat.

This one-temperature model assumes electrons and phonons share one local temperature. It is appropriate when equilibration between them is fast compared with the behavior being studied or when a qualitative thermal suppression model is sufficient.

It cannot resolve ultrafast electron heating, delayed phonon response, phonon escape, quasiparticle diffusion, recombination, wavelength-dependent absorption, or condensation-energy release. Those omissions are central for quantitative SNSPD hotspots and fast optical vortex manipulation. Material properties are currently not a comprehensive temperature-dependent database.

## Applied magnetic field

A uniform perpendicular magnetic field is represented by a vector potential. Symmetric and Landau gauges describe the same magnetic field. Gauge-equivalence tests check that observables agree when the representation changes.

An applied field influences superconducting phase gradients, currents, vortex energetics, and flux through holes. The current interface is built around perpendicular fields on a two-dimensional film. Arbitrary three-dimensional magnets, tilted fields, and time-dependent electromagnetic radiation are outside the present scope.

## Magnetic screening and self-field

Current in the film produces its own magnetic field. SHS converts volume current to sheet current and evaluates a regularized thin-film Biot–Savart response. When self-field feedback is enabled, it updates the induced vector potential and repeats the coupled solve until the induced field and other coupled fields converge.

This captures quasi-static thin-film screening and weak Meissner response. It assumes a thin sheet, an effectively open surrounding space, and instantaneous magnetostatics. The regularization height and finite computational domain affect short-range and edge behavior. Strong screening, thick films, multilayers, magnetic substrates, full three-dimensional fields, retardation, displacement current, and radiation are not quantitatively covered. Strong-screening results require refinement and comparison against an independent solver.

## Vortices

In the full TDGL model, vortices are features of `psi`, not separate particles. They can deform, interact with boundaries and holes, respond to currents and heat, and in principle annihilate with opposite winding when the grid and timestep resolve the event. Diagnostics locate signed winding on grid plaquettes and link detections over saved frames.

Tracking can lose identities when vortices move more than the configured maximum between frames, overlap, pass through regions where `|psi|` is too small, or form dense clusters. Detected counts also depend on the amplitude floor and mesh resolution. Initial vortex seeds are numerical initial conditions.

The separate reduced vortex model treats each vortex as an overdamped point acted on by vortex interactions, Lorentz force, thermal force, and pinning. It is fast and useful for protocol screening. It cannot create, destroy, split, deform, or resolve vortices and should not replace TDGL verification.

## Holes, rings, and fluxoids

Circular and rectangular holes remove cells from the active superconducting, electrical, and thermal domains. Links crossing a hole boundary are omitted, giving no superconducting, normal-current, or thermal leakage through the void.

A closed superconducting contour can carry integer phase winding. SHS reports both the topological winding and a London fluxoid made from magnetic flux plus supercurrent circulation. Agreement is a strong consistency check when `|psi|` remains nonzero along the contour.

Rings and perforated films are implemented. A SQUID is not yet complete because the engine lacks validated Josephson weak-link laws, junction capacitance/resistance models, external circuit coupling, and `Φ0`-periodic critical-current validation.

## Static pinning

Pinning sites locally reduce the linear GL coefficient with Gaussian profiles. This lowers condensation energy and makes the region favorable to a vortex core. Multiple sites add up to a configured maximum.

This is a controllable effective defect landscape. Its strength is dimensionless and is not yet derived from a measured defect composition, thickness change, local `Tc`, or microscopic disorder. Quantitative depinning currents therefore require calibration.

## Moving laser

The laser is an absorbed Gaussian power distribution moving along piecewise-linear waypoints. The power envelope can ramp up or down at each waypoint. Only power landing on active material is deposited; power outside the device or in a hole is lost rather than redistributed.

The laser influences `psi` through the heat equation: it raises temperature, weakens the condensate, and creates a moving low-condensation-energy region that may attract or release vortices. The code records deposited power, beam position, vortex paths, and nearest beam–vortex separation.

This is a thermal optical model. It omits electromagnetic coupling of light, wavelength and polarization, absorption depth, electron–phonon delay, nonequilibrium quasiparticles, stochastic activation, and experimentally calibrated optical efficiency. It can demonstrate a physically sensible mechanism and compare protocols, but cannot yet predict a laboratory capture threshold.

## Coupling among systems

The major feedback loop is:

```text
temperature and pinning change psi
psi and vector potential determine supercurrent
supercurrent enters charge continuity and magnetic screening
normal current produces Joule heat
heat changes temperature
temperature changes psi again
```

The solver iterates this loop within each timestep. Relaxation improves fixed-point stability but does not guarantee that a large physical timestep is accurate. A converged coupled iteration can still be temporally under-resolved.

## Current physics priorities

The most valuable physical additions are a validated semi-implicit TDGL integrator, an implicit or multirate thermal integrator, two-temperature electron–phonon physics, calibrated optical absorption, Josephson weak-link models, and stronger independent magnetic-screening benchmarks. These additions should keep the existing conservation and gauge checks as acceptance criteria.

