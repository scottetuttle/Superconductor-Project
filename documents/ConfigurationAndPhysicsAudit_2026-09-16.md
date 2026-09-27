# Configuration propagation and physics audit

## Outcome

The production builder does load geometry, material, initial temperature, current, thermal settings, electrical settings, TDGL controls, magnetic screening controls, laser paths, pinning sites, and boundary conditions from the simulation JSON. Targeted tests now verify geometry/material construction, initial state, thermal and TDGL model construction, laser position/power/width, pinning landscapes, and independence of separately loaded configurations.

The main source of confusion is the diagnostics layer. `tools/tdgl_diagnostics.py` deliberately resolves an experiment from several files and gives some values higher precedence than the simulation JSON. Every new summary now contains `configuration_provenance` and `effective_simulation_config`, and the runner prints when its requested duration differs from the simulation duration.

## Configuration precedence

From lowest to highest precedence:

1. `src/shs/config/defaults.json` supplies omitted optional settings.
2. The simulation JSON supplies geometry/material references and production physics settings.
3. The diagnostics JSON `run` section supplies timestep and step count for that diagnostic run.
4. The diagnostics `transport` section supplies drive mode and target current.
5. Diagnostic sweep entries override their corresponding values for an individual case.

Applied magnetic field and initial vortices are diagnostic inputs, not fields in the simulation schema. A diagnostics run therefore does not mean “execute the simulation JSON unchanged.” For example, `100` steps at `1e-14 s` runs for `1e-12 s` even if the simulation JSON duration is `1.6e-13 s`.

## Settings that currently affect computation

- Geometry filename, width, height, thickness, mesh, holes, and contacts.
- Material filename and its `Tc`, coherence length, penetration depth, conductivity/resistivity, thermal conductivity, heat capacity, and thickness.
- Initial and bath temperatures, relaxation rate, and thermal integration controls.
- Current or voltage drive, contact names, reference voltage, conductivity model, and electrical solver controls.
- TDGL `u`, `gamma`, normalization, temperature coefficient, scalar potential, timestep limit, and stability factor.
- Laser enable, absorbed power, width, position, time, repetition, and waypoint power fraction.
- Static pinning enable, sites, widths, strengths, and total suppression cap.
- Self-field enable, fixed-point tolerance, iteration limit, relaxation, and source stride.
- Coupling convergence controls, field scales, and common numerical safeguards.

## Parsed settings that do not yet implement their advertised physics

- `tdgl.kappa` is validated and stored but is not used directly by the TDGL evolution. Material `lambda` and `xi` enter scaling and magnetic calculations instead.
- `electromagnetic.include_displacement_current` is stored but no displacement-current evolution exists.
- Material `gl_alpha`, `gl_beta`, and `tdgl_u` are copied into the material map but the normalized TDGL solver uses the simulation-level normalized coefficient and `tdgl.u`.
- The electromagnetic model's older `ground_voltage` setting is not part of the active simulation configuration path.

These fields should eventually be activated or renamed/deprecated. They must not be used to infer capabilities from a configuration file today.

## Highest-impact physical limitations

### 1. Accessible physical time

Explicit gauge-covariant TDGL and explicit thermal diffusion require very small internal steps. Experimentally meaningful optical manipulation can take nanoseconds to microseconds, while current full-TDGL examples cover femtoseconds to picoseconds. This prevents direct qualification of the successful reduced-model protocol. An IMEX/semi-implicit TDGL method, implicit thermal solve, and error-controlled multirate stepping are the highest-value numerical changes.

### 2. One-temperature optical and hotspot physics

Optical and Joule power enter one temperature with constant material coefficients. There is no electron temperature, phonon temperature, electron-phonon transfer, substrate spreading/escape, interface resistance, wavelength-dependent absorption, nonequilibrium quasiparticle diffusion, or order-parameter relaxation heat. This is the largest physical limitation for laser manipulation and hotspot recovery.

### 3. Vortex stability near boundaries

Seeded vortices are imposed phase/core profiles rather than equilibrated states. Normal contacts set the order parameter to zero and can become strong vortex escape channels. Finite edges also attract vortices through image/boundary physics. A vortex placed a few coherence lengths from an edge may leave even with zero current. Experiments must include laser-off equilibration and control runs with zero/reversed current and vortex charge.

### 4. Thin-film magnetic environment and finite domain

Self-field coupling uses a quasi-static free-space sheet-current convolution on the same finite rectangular domain. It does not represent a separately enlarged vacuum domain, arbitrary magnetic media, leads, multilayers, or full nonlocal device/substrate electrodynamics. When the Pearl length is much larger than the device, quantitative vortex energies and edge interactions require special care.

### 5. Normal transport and energy conversion

The constant normal conductivity option is appropriate to the standard generalized-TDGL current equation, but conductivity is not temperature-, frequency-, or nonequilibrium-dependent. Condensate relaxation and vortex motion do not return their dissipated free energy to the heat equation. This limits switching, retrapping, and hotspot predictions.

### 6. Pinning and fluctuations

Gaussian pinning sites are phenomenological reductions of the linear GL coefficient. Their strengths are not calibrated pinning energies. Thermal Langevin terms are absent, so depinning and switching are deterministic rather than statistical.

### 7. Model regime and device breadth

The normalized TDGL model is most defensible near `Tc` in dirty gapless-like regimes. It does not yet contain microscopic low-temperature dynamics, calibrated Josephson weak-link interfaces, RCSJ/circuit coupling, or arbitrary multilayer/drawn outer geometry. Rectangular films with rasterized holes are supported.

## Recommended next work

Implement the electron-phonon/substrate thermal model and an implicit or multirate time integrator together. The thermal model supplies the correct force landscape; the time integrator makes its experimental timescale reachable. Before expanding the UI, also add warnings for inactive configuration fields and a machine-readable resolved-configuration artifact shared by all runners.
