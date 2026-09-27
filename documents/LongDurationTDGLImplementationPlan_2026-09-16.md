# Long-duration TDGL implementation plan — 2026-09-16

## Objective

Increase simulated physical duration without weakening gauge covariance or silently skipping fast vortex, electrical, optical, or thermal dynamics. The current explicit TDGL update is diffusion-limited: refining the mesh makes its stable timestep scale approximately with the square of the cell size, so long experiments require many internal substeps.

## Recommended next implementation

Implement a sparse, gauge-covariant IMEX integrator as a second selectable TDGL backend. Keep the current explicit solver as the reference implementation.

For each accepted step:

1. Assemble the link-variable covariant Laplacian on active superconducting cells, including insulating holes and configured normal contacts.
2. Advance that stiff linear diffusion term implicitly with backward Euler for the first implementation.
3. Advance the local temperature coefficient, cubic order-parameter term, scalar-potential temporal link, pinning suppression, and generalized TDGL correction consistently with the selected IMEX scheme.
4. Estimate temporal error by step doubling. Reject the step and restore all fields when the error exceeds tolerance; grow it conservatively when the solution is smooth.
5. Cache the sparse topology for a fixed mesh and geometry. Reuse a factorization while the vector potential and timestep remain unchanged; otherwise update matrix values without rebuilding connectivity.

Backward Euler is the useful first milestone because it is robust and straightforward to validate. A second-order scheme such as Crank–Nicolson/Adams–Bashforth or an additive Runge–Kutta IMEX method should follow only after the first-order path agrees with the explicit reference.

## Coupled timescale strategy

- Resolve TDGL and electrical changes at the fastest accepted timestep.
- Update the laser position from continuous physical time, including on rejected and subdivided intervals.
- Permit the thermal solver and magnetic screening iteration to run at integer multiples of the TDGL step only when monitored temperature, current-continuity, field, and order-parameter errors stay below configured thresholds.
- Force synchronization near vortex creation/annihilation, contact transients, large optical gradients, phase slips, or failed nonlinear/screening convergence.
- Add checkpoint/restart files containing physical time, all fields, adaptive-controller history, vortex identities, laser protocol state, and configuration provenance. This makes long runs recoverable and reproducible.

## Required validation before making it the default

- Temporal convergence against the existing explicit solver on the same mesh.
- Exact gauge-equivalence regression for the supported gauges.
- Uniform zero-field relaxation compared with the analytic local TDGL solution.
- Vortex position, winding, and fluxoid agreement during stationary and moving-laser tests.
- Current-continuity and terminal-current balance under transport bias.
- Thermal energy balance and agreement between single-rate and multirate runs.
- Step rejection and checkpoint/restart equivalence.
- Measured wall-clock speedup at fixed physical error, rather than at fixed step count.

## Configuration surface

The integrator should be selected explicitly, for example `tdgl.integrator = "explicit" | "imex_be"`. Adaptive controls should use physical error tolerances, minimum and maximum normalized timestep, growth/shrink limits, and a rejection limit. Multirate thermal and magnetic intervals should remain bounded by error monitors rather than fixed skip counts alone.

## Limits this does not remove

Longer stable steps do not add missing microscopic physics. Electron–phonon nonequilibrium, quantitative optical absorption, material-specific pinning, Josephson weak-link boundary laws, full Maxwell retardation, and calibrated noise remain separate model improvements. The integrator must report when a configured experiment exceeds the validity regime of its selected physics model.

