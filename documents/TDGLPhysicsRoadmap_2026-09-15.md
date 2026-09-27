# TDGL physics cleanup and roadmap - 2026-09-15

## Runnable state

The simulator now has a coherent small-step generalized TDGL path. `configs/simulations/nbn_tdgl_small_step.json` runs a 100 fs example as ten 10 fs coupled steps. It uses the electrical scalar potential in the gauge-covariant time derivative, the existing link-variable spatial operator, explicit stability-controlled substeps, normal-contact order-parameter boundaries on the transport sides, and an accurate sparse electrical solve.

The generalized update implements the Kramer-Watts-Tobin `gamma` correction. Setting `gamma = 0` recovers the established SHS relaxation equation. A spatially constant voltage changes only the global phase in regression tests; it does not change the order-parameter amplitude or supercurrent. Normal contacts impose `psi = 0`; insulating sides impose a zero gauge-covariant normal derivative.

The SHS normalization uses the material coherence length, the configured dimensionless GL amplitude, the time scale `pi*hbar/(8*k_B*Tc)`, and the scalar-potential scale `Phi0/(2*pi*time_scale)`. The potential scale follows from gauge covariance for the chosen SHS time scale. It must not be mixed with a conductivity-based TDGL time scale unless all coefficients are converted together.

## Physics that remains incomplete

The current implementation is a useful small-domain TDGL backbone, but it is not yet a complete quantitative transport model.

1. The electrical equation solves charge continuity using a voltage drive, but terminal current injection is not yet a self-consistent TDGL boundary condition. `current.value` remains unused.
2. A normal contact can currently occupy only a complete mesh side. The geometry contains contact masks, but TDGL boundary conditions do not yet apply to arbitrary terminal segments.
3. The vector potential is externally held. There is no self-consistent magnetic-field or Maxwell solve, so vortex self-fields and screening currents are absent; configured `kappa` is consequently inactive.
4. Normal conductivity and phenomenological GL coefficients need calibration against the precise nondimensional TDGL convention chosen for an experiment. The implemented linear coefficient is the existing SHS `1 - T/Tc` convention.
5. Thermal coupling includes Joule heating but does not yet include every possible order-parameter relaxation or substrate/interface heat term.
6. The explicit TDGL integrator is appropriate for small steps and grids. Stiff, long-duration simulations need an implicit or IMEX method with a temporal error estimate.

## Next implementation stages

The next physics change should map normal-contact conditions to actual geometry masks and add the accompanying scalar-potential/current boundary condition. That will make transport terminals internally consistent and allow a real current-driven mode.

After terminal consistency, implement a self-consistent vector-potential solve and validate flux quantization, gauge invariance, vortex equilibrium, and current conservation. Only then should results be treated as quantitative magnetic TDGL predictions.

For speed, cache the sparse electrical matrix topology and factorization whenever conductivity is unchanged. For changing conductivity, reuse the sparsity pattern and add a preconditioned Krylov path. Replace explicit TDGL stepping with an IMEX or semi-implicit covariant Laplacian solve, and control the outer timestep from an embedded temporal error estimate rather than coupling residual alone.

## Validation completed

- Scalar-potential phase rotation and gauge-offset invariance.
- Generalized update equivalence at `gamma = 0` and finite behavior for nonzero `gamma` at small steps.
- Gauge-covariant insulating boundaries and zero-order-parameter normal contacts.
- Sparse electrical solution agreement with SOR and equation-defect verification.
- End-to-end execution of the conservative small-step configuration.

The adaptive benchmark now reports current-continuity, Joule-power, thermal-energy, minimum-temperature, and maximum-order-parameter diagnostics. These should become acceptance metrics for future method comparisons rather than relying on runtime and fixed-point residual alone.
