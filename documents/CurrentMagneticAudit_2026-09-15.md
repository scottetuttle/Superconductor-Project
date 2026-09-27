# Current and magnetic-field audit - 2026-09-15

## Corrections and additions

The previous electromagnetic module contained declarations but no field calculation. Vector-potential arrays remained zero, `magnetic_field_x/y` remained zero, and `current.value` did not control transport. The electrical solver itself consistently solves `div(Jn + Js) = 0` with `Jn = -sigma grad(V)` and reports link currents in A/m². Its face coefficient conversion was checked: the Poisson coefficient contains `sigma/dx²`, and multiplying it by `dx² E` recovers `sigma E` without an extra length factor.

The diagnostics now provide three gauges for a uniform applied perpendicular field. Numerical curl tests recover the requested field in every gauge to approximately machine precision. The Landau-gauge initial order parameter receives the corresponding gauge phase transformation.

A separate 50 nm square transport geometry places current contacts on the complete left and right edges. Current mode integrates `Jx` over a vertical cross section and film thickness, then tunes terminal voltage using rollback-based trials. This is a practical current-controlled wrapper around the existing Dirichlet electrical solve. The more direct future formulation is a Neumann scalar-potential boundary at each normal terminal.

The current-generated magnetic field is calculated at a nonzero height using the thin-film Biot-Savart integral with sheet current `K=Jd`. An FFT convolution evaluates the discrete uniform-grid integral efficiently. It is currently diagnostic-only and is not added back to the TDGL vector potential.

## Research basis

The implementation follows the pyTDGL background equations for the covariant derivative, supercurrent, normal-terminal conditions, magnetic/vector-potential scales, and current conservation: https://github.com/loganbvh/py-tdgl/blob/main/docs/background.rst

The thin-film induced-vector-potential expression and the need to iterate induced field and current to self-consistency follow the pyTDGL screening formulation: https://py-tdgl.readthedocs.io/en/latest/notebooks/screening.html

Those sources prescribe `psi=0` and a normal derivative of scalar potential proportional to applied current at normal terminals. The SHS voltage controller reaches the requested net current but is an intermediate implementation until the sparse electrical operator supports explicit Neumann terminal fluxes.

## Validation experiment

The generated run in `benchmark_results/tdgl_diagnostics_2` compares 0 T and 1 T at a target current of 10 nA for 100 steps of 10 fs.

- Both cases converged and reached the target current within the configured 0.1% tolerance.
- The applied-field curl error was zero at reported precision.
- The zero-field final relative continuity defect was about `3.1e-15`; the 1 T case was about `1.4e-16`.
- The zero-field transverse current-profile nonuniformity was about `0.0027`.
- At 1 T, that metric approached `1.0`, and the diagnostic self-field reached about 7.3 mT even though the net transport current remained 10 nA. This reflects large local field-driven supercurrents and means magnetic screening should not be neglected for quantitative interpretation of this case.

The 1 T trajectory required substantially more coupling work than the zero-field trajectory. This is physically informative, but the result remains qualitative because induced magnetic feedback is absent.

## Remaining work

1. Add Neumann current-flux boundary terms directly to the sparse scalar-potential system and enforce signed terminal-current balance.
2. Store applied and induced vector potentials separately.
3. Iterate the induced vector potential and TDGL current to a configured screening tolerance during every physical step.
4. Add fluxoid diagnostics on user-selected closed contours. A nonzero vortex-free fluxoid is the key test for whether screening can be neglected.
5. Reconcile the complete TDGL time, voltage, normal-conductivity, and current normalization with one selected convention before quantitative comparison with experiment.
