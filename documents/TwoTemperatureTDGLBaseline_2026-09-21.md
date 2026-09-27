# Two-temperature and TDGL time-integration baseline — 2026-09-21

## Implemented

The optional `thermal.model = "two_temperature"` adds a film-phonon temperature
field. `fields.temperature` is the **electron** temperature and remains the
temperature that controls TDGL coefficients and normal conductivity. Joule,
laser, and externally supplied heating enter the electron balance. The
configured material heat capacity and conductivity are partitioned between
electrons and phonons by explicit fractions. Finite-volume diffusion respects
holes and insulating boundaries. Prescribed boundary temperature clamps both
subsystems; convection and prescribed heat flux act on phonons. Phonons relax
toward the substrate bath at `phonon_escape_rate_per_s`.

The initial exchange model is

`Ce dTe/dt = ... - G(Te-Tp)`,

`Cp dTp/dt = ... + G(Te-Tp) - Cp (Tp-Tbath)/tau_escape`.

The local electron–phonon transfer is integrated analytically at each thermal
substep, so it conserves `Ce*Te+Cp*Tp` exactly in the absence of other terms
and does not impose a stiff explicit exchange timestep. Diffusion and sources
retain the existing bounded explicit substeps. Coupled fixed-point iterations
restart both temperatures from the accepted time level and include phonon
temperature in convergence checks. Physics diagnostics now report separate
electron and phonon energies plus their sum. TDGL diagnostics report both
temperatures and their maximum difference.

TDGL now defaults to `time_integrator = "heun"`. It uses a covariant Euler
predictor and averaged old/predicted spatial right-hand sides, parallel
transporting the predicted right-hand side back across the scalar-potential
temporal link. The old `"euler"` option remains for comparison. A uniform
relaxation regression demonstrates lower error against a refined solution.
This is a second-order improvement for the gapless `gamma=0` equation with
frozen coupled inputs. Nonzero Kramer–Watts–Tobin gamma still uses the existing
nonlinear local update with an averaged spatial RHS; its global order has not
been established. The coupled electrical, thermal, and magnetic pass remains
operator-split/fixed-point, so this change alone does not prove second-order
accuracy for the whole multiphysics solver.

## Example

`configs/simulations/nbn_two_temperature_hotspot.json` enables the model with
illustrative coefficients; `tools/config/nbn_two_temperature_hotspot.json`
selects diagnostics. Run:

`python tools/tdgl_diagnostics.py --config tools/config/nbn_two_temperature_hotspot.json`

The ten-step smoke run converged. Its final maximum electron and phonon
temperatures were approximately 15.3131 K and 14.5001 K, respectively. This
is an implementation check, not a calibrated NbN prediction.

## Physical limits and next validation

The current material file has a single effective heat capacity and thermal
conductivity. Their partition fractions and linear transfer coefficient are
therefore **illustrative inputs**, not measured NbN material functions.
Quantitative hotspot work needs independently measured `Ce(Te)`, `Cp(Tp)`,
`ke(Te)`, and `kp(Tp)`, a calibrated electron–phonon transfer law (often a
power difference rather than `G(Te-Tp)`), and substrate-interface conductance.
It also needs order-parameter condensation/recombination energy bookkeeping
and optical absorption calibration. The two-temperature model assumes each
subsystem can be assigned a local equilibrium temperature; it cannot resolve
nonthermal quasiparticle distributions. [Nanowire hotspot study](https://www.nature.com/articles/s41467-022-32719-w)
uses separate electron and phonon balances with nonlinear transfer and escape
laws; its measured timescales are a useful future calibration target.

Stochastic fluctuations remain **disabled and unimplemented**. Adding
arbitrary Gaussian noise now would yield mesh- and timestep-dependent switching
rates. The next implementation should derive discrete TDGL and thermal noise
covariances from the fluctuation–dissipation relation using each control-volume
size and local electron temperature, preserve gauge covariance, and keep the
same random increment across failed fixed-point iterations and adaptive
retries. Tests should check seed reproducibility, variance scaling with cell
volume and timestep, equilibrium distributions, and switching statistics.
[Berger's fluctuation–dissipation analysis](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.75.184522)
shows why the normalization of Langevin terms matters in nonuniform wires.

The long-duration optimization remains a separate step: a gauge-covariant
semi-implicit or IMEX diffusion solve with embedded temporal-error control is
needed before large meshes and long physical durations can be efficient.

