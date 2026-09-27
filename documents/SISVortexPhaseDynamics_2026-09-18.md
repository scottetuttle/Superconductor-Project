# SIS vortex and phase-dynamics investigation — 2026-09-18

## What the focused run represents

The SIS microwave runner is a lumped RCSJ junction with a static London phase
texture imposed by pinned electrode vortices. It is not a spatial sine-Gordon
or TDGL evolution of the junction and electrodes. Vortex positions never
change during this calculation.

The finite-width model first evaluates

`F = width_average(exp(i * theta_vortex(y)))`.

It then replaces the spatial junction by a lumped critical current
`Ic_eff = Ic * abs(F)` and phase offset `arg(F)`. The displayed local current is
reconstructed afterward as `sin(theta_vortex(y) + phi(t))`; it does not feed
back as independent spatial dynamics.

## Focused opposite-vortex result

For `sis_shapiro_focused.json`, the opposite-vortex texture produces:

- bare critical current: 2.000 microamperes;
- effective critical current: 1.76333 microamperes;
- static phase offset: -2.61323 radians;
- microwave frequency: 5.000 GHz;
- mean Josephson frequency at the animated 2 microampere bias: 9.99957 GHz;
- locking ratio: 1.99991, identifying the second Shapiro step;
- small-signal plasma frequency: 26.0498 GHz;
- McCumber parameter: 0.42863.

The lower-current step boundaries primarily follow from the 11.8% reduction
of effective critical current. The static phase offset changes the phase origin
but does not itself create vortex motion.

The voltage spectrum contains the 5 GHz drive, the 10 GHz locked Josephson
rotation, and nonlinear harmonics/mixing components at approximately 20, 25,
and 30 GHz. The two visibly different rotation speeds within a drive cycle are
the acceleration and deceleration of one nonlinear phase variable. On the
second step it advances twice per microwave period. The response near the
26 GHz plasma scale sharpens this nonuniform motion.

## Why the vortex appears to spin

The film animation adds `phi(t)` uniformly to the right electrode and displays
phase modulo 2 pi. A cyclic color map therefore rotates its branch cut and
color contours around a fixed vortex. This is a visualization of global phase
advance superposed on a static winding texture. It is not evidence that the
vortex core or circulating superflow physically rotates through the film.

## New diagnostics

The focused runner now saves, for the selected GIF bias:

- the full steady-state transient as CSV;
- a normalized voltage FFT as CSV;
- separate FFT columns and plots for phase rate, voltage, Josephson current,
  quasiparticle/resistive current, capacitive current, and centerline local
  Josephson current;
- a JSON frequency summary;
- expanded phase, phase-rate, RCSJ current-balance, and FFT plots;
- a time-versus-junction-position local-current kymograph;
- the complete input configuration inside `summary.json`.

These outputs distinguish global phase wrapping, Shapiro locking, nonlinear
harmonics, and static spatial vortex imprinting.

In the RCSJ equations, the quantity labelled `resistive_current_A` is
`V/R = (Phi0 / 2 pi R) * phase_rate`. Its spikes must therefore coincide with
the phase-rate and voltage spikes; they are the same dynamical variable scaled
by constants. The resistance represents a phenomenological dissipative
quasiparticle/shunt channel across the junction. It does not mean that a region
of either superconducting electrode temporarily becomes normal. The present
model has no spatial gap suppression, barrier breakdown, microwave heating, or
vortex-core dissipation capable of making that claim.

## Model boundary

Claims about vortex motion, fluxon propagation along the barrier, spatial
plasma modes, or vortex-junction back-action require a long-junction
sine-Gordon model or a junction-resolved TDGL/electromagnetic calculation.
The present reduced model can study how a prescribed pinned-vortex phase
texture shifts effective critical current and lumped Shapiro response, but it
cannot predict dynamical vortex rotation.
