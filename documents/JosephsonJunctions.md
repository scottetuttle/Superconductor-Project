# Josephson junctions

## SIS RCSJ backend

SHS includes a lumped SIS junction backend implementing

```math
I(t)=I_c\sin\varphi+\frac{V}{R_n}+C\frac{dV}{dt},
\qquad
V=\frac{\Phi_0}{2\pi}\frac{d\varphi}{dt}.
```

It supports current and voltage bias, DC plus sinusoidal AC drive, capacitance, quasiparticle resistance, phase dynamics, dissipated power, and uniform-junction Fraunhofer modulation. `critical_current_A` may be supplied from measurement or derived from the weak-coupling Ambegaokar–Baratoff relation using temperature, `Tc`, and normal resistance.

Run the supplied current sweep with:

```powershell
python tools/sis_junction_diagnostics.py --config tools/config/sis_junction_sweep.json
```

The example has `Ic = 2 µA`. Its ±1.5 µA cases settle to numerically zero average voltage, while ±2.5 µA and ±4 µA enter symmetric finite-voltage running-phase states. Output includes phase, voltage, Josephson current, component-resolved CSV histories, and an I–V curve.

This is the standard complete classical lumped model for an SIS junction. It is not a microscopic spatial tunnel calculation: it does not solve Usadel/BCS Green functions in the electrodes, quantum phase escape, frequency-dependent environmental impedance, self-heating, or nonequilibrium quasiparticles.

The electron–phonon two-temperature model should be added next for switching, retrapping, hysteresis caused by heating, and optical response. It was not required before implementing the Josephson relations themselves. Junction Joule power is already exported so it can later enter the electron energy equation.

## Vortices near an SIS junction and microwave response

Run:

```powershell
python tools/sis_vortex_microwave_diagnostics.py --config tools/config/sis_vortex_shapiro.json
```

The finite-width junction is sampled across its aperture. Applied perpendicular field supplies a linear phase gradient. Configured vortices in the left or right electrode add signed London phase textures. The complex aperture average determines the effective critical-current suppression and phase offset used by the RCSJ dynamics. This captures interference from stationary vortices but does not resolve their cores, forces, creation, or motion.

The microwave is represented as an AC current added to the DC bias. This is appropriate for a junction driven by a coupled microwave line after its delivered current amplitude is known. It is not a Maxwell simulation of an incident electromagnetic wave, antenna, or package.

At frequency `f`, Shapiro plateaus are expected at

```math
V_n=n\Phi_0 f.
```

The diagnostic plots the ordinary I–V curve and `V/(Phi0 f)` with integer reference lines. It also records the nearest step and locking error for every bias, plots unwrapped phase, and shows the vortex-induced current texture across the aperture. The supplied 5 GHz example has a step spacing of 10.339 microvolts and produces multiple integer-locked plateaus.

In the supplied comparison, no vortex gives `Ic = 2 µA`; a nearby +1 vortex gives an interference amplitude near 0.499, phase offset near 0.404 rad, and effective `Ic` near 0.998 µA. The configured opposite-vortex pair gives amplitude near 0.882, phase offset near -2.61 rad, and effective `Ic` near 1.76 µA. These are predictions of the configured London aperture approximation, not universal vortex-junction values.

The expanded configuration sweeps four delivered microwave amplitudes (`0`, `0.7`, `1.4`, and `2.1 µA`) for every vortex layout and every DC bias. `shapiro_sweep_matrix.png` places all twelve normalized I–V responses in one grid. The zero-microwave column is the required control: isolated crossings of integer voltage are not Shapiro plateaus unless a finite bias interval remains locked.

Each vortex layout receives a static `*_film_phase.png` image showing the idealized phase of the two electrodes, the insulating barrier, and vortex positions. At the selected microwave amplitude and DC bias, an optional GIF animates global electrode phase and the local Josephson-current distribution across the aperture. These images deliberately omit `|psi|` and vortex cores because the lumped SIS model does not calculate them. Use full TDGL when core deformation or vortex motion is the research question.

GIF generation is the most expensive output operation. Set `output.gif_enabled` to false for broad sweeps, then enable it after choosing a scenario. `gif_microwave_amplitude_A` and `gif_bias_current_A` select the single operating point animated for each vortex layout.

For detailed plateau measurement, run `tools/config/sis_shapiro_focused.json`. It simulates 40 ns, or 200 periods at 5 GHz, uses 12,000 saved points and 0.1 microamp DC spacing, and averages only the latter half. The phase-aware GIF shows the final two microwave periods with 60 frames per period instead of compressing the complete run. The example resolves locked ranges of roughly 1.1–1.6 microamps for `n=1`, 1.7–2.2 microamps for `n=2`, and 2.3–2.6 microamps for `n=3`.

High-frequency integration is protected by `run.microwave_steps_per_cycle`. The ODE solver uses the smaller of `max_step_s` and one microwave period divided by this value. Increasing frequency therefore decreases the internal maximum step automatically. Saved output sampling remains separately controlled by `run.samples`.

## Implemented model

SHS now supports a spatially resolved vertical TDGL weak link. The configured strip subtracts `coefficient_suppression` from the local linear GL coefficient. A value larger than the surrounding positive coefficient makes the strip locally normal-favoring while superconducting banks remain on each side. This represents an effective Dayem bridge or resolved SNS-like weak region.

It is not yet a microscopic tunnel barrier. There is no standalone sinusoidal current-phase boundary condition, capacitance, RCSJ dynamics, quasiparticle tunneling model, or external microwave circuit.

## Measurements

The diagnostic system records:

- terminal voltage across the source and sink contacts;
- transported current;
- gauge-invariant phase difference between the two junction banks;
- mean order-parameter amplitude inside the weak link;
- left- and right-bank amplitudes;
- ordinary temperature, current, field, and convergence diagnostics.

The phase measurement subtracts the vector-potential line integral, so changing gauge should not change its physical value.

## Current sweep

Run:

```powershell
python tools/josephson_diagnostics.py --config tools/config/josephson_current_sweep.json
```

The tool creates `josephson_sweep.png` containing the I–V curve, numerical differential resistance, phase response, and weak-link amplitude. `josephson_transients.png` shows voltage and phase versus time for every bias. Standard field maps, kymographs, current maps, CSV histories, and summaries are also retained.

The supplied short validation sweep converged at all five currents from -10 nA to +10 nA. It produced an odd, symmetric voltage and phase response and recovered every requested current. The response was nearly linear over 200 fs; this validates measurement symmetry and wiring but does not yet demonstrate a nonlinear Josephson critical current. Longer relaxation and continuation from one bias point to the next are required for a defensible switching curve.

The repeated +10 nA case carried approximately 0.973 nA as supercurrent and 9.027 nA as normal current at the center section: about 9.7% superconducting and 90.3% normal. This occurs because the current model retains the full normal-state conductivity throughout the device and uses normal TDGL contacts. The weak-link amplitude decreased only about 2.8% from its initial value during the short run. The result is therefore consistent with a resistive metallic weak link that has not equilibrated, rather than an SIS tunnel junction.

A longer 1 ps experiment confirms that time is the dominant limitation in this setup. At 10 nA magnitude, the supercurrent share increased from about 9.7% at 200 fs to 26.4% at 1 ps, while the phase-drop magnitude increased from approximately `3.70e-5` to `1.41e-4` rad. Doubling the device dimensions over the short duration produced much less change because added space does not immediately create the required phase gradient. These values still describe a transient, not a converged DC branch.

Field and current images mark the configured weak-link strip with translucent cyan boundaries. Each case also includes `junction_profile.png`, which plots the cross-device mean `|psi|`, supercurrent, and normal current. `josephson_transients.png` includes the normal and superconducting current buildup versus time. The marker represents a material weak region, not an empty geometric gap.

## Research interpretation

Use current reversal to test symmetry and inspect whether voltage remains near zero below a stable supercurrent threshold. A transition to time-varying phase and nonzero average voltage is the expected signature of a running Josephson state. Estimate critical current only after timestep, mesh, link-width, suppression-strength, and run-duration convergence.

The current sweep initializes every bias independently. That is appropriate for equilibrium comparison but does not capture hysteresis. A future continuation mode should reuse the final state from the preceding bias, sweep upward and downward separately, discard a configurable transient, and time-average voltage and phase velocity. AC bias is also required to study the AC Josephson relation and Shapiro steps.

## Configuration controls

- `josephson.enabled`: activates the weak link.
- `josephson.center_x_m`: horizontal center of the vertical strip.
- `josephson.width_m`: physical strip width; it should span several mesh cells.
- `josephson.coefficient_suppression`: dimensionless reduction of the local GL coefficient.

The implementation currently supports one vertical full-width junction. Multiple, angled, shaped, or material-resolved junctions remain future geometry work.

An SIS implementation requires a different model: strongly reduced normal conductance in the barrier together with an interface Josephson coupling or calibrated tunnel-boundary condition. Simply setting conductivity to zero would also block the present finite-volume current path and would not by itself create Josephson tunneling.

## Stateful current-reversal protocol

Run the zero-current relaxation, positive-bias stabilization, and negative-current reversal experiment with:

```powershell
python tools/josephson_reversal_protocol.py --config tools/config/josephson_reversal_protocol.json
```

The supplied configuration runs independent 5 nA and 10 nA experiments. Every experiment starts from a fresh device, stabilizes at zero current, continues the same state at positive current, and then reverses current without reinitializing `psi`. This preserves the phase and condensate history needed to observe reversal dynamics.

A stage may end after its minimum duration when weak-link amplitude and voltage drift are small and either the phase is stationary or its phase-velocity has stabilized. The latter condition prevents a valid running Josephson state from waiting forever for a constant phase. Every stage also has a hard maximum duration, and its recorded `stop_reason` distinguishes stabilization from exhaustion or coupling failure.

`protocol.csv` records every physics step. `reversal_protocol.png` shows commanded current, phase difference, voltage, and the normal/superconducting current split. Optional GIF frames show order-parameter amplitude and phase during both switches.

The expensive field snapshot, Biot–Savart diagnostic, and GIF capture are sampled every `output.sample_every` steps. The coupled physics still advances at every `dt_s`. For preliminary searches, set `gif_enabled` to false and increase `sample_every`; after selecting a useful current, enable the GIF for that case. This reduces visualization cost without changing the physical integration.

The junction uses a dedicated 100 nm × 50 nm geometry (`20xi × 10xi`) with a 201×101 mesh, giving approximately `0.1xi` spacing. Increasing only `nx` and `ny` refines the same physical device; it does not enlarge it. A 200×200 mesh on the former 50 nm square had `0.05xi` spacing, about 40,000 electrical unknowns, and required roughly 69 explicit TDGL substeps per 10 fs physical step. That combination explains hour-scale protocol runs.
