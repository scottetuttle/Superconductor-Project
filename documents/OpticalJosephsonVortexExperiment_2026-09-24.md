# Optical Josephson vortex experiment

Run:

```powershell
python tools/optical_transport_overnight.py --config tools/config/optical_junction_vortex_transport.json
```

The experiment uses a 150 nm by 100 nm, 1 nm-grid NbN film with a vertical
resolved weak link centered at 75 nm. A `+1` vortex begins on each bank. The
feedback laser follows only the left-bank vortex and pulls it in `+x` toward
the junction. The right-bank vortex evolves solely through TDGL, magnetic,
thermal, boundary, junction, and vortex-vortex interactions. Add further
objects to `initial_vortices` to study collective motion; change an entry's
charge to `-1` for a vortex-antivortex arrangement.

The live preview shades the junction and marks all detected vortices. Scalar
diagnostics include imposed current, junction voltage, gauge-invariant phase
difference, weak-link amplitude, and voltage normalized by `Phi0*f` when the
microwave drive is active. `junction_dynamics.png` summarizes these signals at
termination.

`electrical_drive.enabled` is false for the first interpretable run. Set it to
true to apply

`I(t) = I_dc + I_ac sin(2*pi*f*t + phase)`

through the left and right contacts. The supplied 50 GHz frequency has a 20 ps
period and is highly resolved by the 0.01 ps timestep. Vary `dc_current_A`
across separate runs while holding the other controls fixed to test for
plateaus near `V/(Phi0*f) = integer`. A proper Shapiro sweep requires multiple
DC biases, enough cycles to discard transients, and comparison against both a
zero-microwave control and the existing SIS RCSJ model.

## Physical interpretation limits

The spatial junction is a phenomenological TDGL weak-link strip formed by
reducing the local GL coefficient. It supports order-parameter suppression,
phase difference, current redistribution, vortex motion, and electrothermal
feedback. It is not a microscopic insulating tunnel barrier and does not
include quasiparticle tunnelling theory, junction capacitance as a localized
Maxwell element, or propagation of an applied electromagnetic wave. The
microwave option is a sinusoidal terminal-current drive. Results can reveal
resolved weak-link phase locking and vortex back-action, but should not be
described as quantitative SIS Shapiro physics without cross-validation.
