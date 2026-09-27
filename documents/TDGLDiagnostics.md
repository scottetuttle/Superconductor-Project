# TDGL configuration diagnostics

Run the diagnostic sweep from the repository root:

```powershell
python tools/tdgl_diagnostics.py --config tools/config/tdgl_diagnostics.json
```

Build the human-readable results browser after any diagnostic runs:

```powershell
python tools/build_results_dashboard.py
```

Open `benchmark_results/index.html` in a browser. It discovers every top-level run summary, displays final physical and solver metrics, highlights model-regime warnings, embeds all generated PNG/GIF media, and links the underlying CSV and JSON data. Its SQUID-readiness panel distinguishes currently validated ring physics from missing junction-level SQUID physics.

The JSON file documents every diagnostic control inline. Its `sweeps` object accepts lists under dotted `SimulationConfig` attribute names and forms their Cartesian product. For example, testing two `gamma` values and three voltages creates six independent cases:

```json
"sweeps": {
  "tdgl.gamma": [0.0, 2.0],
  "electrical.voltage_left": [1e-7, 1e-6, 1e-5]
}
```

Each case starts from a newly built simulation. This prevents the final state of one configuration from contaminating the next configuration.

## Outputs

Each case directory contains:

- `diagnostics.csv`: physical time histories, conservation defects, vortex counts, temperatures, order-parameter magnitudes, electrical residuals, and internal solver work.
- `summary.json`: settings, convergence state, runtime, total coupling work, and final diagnostics.
- `timeseries.png`: order parameter, temperature, signed vortex counts, and current continuity.
- `kymographs.png`: centerline temperature, amplitude, and phase versus time.
- `final_state.png`: final amplitude, phase, temperature, voltage, and detected vortex locations.
- `current_and_field.png`: total, superconducting, and normal-current maps with current streamlines, plus applied and current-generated magnetic fields.
- `order_parameter.gif`: animated amplitude, absolute phase, and spatial phase change when GIF output is enabled.
- `vortex_trajectories.csv`: persistent signed vortex identities and positions in cells and metres.
- `vortex_trajectories.png`: signed vortex tracks and the commanded laser path over final `|psi|`.

The absolute-phase panel always uses the cyclic `[-pi, pi]` range so colors are comparable between frames. Small voltage-driven changes can therefore look stationary on that panel. The spatial-phase-change panel subtracts the uniform global phase rotation, which is gauge dependent, and dynamically scales the remaining change relative to the first frame. Its frame title reports the actual maximum change in radians. Phase is masked wherever `|psi|` is effectively zero because it is undefined there.

The parent directory contains `comparison.png` and a combined `summary.json`. Existing runs are preserved by adding a numeric suffix to the requested output directory.

Each saved row also evaluates a rectangular fluxoid contour configured by `fluxoid.contour_inset_fraction`. It reports magnetic-flux and supercurrent contributions, the physical London fluxoid in units of `Phi0`, the topological phase winding, and their difference. A `NaN` result means the contour crosses a point below the configured amplitude floor, where phase and the London current term are not reliable.

## Vortex interpretation

The detector integrates the gauge-covariant phase around each plaquette and restores the vector-potential circulation to obtain signed topological winding. A positive count identifies counterclockwise `+2*pi` winding and a negative count identifies clockwise winding. Plaquettes touching points below `vortex_detection.amplitude_floor` are ignored because phase is undefined where the order parameter is zero.

A seeded vortex is an initial-condition probe. Its survival, motion, or disappearance can compare algorithms, but it is not proof that the present model generated a physical vortex. Quantitative vortex dynamics still require self-consistent magnetic-vector-potential evolution, realistic contacts, and calibrated material parameters.

## Choosing sweeps

- Sweep `tdgl.gamma` to measure how amplitude/phase relaxation affects vortex motion and solver cost.
- Sweep `run.dt` in separate configuration files to check time-step convergence. Results should approach a stable trajectory as `dt` decreases.
- Sweep `tdgl.max_normalized_timestep` to distinguish internal-integration error from outer coupling behavior.
- Sweep `electrical.voltage_left` to locate transport-driven suppression or vortex-entry thresholds.
- Sweep `coupling.tolerance` and `electrical.solver.tolerance` to confirm that physical observables stop changing before spending more solver work.
- Sweep hotspot temperature and radius by creating separate diagnostic JSON files; these values describe initial conditions rather than `SimulationConfig` attributes.

Compare physical diagnostics first. A faster case is useful only when its current-continuity defect, temperature, energy, order parameter, and vortex trajectory remain within the required accuracy.

## Current and magnetic-field modes

The supplied transport configuration uses full left and right edge contacts. In `current` drive mode, the production finite-volume electrical solver applies `target_current_A` directly as balanced Neumann flux at the source and sink terminals. Its augmented sparse system fixes only the arbitrary mean-voltage gauge and retains every local charge-continuity equation. `voltage` mode uses fixed terminal potentials. `voltage_search` retains the older rollback-based voltage controller for controlled comparisons, but is no longer the preferred current implementation.

A uniform perpendicular field is represented through a symmetric or Landau-gauge vector potential. The order parameter is phase-transformed when selecting a Landau gauge so the gauges describe the same initial physical state. Every saved sample reports the error between the numerical curl of the vector potential and the requested magnetic field.

The current-generated field uses the thin-film Biot-Savart law at a configured positive observation height. It is a post-processing diagnostic. It does not yet feed the induced vector potential back into TDGL, so cases with large self-fields or strongly nonuniform currents require a future self-consistent screening solve.
