# Current physics demonstrations — 2026-09-16

These configurations emphasize relative behavior and conservation checks supported by the current solver. Seeded vortices are initial conditions; they must not be presented as predicted nucleation events.

| Demonstration | Configuration | Primary observation | Strongest interpretation |
|---|---|---|---|
| Vortex dipole relaxation | `tools/config/demo_vortex_dipole_relaxation.json` | Opposite winding cores attract, deform, and may annihilate | TDGL phase and core dynamics without external forcing |
| Current-driven vortex | `tools/config/demo_current_driven_vortex.json` | Reversing current reverses transverse vortex motion | Direction and symmetry of current-induced motion |
| Fluxoid field bias | `tools/config/demo_fluxoid_field_bias.json` | Ring currents redistribute while winding and fluxoid are monitored | Multiply connected, gauge-covariant response; not a SQUID |
| Laser capture/release | `tools/config/demo_laser_capture_release.json` | Laser-vortex separation with zero and finite transport current | Qualitative thermal pinning and competition with transport force |

Run any demonstration from the repository root:

```powershell
python tools/tdgl_diagnostics.py --config tools/config/demo_vortex_dipole_relaxation.json
```

Replace the filename for another case, then rebuild the combined browser with:

```powershell
python tools/build_results_dashboard.py
```

For scientific comparison, first inspect `summary.json` for convergence and model-regime warnings. Then compare vortex trajectories, current/field maps, kymographs, and the diagnostic CSV. Repeat important conclusions at half the timestep and finer mesh spacing. Absolute vortex speed, optical power threshold, nucleation threshold, and escape time are not yet experimentally calibrated.

## Initial end-to-end check

The laser capture/release comparison completed all 16 steps in both cases and passed the configured model-regime checks. The terminal-current solver recovered 5.000000000013 nA in the driven case with a relative continuity residual near `1.2e-14`. The final laser-to-nearest-vortex separation was about 7.08 nm in both cases. At this short 160 fs duration, 5 nA did not produce a resolved difference in the tracked vortex position; that is a useful null result rather than evidence of strong current manipulation. Longer-time integration is needed before increasing the force simply to make the plot move.
