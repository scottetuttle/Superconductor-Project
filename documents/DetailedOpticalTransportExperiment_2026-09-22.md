# Detailed 5 nm optical transport experiment

## Higher-field core-localization correction (2026-09-23)

The `optical_transport_fast_follow_subcritical_6` run stopped at step 607 with
three detected `+1` windings. The transported central winding was still in its
preceding plaquette; the other two were more than 60 nm away. The failure came
from the bilinear complex-zero interpolation becoming singular in the broad,
heated central core and returning `NaN`, rather than from identity selection.

Core localization now uses the bilinear complex zero when it is
well-conditioned, a quadratic local minimum of `|psi|^2` when it is not, and
the center of the already validated winding plaquette as a bounded final
fallback. The configured maximum displacement continues to reject identity
swaps. This allows tracking through multiple-vortex states while preserving
the winding detector as the topological source of truth.

## Vortex identity tracking

The transport runner supports a continuation tracker through the optional
`tracking` configuration block. With `require_unique_topology: false`, it
enumerates all detected `+1` windings and follows the complex zero nearest the
previous accepted vortex position. `maximum_jump_m` rejects an implausible
one-step displacement, so a newly detected edge vortex cannot steal the
identity of the central transported vortex. Global positive and negative
counts are still written to diagnostics, together with `topology_is_unique`;
extra vortices are recorded rather than hidden.

This addresses the failure in `optical_transport_fast_follow_subcritical_2`.
At its failure frame the original central vortex remained present while a
second `+1` winding appeared near the lower boundary. The former unique-count
rule stopped the run even though the transported vortex had not been lost.
That boundary winding may represent physical entry or a low-amplitude boundary
detection event and should still be examined. If no nearby `+1` winding exists,
the runner stops with `tracked_vortex_unresolved` and saves a field frame. New
field frames include both vector-potential components so gauge-invariant
winding can be reconstructed exactly after a run.

Run `python tools/optical_transport_overnight.py --config tools/config/optical_transport_detailed_large.json`. This uses the 150 nm uniform reference film while adaptive two-level evolution remains unqualified. The vortex begins near the center and the target is 5 nm in +y. The feedback beam remains nominally 15 nm ahead and is speed limited. Pinning is disabled so the first run isolates optical, current and boundary effects.

Laser hysteresis reads the interpolated electron temperature at the tracked vortex, not the maximum temperature anywhere in the film. It switches off at 15.15 K and can restart below 14.4 K after a 50-step dwell. The NbN reference Tc is 15.5 K. There is no global-temperature hard stop, so the remote beam center may cross Tc. The run can still terminate if the unique vortex topology is lost, which is necessary because a normal region passing through the vortex makes its position undefined. This protocol permits a larger surrounding gradient without deliberately heating the tracked core through Tc.

`diagnostics.csv` is flushed during execution and records position, target distance, sampled velocity, electron and phonon temperature at the vortex, global peak temperature, target-directed thermal and GL gradients, laser state, beam position, local current-derived Lorentz estimate, London thermal-force estimate, and nearest-boundary distance. The London estimate uses the reference penetration depth, coherence length and Tc; it is not a complete calibrated force balance. Boundary/image forces, vortex viscosity, pinning and nonequilibrium quasiparticle forces are not separately reconstructed. Graphs `transport_diagnostics.png` and `motion_and_force_proxies.png` are generated at termination.

Compressed `field_frames/frame_*.npz` files are written atomically every 50 steps and at start/termination. Each contains electron and phonon temperature, complex psi, vortex and beam position, laser state, step and time. They survive a later crash independently of the 200-step restart checkpoint. After a completed or interrupted run, execute `python tools/visualize_optical_transport.py <result-directory>` to generate the trajectory GIF, position and thermal kymographs, checkpoint plots, and `transport_fields.gif` showing temperature, amplitude and phase. With 10,000 steps, the default cadence produces about 201 frames and may consume roughly 100–150 MB depending on compression.

The run stops on target arrival within 0.25 nm, topology loss, solver failure, three consecutive 1,000-step windows with less than 0.05 nm targetward progress each, 10,000 steps, or ten wall-clock hours. It is restartable from the most recent committed checkpoint with the exact same configuration.

## Demonstrated 10 nm run and fast-follow control

`optical_transport_detailed_large_3` reached 9.768 nm targetward displacement in 11.24 ps and stopped 0.241 nm from the 10 nm target. Its mean finite-window velocity was about 869 m/s and its peak sampled velocity was 1,183 m/s. Velocity and target-directed temperature gradient had correlation 0.952 over finite velocity samples. The gradient rose to 0.101 K/nm, while the local vortex temperature reached 15.146 K. The global hotspot reached 15.540 K, as allowed by local-core thermal control. The laser was off for only two saved samples near the end.

The beam moved at an average 280 m/s, close to its configured 300 m/s cap, while the vortex became much faster. Consequently, beam–vortex y separation collapsed from about 15 nm to 8.35 nm. The London thermal-force estimate was about 24 times the largest local Lorentz estimate in this zero-applied-current experiment. These observations identify feedback bandwidth and the thermal drive as the first acceleration variables; they do not establish a calibrated physical NbN terminal velocity.

`tools/config/optical_transport_fast_follow.json` changes only the beam speed cap from 300 to 2,000 m/s and the update interval from five steps to one. Power, spot width, thermal thresholds and 10 nm target remain identical. This is the correct first comparison because it tests whether preserving the requested lead sustains the gradient without confounding the result with more heating. New runs also record sampled acceleration. If fast-follow improves the late velocity while maintaining topology, the next controlled sweep should vary absorbed power and offset separately.
