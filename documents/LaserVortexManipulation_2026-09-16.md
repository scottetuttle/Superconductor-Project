# Moving-laser vortex manipulation checkpoint

## Physical purpose

Focused optical heating can lower the local superconducting condensation energy and create a transient attractive region for an Abrikosov vortex. Experiments have demonstrated manipulation of individual flux quanta by local laser heating, while theoretical work emphasizes competition among thermal attraction, background pinning, beam size, power, and translation speed.

Primary references:

- I. S. Veshchunov et al., *Optical Manipulation of Single Flux Quanta*, arXiv:1608.07453, https://arxiv.org/abs/1608.07453
- *Optical vortex manipulation for topological quantum computation*, arXiv:2104.15012, https://arxiv.org/abs/2104.15012

## Implemented model

Simulation configuration now accepts a `laser` section with `enabled`, `absorbed_power_W`, Gaussian standard deviation `sigma_m`, `repeat_path`, and piecewise-linear waypoints containing `time_s`, `x_m`, and `y_m`.

The beam uses the areal absorbed-power profile

`I(r) = P_abs/(2*pi*sigma^2) exp(-|r-r_laser|^2/(2*sigma^2))`.

It is converted to volumetric power through local film thickness and added to external and Joule heating. Power falling outside the active film or into a hole is not reassigned to the remaining material. The coupled solver evaluates the beam at the midpoint of each physical timestep. Coordinate grids and repeated coupling-iteration profiles are cached.

Pinning emerges through TDGL temperature dependence: local heating suppresses the equilibrium condensate amplitude and reduces the energy cost of placing a vortex core there. No empirical force is applied directly to a vortex.

## Vortex paths

Diagnostics match same-sign vortices between saved frames using bounded nearest-neighbor assignment. They save `vortex_trajectories.csv` and plot `vortex_trajectories.png`, overlaying signed tracks and the commanded laser path on final order-parameter amplitude. `vortex_tracking.max_displacement_cells` controls the largest allowed motion between saved frames.

## Validation

Automated tests verify exact path interpolation, Gaussian deposited power, loss of power falling into an insulating hole, zero deposition in inactive geometry, thermal energy gain equal to deposited optical energy, and persistent sign-preserving vortex identities.

The reference run uses a `500 nW`, `3 nm` spot moving from the right side of the ring toward its upper side over `150 fs`. It converged, deposited `497.4 nW` at the final sampled position, preserved the seeded `+1` vortex, and generated a trajectory plot. The vortex remained near its initial location and shifted only one grid cell. This is a useful negative result: the reference route is too fast and/or too weak to drag the vortex under the current one-temperature dynamics.

## Limits

The model treats absorbed optical energy as instantaneous effective-temperature heating. It does not yet include wavelength-dependent absorption, reflection, optical penetration depth, separate electron and phonon temperatures, electron-phonon transfer delay, nonequilibrium quasiparticles, static material pinning, stochastic thermal activation, or calibrated experimental optical efficiency. Quantitative capture thresholds and dragging velocities require those additions and experimental comparison.

## Feasibility program

Laser dragging must remain a tested hypothesis rather than an assumed feature. A useful result may be a map showing that no stable capture-and-transport window exists for a particular material and geometry. Development should proceed in this order:

1. Sweep absorbed power, spot size, path speed, transport current, field, and initial laser-vortex separation. Record capture distance, vortex lag, displacement, release, annihilation, peak temperature, and topological-charge conservation.
2. Add a configurable static pinning-energy landscape. **Implemented:** Gaussian defects now reduce the local linear GL coefficient, with additive overlap and a configured cap. This permits depinning competition tests, although strengths still require calibration.
3. Replace instantaneous one-temperature heating with electron and phonon temperatures, electron-phonon energy transfer, substrate escape, and a measured or wavelength-dependent absorption model.
4. Add fluctuation terms consistent with the chosen TDGL normalization and timestep so thermally activated depinning can be studied statistically across repeated seeds.
5. Calibrate material parameters and optical coupling against experiments before claiming quantitative thresholds. Compare stationary-beam capture, release after beam removal, maximum translation velocity, and failure modes.

The first-stage model is still useful for rejecting implausible settings and exposing numerical or topological failures. It is not sufficient by itself to establish that a predicted dragging protocol will work in a device.

Diagnostics now record the nearest laser-vortex separation at every saved frame and create `laser_vortex_distance.png`. Summaries deliberately leave `capture_claimed` false: a future sweep must define a distance threshold, require sustained following for a minimum duration, and verify controlled release before labeling a run successful.
