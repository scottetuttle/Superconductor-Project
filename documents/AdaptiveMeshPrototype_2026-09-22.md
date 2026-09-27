# Isolated adaptive-mesh feasibility prototype

The production uniform-grid solver is unchanged. Experimental code lives in `src/shs/experimental/adaptive_mesh.py`, and `tools/nested_mesh_feasibility.py` reads a saved checkpoint without evolving or modifying it. Run `python tools/nested_mesh_feasibility.py --config tools/config/nested_mesh_feasibility.json`.

## What the existing solver already reuses

Every physical `coupled_step` copies the previously accepted fields and starts its fixed-point trial from that state. TDGL and thermal integration also start from the accepted `psi`, electron temperature, and phonon temperature. Heun supplies a predictor inside one TDGL substep. The sparse electrical backend receives the previous voltage, but a direct solve does not gain iterations from that initial value. The current solver does not extrapolate from several accepted physical time levels.

The experimental history predictor uses three accepted states. It estimates the local first and second temporal differences and damps linear extrapolation when the recent change is curved or reverses. Quiet cells naturally receive nearly zero correction because their latest difference is nearly zero. Complex `psi` is predicted through amplitude and wrapped temporal phase increments rather than raw phase. This remains experimental: a changing scalar-potential gauge, vortex crossing, topology change, or abrupt laser switch must reject or strongly damp a prediction. A predictor can reduce nonlinear or coupled iterations, but cannot avoid evaluating explicit TDGL/thermal updates on every active cell.

## First coarse-exterior result

The checkpoint at step 2200 of `optical_transport_overnight_6` was represented with a factor-four coarse exterior, an exact 15 nm-radius fine patch around the vortex, and a 4 nm transition band. The estimated composite node count was 1,737 instead of 10,201, a nominal 83.0% reduction before interface overhead.

The reconstructed exterior had 0.00374 K RMS electron-temperature error and 0.000969 RMS order-parameter-amplitude error. Maximum errors were 0.0338 K and 0.00252 respectively. These scalar/amplitude results are promising for a fixed nested-patch evolution test. The gauge-covariant x-link phase RMS error was 0.0938 rad, while the y-link error was 0.000511 rad. The large x-link error means ordinary interpolation of real and imaginary `psi` is not acceptable as the production coarse/fine transfer operator. It likely reflects the strong global phase variation and normal contacts along x.

## Required next qualification

Implement gauge-covariant restriction and prolongation by parallel-transporting complex `psi` to a common reference with the link variables before averaging or interpolation. Preserve winding and flux explicitly. Then evolve a fixed fine patch coupled to a coarse exterior and compare it at identical physical times with a uniform-grid reference. Qualification must include vortex position and count, covariant current, temperature and energy balance, current continuity, phase-link error, and behavior at the coarse/fine interface. A moving patch should only be attempted after the fixed interface passes those checks.

## Centered larger-film qualification

The follow-up reference uses `configs/simulations/nbn_optical_vortex_large_unpinned.json`: a 150 nm square film with the original 1 nm reference spacing. The vortex begins near (75.5, 75.5) nm and a stationary 250 nW, 10 nm-sigma beam is 15 nm ahead in +y. `tools/prepare_nested_mesh_reference.py` generated `benchmark_results/nested_mesh_simple_reference_2` after 20 dark and 30 illuminated steps. This keeps the vortex near the center and ends at a peak electron temperature of 13.0583 K. An inherited inconsistency was found in the source simulation, whose film and bath were 13 K while fixed thermal contacts were 14 K; the isolated large-film configuration sets its contacts to 13 K. The source configuration was not changed.

The experimental transfer now includes coarse restriction by injection of `psi` and composition of fine link holonomies into coarse vector-potential links. Prolongation transports every coarse-corner `psi` to the target node through link variables before bilinear interpolation. A random-gauge test confirms that prolongation commutes with a gauge transformation to numerical precision. This is separate from the production solver.

The simple feasibility case uses a factor-five exterior, a 20 nm exact fine radius, and a 5 nm transition. It estimates 2,858 composite nodes versus 22,801 uniform nodes, an 87.5% nominal reduction. Exterior electron-temperature RMS error is 6.81e-5 K and amplitude RMS error is 0.00278. Ordinary complex interpolation gives x/y covariant link-phase RMS errors of 0.00301/0.000544 rad. Parallel transport reduces these to 0.000586/0.000528 rad: about a five-fold improvement in x and a small improvement in y. This is a much clearer result than the late overnight checkpoint, though the amplitude error and interface conservation still require qualification during actual evolution.

The next stage remains a fixed two-level evolution. Coarse and fine updates must exchange heat and current flux conservatively (“refluxing”), and the fine patch must receive boundary data through the gauge-covariant prolongation. The current result validates representation and gauge transformation behavior; it does not yet establish stable time evolution or a speedup.
