# SHS human guide

SHS is a research simulator for thin superconducting films. It evolves the superconducting order parameter together with electrical current, heat, magnetic screening, optical heating, and device geometry. Its strongest current use is controlled two-dimensional TDGL experiments near the critical temperature. It is an experimental physics engine, not yet a general-purpose or experimentally calibrated superconductor design package.

## Where to begin

1. Read [PhysicsGuide.md](PhysicsGuide.md) to understand what each subsystem represents, what behavior it can produce, and where its physical limits are.
2. Read [MathematicsAndSoftwareGuide.md](MathematicsAndSoftwareGuide.md) for equations, normalization, discretization, solver ordering, convergence, and code locations.
3. Use [ConfigurationReference.md](ConfigurationReference.md) when editing JSON inputs.
4. Use [CurrentPhysicsDemonstrations_2026-09-16.md](CurrentPhysicsDemonstrations_2026-09-16.md) for runnable examples.
5. Open `benchmark_results/index.html` after a run to inspect plots, animations, warnings, and raw data.
6. Read [JosephsonJunctions.md](JosephsonJunctions.md) for the resolved weak-link model and I–V/phase sweep workflow.

## What a simulation contains

```text
material + geometry + contacts + operating conditions
                         |
                         v
              mesh and active-domain mask
                         |
                         v
 laser/pinning -> TDGL <-> electrical current -> Joule heat
                    ^             |                 |
                    |             v                 v
                    +------ magnetic field <---- thermal solver
                         |
                         v
             diagnostics and visual output
```

The arrows represent feedback within a physical timestep. The coupled solver repeats the enabled subsystems until their shared fields stop changing within configured tolerances. If convergence fails, that timestep is not a trustworthy physical state.

## Confidence levels

| Capability | Present confidence | Appropriate use |
|---|---|---|
| Geometry, holes, masks, and contacts | Strong software validation | Rectangular thin films with circular or rectangular perforations |
| Gauge-covariant TDGL | Strong numerical qualification; limited experimental calibration | Near-`Tc`, dirty-film, small-timestep studies |
| Current-controlled normal transport | Strong conservation validation | Two-terminal thin-film current redistribution |
| Single-temperature heat flow | Strong conservation validation; simplified material physics | Qualitative hotspot and recovery studies |
| Applied magnetic field | Gauge-equivalence tests available | Uniform perpendicular fields in thin films |
| Thin-film magnetic screening | Early qualified implementation | Weak-to-moderate screening studies with convergence checks |
| Fluxoid and vortex diagnostics | Useful and regression tested | Winding, vortex paths, perforated rings |
| Moving laser heating | Physically sensible exploratory model | Qualitative thermal manipulation and protocol comparisons |
| Reduced vortex dynamics | Parameter-screening model | Rapid force competition and protocol exploration |
| Spatial junctions and SQUIDs | Incomplete | Resolved TDGL weak link and ring groundwork; no complete spatial SQUID yet |
| Lumped SIS junction | Classical RCSJ implemented; microscopic/thermal calibration incomplete | DC/AC phase dynamics, I–V curves, flux modulation |
| SIS vortex/microwave diagnostics | Stationary London phase texture plus AC RCSJ drive | Vortex-modified critical current and Shapiro-step studies |
| Electron–phonon nonequilibrium | Not implemented | Cannot yet predict ultrafast detector or optical thresholds |

## How to judge a result

A plot that looks plausible is not sufficient. Check, in this order:

1. The case completed every requested step and reports coupled convergence.
2. The model-regime section has no unresolved warning about temperature, thickness, or mesh spacing.
3. Terminal current balance and current-continuity residuals are small.
4. Temperature, deposited power, and Joule power remain finite and physically plausible.
5. Vortex winding and fluxoid diagnostics agree where the contour is valid.
6. The behavior persists when the timestep is halved and the mesh is refined.
7. Quantitative claims are compared with an analytic result, published solver, or experiment.

Changing the timestep, mesh, or solver tolerance and obtaining a materially different answer means the original result was not converged.

## Common interpretation errors

- A seeded vortex demonstrates evolution of an imposed initial state. It does not demonstrate spontaneous vortex nucleation.
- A warm Gaussian spot represents absorbed heat. It does not resolve photon absorption, quasiparticle cascades, or separate electron and phonon temperatures.
- A hole with trapped winding is a superconducting ring. It becomes a SQUID only after physically valid Josephson weak links and their phase dynamics are present.
- A calculated self-field is useful only when the screening iteration converges and the thin-film approximation is valid.
- The reduced vortex solver predicts motion from an assumed force law. It does not reproduce core deformation, phase slips, creation, or annihilation.
- Longer runtimes do not improve physical accuracy unless temporal and spatial convergence are maintained.

## Running a documented experiment

```powershell
python tools/tdgl_diagnostics.py --config tools/config/demo_current_driven_vortex.json
python tools/build_results_dashboard.py
```

The diagnostics tool records the effective configuration in its summary. This is the authoritative record of what was actually run, including sweep overrides.
