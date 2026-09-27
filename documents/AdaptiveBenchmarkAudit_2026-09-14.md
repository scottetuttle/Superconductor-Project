# Adaptive coupled benchmark audit — 2026-09-14

## Conclusion

`tools/benchmark_adaptive_coupled_updated.py` is now a substantially safer experimental coupling harness, but the underlying simulator is not yet a complete electrothermal TDGL model. The corrected benchmark can compare solver ordering, relaxation, tolerances, substep limits, cost, and retry behavior without advancing a state multiple times or leaking a failed attempt. It cannot yet determine a universally “best” strategy from physical error because it has no local truncation-error estimator, no electromagnetic solve, and no voltage/scalar-potential term in TDGL.

The recommended near-term role is an **offline policy laboratory**: run controlled scenario matrices, record accuracy and cost, and promote only well-validated controller rules into `CoupledSolver`. It should not yet be treated as an autonomous production controller.

## Corrections made during this audit

- Full-state checkpoints now include derived current and heat fields. Failed attempts restore the accepted state completely.
- TDGL and thermal differential solves now start from the accepted time level while consuming the current endpoint guesses. This makes coupling iterations solve one physical endpoint instead of repeatedly advancing time.
- All six operation orderings now exchange current guesses consistently.
- The electrical adapter now honors configured SOR omega and iteration limits.
- Adapted TDGL models preserve configured `u`, `gamma`, `kappa`, and stability safety factor.
- Coupling residuals use mesh-independent RMS values and configured physical scales.
- An unrecoverable time window stops the run and rolls back instead of continuing later requested steps.
- The benchmark can set its inner electrical tolerance from the outer coupling target, avoiding accidental oversolving.
- The heat source is centered in physical coordinates and its Gaussian radius is specified in meters.
- Added safeguarded Aitken relaxation, complete solver-work counters, direct adaptation outcome accounting, unique timing plots, and focused regression tests.

## What is working

### Time and rollback semantics

Every coupling attempt represents the same interval from accepted state \(y_n\) to candidate endpoint \(y_{n+1}\). Only a converged candidate advances physical time. Reducing `physical_dt` subdivides the requested outer interval, and failed attempts do not advance time.

### Thermal implementation used by the benchmark

The benchmark calls the production conservative thermal solver. That implementation supports insulating, fixed-temperature, fixed-flux, and convection boundaries; separates external and Joule heating; applies mesh/material stability limits; and evolves each coupling candidate from the accepted temperature.

### Electrical continuity implementation

The benchmark calls the common production electrical adapter. It solves the discrete continuity equation with harmonic face conductivity, contact Dirichlet values, insulating exterior faces, and supercurrent divergence as a source. The benchmark now records total SOR work, which is essential because coupling iteration count alone concealed most runtime.

### Ordering consistency after correction

A 100×100, one-step test at `dt=1e-14 s`, initial `sigma=0.5`, and benchmark tolerance `1e-3` converged for all six orderings. With Aitken acceleration, the cases required 3–4 coupling iterations. Mean current ranged from approximately `4.29943e8` to `4.29952e8 A/m²`; final residuals were about `6.19e-5`. Exact results are in `documents/validation/2026-09-14/adaptive_ordering_aitken_smoke.json`.

A strict 9×7 test with coupling tolerance `1e-8` also converged for all orderings in 3–4 iterations. Exact results are in `documents/validation/2026-09-14/adaptive_strict_small_grid.json`.

### Aitken acceleration

For the representative 100×100 ordering `tdgl_electrical_thermal`, fixed `sigma=0.5` previously required 12 coupling iterations and about 5.5 seconds. Safeguarded Aitken relaxation required 3 iterations and about 3.6–3.9 seconds. Electrical work fell from 15,626 to 10,242 SOR iterations. The wall-time improvement is smaller than the outer-iteration improvement because the initial electrical solves dominate cost.

## What is not working as a complete physical model

### Voltage does not enter TDGL dynamics directly

The implemented TDGL equation is simplified relaxational TDGL. It omits scalar/electrochemical potential dynamics, generalized-gamma dynamics, and transport-contact TDGL boundary conditions. Electrical bias affects `psi` mainly through Joule heating; it does not presently produce the complete phase evolution, current-driven suppression, or phase-slip physics expected from an electrostatic TDGL formulation. Although `gamma` is preserved in configuration, the evolution equation does not use it.

### Electromagnetics is dormant

No operation updates vector potential or magnetic field. `Ax` and `Ay` residuals are therefore always zero. There is no self-field, screening-current, or displacement-current evolution. Ordering results do not test electromagnetic coupling.

### Conductivity is phenomenological

Normal conductivity is multiplied by `clip(1-|psi|², 0, 1)`. This is a useful development closure, not a validated quasiparticle transport model. The perfectly superconducting limit creates an underdetermined voltage field and is handled by forcing electric field to zero. Benchmark success under this closure does not validate superconducting transport physics.

### The configured “current” is not a current boundary

The simulation configuration contains `current.value`, but the electrical solve is voltage-driven. Contacts named `left_current` and `right_current` receive imposed voltages. A future adaptive backbone must distinguish current-controlled and voltage-controlled experiments and monitor terminal current/power consistently.

### Contact placement may not describe full-film transport

The film width is 50 nm, while `right_current.x` is 4.8 nm and `left_voltage.x` is 2.4 nm. All contacts occupy the left part of the film. If these coordinates were intended to scale from the earlier 5 µm geometry, likely values would be 48 nm and 24 nm. This audit preserved the user-authored geometry, but current benchmark results should not be interpreted as transport across the full 50 nm film until the intended layout is confirmed.

### Energy and charge validation are incomplete

The electrical equation has manufactured-solution tests and the thermal discretization has conservation tests, but the fully coupled benchmark does not yet report:

- terminal electrical power versus integrated Joule power;
- current-continuity defect at every accepted step;
- deposited optical energy versus thermal energy change and boundary loss;
- free-energy behavior in applicable unbiased TDGL cases;
- gauge-invariant equivalence under gauge transformations.

These should become acceptance gates, rather than diagnostic plots only.

## Controller limitations

### Coupling convergence is not time-integration accuracy

The controller observes fixed-point endpoint changes. Reducing `thermal_max_substep` or `tdgl_max_normalized_timestep` changes discretization error, but a smaller coupling residual does not prove that the time-discretized physical solution is accurate. The benchmark needs step doubling or an embedded method to estimate local temporal error independently.

### Timestep adaptation only shrinks

Failed attempts can reduce physical `dt`, but successful easy steps never grow it. Long simulations will retain an unnecessarily small step after a transient. A future controller should use accepted-step error history with safety factors, rejection factors, growth limits, and event caps.

### One global tolerance obscures weak signals

The current normalized maximum residual is useful for initial experiments, but a global tolerance can accept an ordering lag that is tiny relative to the 9 K absolute temperature while still being large relative to a microkelvin temperature increment. Production-quality acceptance needs per-field absolute and relative tolerances and separate physical invariant limits.

### Adaptation choices are heuristic

Tightening an electrical tolerance helps only when algebraic error limits coupling. Reducing an explicit integrator cap helps temporal error, not necessarily nonlinear coupling. Reducing physical `dt` cannot repair a poor relaxation factor in an algebraic fixed point. The controller should classify failure into algebraic solve error, nonlinear coupling error, temporal error, physical invariant violation, and invalid configuration before choosing an action.

### No higher-order fixed-point accelerator

Aitken relaxation is now available and materially helps the tested case. Strong nonlinear regimes should also evaluate Anderson/IQN-ILS acceleration with bounded updates, history filtering, and fallback to damped iteration.

## Performance findings

On the current 100×100 case, a single initial electrical solve measured:

| Electrical tolerance | SOR iterations | Wall time |
| --- | ---: | ---: |
| `1e-6 V` | 48 | 0.018 s |
| `1e-8 V` | 6,095 | 2.28 s |
| `1e-10 V` | 23,809 | 8.83 s |

Repeated red-black SOR dominates the coupled benchmark. The benchmark’s coupling-scaled inner-tolerance policy avoids solving many orders more accurately than the outer experiment requires. The longer-term solution is a sparse Krylov or multigrid electrical backend with reusable topology/preconditioning, equation-defect verification, and warm starts. Plotting and history formatting are not the principal bottleneck.

The 100×100 grid resolves the 5 nm coherence length with roughly ten intervals after the geometry changed to 50 nm. This is much better than the former 5 µm layout, but convergence under further mesh refinement is still untested.

## Recommended adaptive backbone

1. **Establish acceptance gates.** Check finite values, current continuity, power balance, thermal energy balance, physical bounds, and gauge-invariant observables after every candidate step.
2. **Separate error layers.** Track inner algebraic defect, outer nonlinear coupling defect, temporal truncation error, and physical invariant error independently.
3. **Use weighted per-field tolerances.** Configure absolute and relative tolerances for temperature, `psi`, voltage, current, and heat; do not infer all accuracy from one scalar scale.
4. **Add step-doubling validation.** Compare one step of `dt` with two steps of `dt/2`; use the weighted error to accept/reject and drive a PI or DSP timestep controller.
5. **Retain Aitken as the low-cost accelerator.** Add Anderson/IQN-ILS as an experimental option with history filtering, physical bounds, and automatic fallback.
6. **Make inner solves inexact but coordinated.** Loosen early electrical solves according to the current outer residual and tighten them near convergence, while always verifying the final equation defect.
7. **Allow safe timestep recovery.** Grow `dt` after several low-error accepted steps, with maximum growth, safety, and event-specific caps.
8. **Optimize the dominant kernel.** Replace or supplement Python SOR with sparse Krylov/multigrid methods and record setup, iteration, and solve cost separately.
9. **Build a scenario matrix.** Include equilibrium relaxation, pure diffusion, pure Joule heating, hotspot onset/recovery, near-`Tc` transitions, normal-state intervals, strong bias, mesh refinement, and eventual vortex/field cases.
10. **Select policies by accuracy-adjusted cost.** Compare runtime and work only among runs that pass the same physical error gates. Store the configuration, code revision, environment, and deterministic seed with every result.

## Guidance from related solver projects

- [preCICE implicit coupling documentation](https://precice.org/configuration-coupling) treats multiphysics coupling as a fixed-point problem with explicit convergence measures and checkpoint/repeat semantics.
- [preCICE acceleration documentation](https://precice.org/configuration-acceleration) provides constant relaxation, Aitken, and IQN/Anderson-family acceleration, and recommends participant solves materially tighter than the outer coupling target.
- [PETSc SNES](https://petsc.org/main/manual/snes/) exposes nonlinear Richardson, line searches, nonlinear GMRES, Anderson mixing, quasi-Newton, and trust-region methods behind a common monitored interface. This supports treating coupling strategy as a selectable solver rather than a chain of special cases.
- [PETSc TSAdapt](https://petsc.org/main/manualpages/TS/TSAdaptSetSafety/) separates timestep error/stability goals from nonlinear convergence and uses distinct accepted/rejected-step safety factors.
- [SUNDIALS ARKODE](https://sundials.readthedocs.io/en/v7.2.1/arkode/Introduction_link.html) supports adaptive explicit, implicit, IMEX, and multirate integration. Its architecture is a useful model for separating fast TDGL/electrical work from slower thermal evolution.
- [MOOSE executioners and timesteppers](https://mooseframework.inl.gov/releases/moose/v1.0.0/syntax/Executioner/index.html) combine Picard multiphysics iterations with selectable timestep controllers, including iteration-count-based adaptation.
- [pyTDGL numerical background](https://py-tdgl.readthedocs.io/en/latest/background.html) uses adaptive TDGL steps based on recent `|psi|²` dynamics, retries invalid updates with reduced steps, and cautions that this heuristic needs a strict maximum or disabling in the fully normal state.

These projects point to the same architecture: explicit state checkpoints, independent error measures, interchangeable nonlinear acceleration, guarded accept/reject decisions, and timestep control based on integration error or relevant dynamics rather than coupling iteration count alone.

## Validation

`python -m pytest -q -p no:cacheprovider`: **189 passed in 113.79 seconds**. This includes eight focused adaptive-benchmark tests for checkpoint completeness, accepted-time thermal evolution, TDGL parameter preservation, electrical configuration propagation, mesh-independent residuals, physical Gaussian placement, invalid source width, and Aitken relaxation.
