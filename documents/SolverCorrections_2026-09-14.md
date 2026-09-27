# SHS solver corrections and validation — 2026-09-14

## Baseline

The working tree initially contained only untracked benchmark output under
`benchmark_results/coupled_convergence/`; it was left untouched.

`python -m pytest --collect-only -q` collected **159 tests**. The original full
run was interrupted in TDGL evolution: tests supplied `dt=0.001` as though it
were normalized time, while the solver interpreted seconds. With the NbN time
scale this requested roughly 500 billion explicit substeps per call.

A bounded original-code baseline excluded these entire files:

- `tests/tdgl/test_tdgl_solver.py`
- `tests/tdgl/test_tdgl_thermal_coupling.py`
- `tests/test_thermal_solver.py`
- `tests/test_coupled_solver.py`

That run executed 145 tests: **137 passed, 8 failed** in 195.88 seconds.
Fourteen tests were not executed in that baseline. This is not a claim that the
original complete suite finished. The thermal and legacy coupled files request
up to ten million thermal substeps per call with their original default cap.
The pytest cache also emitted access warnings; verification uses
`-p no:cacheprovider`, without requiring cache-directory permission changes.

| Failing test | Classification and diagnosis |
| --- | --- |
| `test_load_boundary_json` | Stale example expectation: 4.2 K versus JSON's 5 K. |
| `test_builder_boundary` | Same stale boundary expectation. |
| `test_builder_fields` | Stale 3 K initial temperature versus JSON's 9 K. |
| `test_simulation_config` | Stale temperature and current expectations: example now specifies 9 K and 0.01 A. |
| `test_boundary_temperature` | Same stale boundary expectation. |
| `test_monitor_detects_convergence` | Code/policy mismatch: hidden five-iteration minimum contradicted convergence-at-tolerance behavior. Minimum is now explicit, default one. |
| `test_tdgl_suppression_generates_joule_heating` | Test logic: heating need not be strictly positive at every node; assigning E does not prescribe the subsequent voltage solve. |
| `test_supercurrent_is_not_dissipative` | Test logic: the voltage-biased coupled step generates normal current and overwrites manually assigned E/Js. Replaced with a controlled divergence-free current-loop comparison. |

Other inspected test defects: normalized TDGL time was not converted to seconds;
thermal experiments inherited unrelated initial/boundary settings; one uniform
thermal test compared two references to the same mutated field; a 100x100
Python Gauss–Seidel comparison was an unnecessarily expensive correctness test.
Tests now specify their experimental conditions explicitly. The uniform TDGL
and solver-comparison cases use smaller grids without changing their mathematical
purpose. Example-configuration assertions match the shipped example.

## Implementation changes

### Coupling and time

Every candidate starts from the same accepted temperature and order parameter.
Trial endpoint fields supply coupling inputs; they are not repeatedly advanced
as new initial states. Candidate states are private copies, and acceptance
updates the existing Fields container only after convergence. Exceptions and
iteration exhaustion preserve the accepted state. The config=None compatibility
path remains an explicitly sequential electrothermal step without TDGL.

Coupling uses RMS field changes normalized by field-specific characteristic
scales and relative field magnitudes. JSON's optional `coupling` subsection
controls tolerance, iteration budgets, minimum iterations, prediction/checkpoint
settings, and primary-field scales. The default scales are 1 K, 1 for psi,
1 mV, and 1e-7 T m for each vector-potential component. The electrical equation
residual is capped below the outer voltage criterion where necessary; the
user's electrical tolerance is never loosened and iteration budgets are not
silently enlarged. Full-duration runs take a final partial step rather than
rounding the requested duration. Results include accepted elapsed time.

### Electrical transport and units

All coupled entry points share the same configuration-aware electrical adapter,
including the existing phenomenological normal fraction `clip(1-|psi|^2, 0, 1)`.
Direct standalone electrical calls retain their normal-state default for
compatibility; they may pass an explicit superconducting fraction.

The source sign follows `div(-sigma grad(V) + Js)=0`, hence
`div(sigma grad(V))=div(Js)`. Electrical differences, supercurrent divergence,
and normal current use outgoing mesh links. Harmonic face conductivity prevents
normal conduction across a zero-conductivity interface. Unspecified exterior
faces are insulating, not frozen at the initial voltage. Gauss–Seidel and
red-black SOR now solve the same discrete equation. Convergence measures the
maximum diagonal-scaled equation defect in volts, not the iterate change.
Electrical failure raises before any field commit and propagates through the
transactional coupling layer. Disconnected zero-conductivity cells with nonzero
source are rejected.

Shared vector potential is SI (T m); TDGL converts it using TDGLScales. The
normalized supercurrent is converted to A/m^2 before electrical coupling.
Physical time retains SHS's existing tau=pi*hbar/(8*k_B*Tc) convention. This
change does not substitute another TDGL time normalization or add missing
scalar-potential dynamics.

### Thermal, boundaries, and sources

Thermal diffusion is a conservative harmonic-face discretization, with half-width
control volumes at boundary nodes. Fixed temperature, insulation, prescribed
inward heat flux, and convection are implemented. Conflicting fixed-temperature
corners and unsupported or invalid boundary descriptions are rejected. Missing
thermal boundaries on programmatic simulations mean insulation.

Explicit thermal steps respect both the configured cap and the mesh/material/
bath/convection stability bound. TDGL steps similarly respect a bound including
spatial diffusion and the local reaction slope. Requests exceeding one million
internal substeps fail explicitly with a timestep/units message rather than
appearing to hang indefinitely; this is a work-limit guard, not an implicit
change to requested physical time.

TDGL insulating projections use link phases, and right/top Laplacian boundary
links use the actual inward link index. The boundary tests now check the link
derivative rather than a different continuum finite-difference approximation.

Fields separates external and Joule heat sources; heat_source represents their
total. Legacy heat_source input is preserved on the first electrical call.
Subsequent external-source updates should use external_heat_source. Benchmark
thermal adapters have been updated for this ownership contract. The Gaussian
source now respects `(ny, nx)` on rectangular meshes. Contact mapping admits
roundoff at geometric endpoints without expanding contacts by a mesh cell.

NumPy is declared as a runtime dependency; pytest and matplotlib have optional
test/benchmark dependency groups. Benchmark adaptation remains experimental and
separate from production coupling.

## Verification

See the saved pytest logs in `documents/validation/2026-09-14/` for exact output.
The first corrected full run passed 175 tests. Later regression additions cover
all four benchmark heat adapters, the actual shipped configuration, and an
analytical cosine thermal-diffusion case. The final full-suite result is recorded
below after completion.

The added tests exercise: physical-time independence from coupling iteration
count; analytical cooling; failed-step and exception rollback; duration
preservation; SI current conversion; manufactured current continuity/source sign;
insulating electrical edges; equation residual detection under tiny SOR omega;
electrical nonconvergence; external-source persistence; fixed-temperature
boundaries; insulated thermal energy conservation with variable conductivity;
heat-flux/convection energy balance; rectangular hotspot orientation; invalid
physical TDGL timestep workload; benchmark source adapters; coordinated inner/
outer tolerances; and analytical diffusion.

A direct two-step run of the shipped 100x100 NbN configuration with dt=1e-14 s
converged in four total coupling iterations and accepted 2e-14 s. This is a
small integration check, not a physical hotspot benchmark.

## Remaining scientific limitations

Passing these tests does not establish research-grade physical validation.
The model still uses simplified relaxational TDGL and a phenomenological
normal-conductivity suppression law. Generalized gamma dynamics, scalar-potential
TDGL evolution, self-consistent magnetism, transport-contact TDGL boundary types,
multimaterial TDGL, vortex dynamics, and experimental validation remain outside
this correction pass. The perfectly zero-normal-conductivity electrical limit
retains the project's idealized zero-E behavior; it is not a model of a
voltage-biased superconducting terminal.

A corner with nonzero plaquette flux cannot generally satisfy two independently
projected inward-link constraints simultaneously; the y-edge projection takes
precedence at corners. Gauge consistency over arbitrary magnetic configurations
still needs dedicated boundary validation. Joule quantities are represented on
outgoing links; a quantitative coupled electrical/thermal power-balance benchmark
and a fully specified link-to-thermal-volume deposition scheme remain necessary
before claims of full coupled energy conservation. The thermal conservation
tests concern the thermal discretization itself.

The 100x100 mesh on a 5 micrometre film remains too coarse to resolve a 5 nm
coherence length. The original simulation JSON's long duration and very small
internal timestep caps also remain expensive. No default experiment was silently
retuned, and no long benchmark sweep was rerun.

For context on the distinction between insulating and transport-contact TDGL
boundaries and on convention-dependent normalizations, see the primary
[pyTDGL theoretical documentation](https://py-tdgl.readthedocs.io/en/latest/background.html).
Its complete formulation and normalization are not assumed interchangeable with
SHS's current simplified model.

## Final result

`python -m pytest -q -p no:cacheprovider --tb=short`: **181 passed in 53.01 seconds**.
All 159 originally collected tests now execute, with 22 additional regression
cases. No tests were skipped or marked expected-failure in the final run.
