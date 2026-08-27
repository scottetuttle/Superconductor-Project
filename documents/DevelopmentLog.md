# SHS Development Log

**Project:** Superconducting Hotspot Simulator (SHS)

This document records significant development milestones and scientific decisions.

It exists to preserve project continuity independently of conversation history.

---

# August 10, 2026

## TDGL Foundation Milestone

A major development milestone was completed in the transition from classical electrothermal transport toward superconducting physics.

The first functional dimensionless TDGL framework was implemented and validated.

### Implemented

The following TDGL infrastructure was developed or substantially expanded:

* `tdgl/parameters.py`
* `tdgl/gl.py`
* `tdgl/model.py`
* `tdgl/operators.py`
* `tdgl/initialization.py`
* `solvers/tdgl_solver.py`
* `tdgl/__init__.py`

The implementation introduced or substantially developed:

* dimensionless TDGL parameters
* complex superconducting order parameter
* equilibrium superconducting states
* temperature-dependent equilibrium amplitude
* equilibrium order-parameter phase
* ordinary spatial derivatives
* gauge-covariant gradient
* gauge-covariant Laplacian
* TDGL temporal evolution

---

## Equilibrium Validation

The homogeneous equilibrium relationship

[
|\psi|_{\mathrm{eq}}
====================

\sqrt{1-\frac{T}{T_c}}
]

for (T<T_c), with zero order parameter above (T_c), was implemented and tested.

Validation included:

* zero temperature
* intermediate reduced temperature
* near-(T_c) temperature
* exactly (T_c)
* above (T_c)

The equilibrium order parameter was also tested with nonzero phase.

---

## Differential Operator Validation

The TDGL operator layer was tested against simple analytical cases.

Validated behavior included:

* zero gradient for constant fields
* unit derivative for a linear field
* zero Laplacian for a constant field
* expected covariant gradient under constant vector potential
* expected covariant Laplacian under constant vector potential

This established the initial gauge-aware spatial differential infrastructure.

---

## TDGL Dynamic Validation

The TDGL solver was tested using a spatially uniform system.

The numerical solution was compared against the analytical solution of:

[
u\frac{d\psi}{dt}
=================

a\psi-\psi^3.
]

The numerical transient reproduced the analytical solution within the tested tolerance.

Additional validation established:

* stable low-temperature evolution
* decay of an initially superconducting state above (T_c)
* convergence toward the expected superconducting equilibrium

---

## Thermal-to-TDGL Coupling

A localized thermal hotspot was introduced into an initially superconducting NbN simulation.

The hotspot temperature was set to approximately:

[
0.95T_c.
]

TDGL evolution produced a lower local order-parameter amplitude inside the hotspot than in the surrounding colder region.

This established the first demonstrated spatial coupling:

[
T(\mathbf r)\rightarrow\psi(\mathbf r).
]

The current validation establishes qualitative suppression of superconductivity.

It does not yet establish the complete electrothermal feedback loop.

---

## Validation Result

The associated TDGL tests were passing at the completion of the milestone.

The Git commit associated with this milestone modified:

* 14 files
* approximately 1,339 lines added
* approximately 318 lines removed

This represents the first major transition of SHS from superconductivity infrastructure planning into an actively functioning TDGL implementation.

---

## Scientific Significance

SHS can now demonstrate the fundamental temperature-dependent behavior expected from the dimensionless TDGL model:

[
T<T_c
\Rightarrow
|\psi|>0
]

[
T\rightarrow T_c
\Rightarrow
|\psi|\rightarrow0
]

[
T>T_c
\Rightarrow
\psi\rightarrow0.
]

It can also demonstrate local suppression of superconductivity caused by a thermal perturbation.

This provides the foundation required for future current, magnetic, vortex, and electrothermal superconducting physics.

---

## Known Limitations After This Milestone

The following remain incomplete:

* full gauge-invariant discretization
* link variables
* gauge-transformation validation
* complete TDGL boundary conditions
* supercurrent calculation
* current-driven dynamics
* electromagnetic self-consistency
* magnetic screening
* vortex dynamics
* full electrothermal feedback
* optical excitation
* experimental validation

---

## Next Development Target

The next major target is a scientifically robust gauge-covariant TDGL implementation.

The immediate sequence is:

1. Review current covariant discretization.
2. Introduce link variables where appropriate.
3. Validate gauge transformations.
4. Establish TDGL boundary conditions.
5. Develop supercurrent calculations.
6. Validate phase and current dynamics.

---

# Documentation Rule Going Forward

Significant development milestones should be recorded here during development rather than waiting until the end of the day.

If development is interrupted unexpectedly, the most recent completed milestone should still provide enough information to reconstruct the project state.

## August 13, 2026 — TDGL Supercurrent

### Added

Implemented the analytical TDGL supercurrent density in:

- `src/shs/tdgl/model.py`

The model now provides the normalized supercurrent density derived from the complex order parameter and gauge-covariant gradient.

For the current normalized formulation:

`j_s = Im(psi* D psi)`

where:

`D = ∇ - iA`

and `psi*` is the complex conjugate of the superconducting order parameter.

The implementation uses the existing `covariant_gradient()` operator rather than introducing a separate numerical derivative.

### Tests

Added analytical supercurrent tests to:

- `tests/test_tdgl_model.py`

Tests verify:

- zero current for a uniform order parameter with zero vector potential
- expected current from a phase gradient
- expected current response to a vector potential
- correct current behavior for a complex order parameter

All TDGL model tests pass.

### Result

SHS can now calculate the local superconducting current directly from the TDGL order parameter.

This establishes the first explicit connection between:

`psi → phase/amplitude → supercurrent`

and provides the foundation for future current redistribution, electromagnetic coupling, and vortex physics.


added TDGL supercurrent calculation and integrated superconducting current fields into simulation state.


Replaced the TDGL covariant Laplacian with a gauge-consistent link-variable discretization and updated operator tests to validate the discrete formulation. Full test suite passes.

2026-08-27 — Electrical Driving, TDGL/Electrical Coupling, and Test Cleanup

Completed:

Investigated the electrical solver's imposed-voltage behavior after identifying that the solver had previously been using an inappropriate/default voltage value.
Changed the electrical solver's left-contact voltage for the current test case to:
voltage_left = 1e-3 V
voltage_right = 0.0 V
Re-ran the fully coupled thermal/TDGL/electrical test.
Confirmed that the electrical solution now produces a physically plausible response to the imposed voltage:
left/right voltage difference remains approximately 1 mV
electric field is finite and stable
normal current develops as superconductivity is suppressed
Joule heating increases as the normal current develops
|psi| decreases during evolution
fields remain finite throughout the tested evolution.
The coupled test therefore appears to be functionally behaving correctly at this development stage, although this test does not yet constitute rigorous physical validation of the complete coupled model.
Reviewed electrical_solver.py and physics/electrical.py with particular attention to:
imposed voltage boundary conditions
normal versus superconducting current
conductivity suppression by superconducting fraction
electric-field calculation
Joule heating.
Confirmed that the electrical model consistently separates:
J = J_s + J_n
J_n = sigma_eff E
Q_J = J_n · E
Identified a broader architectural issue for future cleanup: physical/control values that are currently embedded directly in solvers, tests, or Python modules should eventually be moved into the configuration system. Voltage is the first clear example.
Decided that solver/model code should generally not contain large amounts of diagnostic print() output. Diagnostics should eventually be represented through structured fields/results or dedicated debugging/analysis tools instead.
Converted the previously diagnostic-heavy coupled heating test into the direction of an automated test:
verify finite fields
verify superconducting suppression
verify imposed voltage drop
verify dissipative heating.
The existing test suite remains useful, but the current coupled tests should be regarded primarily as sanity/integration tests, rather than complete physical validation tests.

Important result:

The electrical/TDGL coupling is now behaving coherently enough to move forward. The remaining uncertainty is primarily how rigorously the numerical and physical behavior is being validated, rather than an obvious failure of the electrical coupling itself.

Next development step:

Perform a systematic audit of the entire test suite. For each test, determine whether it is:

a unit/component test,
a numerical solver test,
an integration/coupling test,
a physical validation test, or
primarily a diagnostic test that should eventually be replaced or removed.

Also begin identifying all hard-coded physical/control parameters in the solver stack and migrate them toward the configuration system where appropriate.