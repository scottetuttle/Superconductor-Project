# SHS TDGL Test Suite Context

Purpose:
This file documents the tests contained in tests/tdgl/. It is intended to preserve the exact role, assumptions, controls, limitations, and scientific meaning of the current TDGL test suite so development can continue without re-auditing these files.

The TDGL tests cover:
- TDGL boundary-condition infrastructure
- Gauge-covariant operators
- Covariant Laplacian
- TDGL model equations
- Supercurrent density
- Order-parameter initialization
- TDGL scaling and normalization
- TDGL solver behavior
- Thermal coupling of the order parameter

Important distinction:
These tests primarily establish mathematical correctness, numerical behavior, and limiting-case physics. They are not yet a complete validation of the physical TDGL model against experimental or published results.


============================================================
## tests/tdgl/test_tdgl_boundary.py
============================================================

Purpose:
Tests the TDGL boundary-condition data structures and the implementation of insulating boundaries.

Tests:
- test_tdgl_boundary_creation
- test_tdgl_boundary_set
- test_insulating_boundary_zero_vector_potential
- test_insulating_boundary_preserves_shape
- test_insulating_boundary_enforces_covariant_normal_derivative


### test_tdgl_boundary_creation

Creates a TDGLBoundaryCondition using:
- side = LEFT
- type = INSULATING

Proves:
- TDGL boundary conditions can be constructed.
- The supplied side and boundary type are stored correctly.

Dependencies:
- TDGLBoundaryCondition
- TDGLBoundarySide
- TDGLBoundaryType

Hard-coded controls:
- LEFT boundary — TEST CONDITION
- INSULATING boundary type — TEST CONDITION

Classification:
TDGL software/data-structure unit test.

Scientific significance:
Low directly. This verifies the boundary-condition infrastructure rather than physical behavior.

Assessment:
Keep. Simple structural test.


### test_tdgl_boundary_set

Creates a TDGLBoundarySet and adds an insulating left boundary.

Proves:
- A TDGL boundary condition can be added to the boundary set.
- The boundary set recognizes the specified boundary side.

Dependencies:
- TDGLBoundarySet
- TDGLBoundaryCondition

Hard-coded controls:
- LEFT side — TEST CONDITION
- INSULATING type — TEST CONDITION

Classification:
TDGL boundary infrastructure test.

Assessment:
Keep.


### test_insulating_boundary_zero_vector_potential

Creates a uniform complex order parameter:
- psi = 1
- Ax = 0
- Ay = 0

Applies the insulating boundary operator.

Proves:
- Applying the insulating boundary operation to a uniform state with zero vector potential does not alter that state.

Classification:
TDGL boundary limiting-case test.

Scientific significance:
Moderate as a numerical invariant.

Assessment:
Keep. Important simple limiting case.


### test_insulating_boundary_preserves_shape

Uses a random complex 20 x 20 order parameter and zero vector potential.

Proves:
- Applying the boundary operation preserves the shape of psi.

Classification:
Numerical/data integrity test.

Scientific significance:
Low physically, but useful for preventing accidental array-shape corruption.

Assessment:
Keep.


### test_insulating_boundary_enforces_covariant_normal_derivative

Uses:
- Random complex psi
- Ax = 0.1
- Ay = 0.2
- dx = 1
- dy = 1

Applies the insulating boundary condition and explicitly evaluates the covariant normal derivative at:
- Left boundary
- Right boundary
- Bottom boundary
- Top boundary

The test expects each normal covariant derivative to be zero.

Proves:
- The insulating boundary implementation enforces the intended covariant Neumann condition on all four sides.

The tested form is approximately:

D_n psi = 0

where the ordinary derivative is corrected by the vector potential.

Classification:
TDGL physics/numerical boundary-condition test.

Scientific significance:
High relative to the current TDGL infrastructure. This is one of the more meaningful tests in this file because it checks the actual mathematical boundary condition rather than only data structure behavior.

Important limitation:
The test verifies the implementation against the same finite-difference form used in the test. It does not independently validate that this discretization is the optimal or physically complete representation of the continuum insulating boundary condition.

Assessment:
Strong test. Keep.


Overall assessment of test_tdgl_boundary.py:

This file establishes both the boundary-condition infrastructure and the basic physical requirement for insulating boundaries.

The strongest test is test_insulating_boundary_enforces_covariant_normal_derivative.

Future improvements:
- Test each boundary independently rather than relying on one combined test.
- Add tests for nonzero and spatially varying vector potentials.
- Eventually verify consistency between boundary treatment and the gauge-covariant operator discretization.


============================================================
## tests/tdgl/test_tdgl_diagnostics.py
============================================================

Purpose:
Tests diagnostic calculations derived from the complex superconducting order parameter, current density, covariant gradients, and free energy.

Tests:
- test_order_parameter_amplitude
- test_order_parameter_amplitude_squared
- test_order_parameter_phase
- test_supercurrent_magnitude
- test_covariant_gradient_zero_for_uniform_zero_field
- test_covariant_gradient_detects_vector_potential
- test_free_energy_density_equilibrium
- test_free_energy_density_normal_state
- test_total_free_energy
- test_amplitude_statistics


### test_order_parameter_amplitude

Uses:
- psi = 3 + 4i
- psi = 1

Expected amplitudes:
- 5
- 1

Proves:
- The order-parameter amplitude is calculated as |psi|.

Classification:
Mathematical unit test.

Scientific significance:
Low individually, but fundamental to essentially all TDGL diagnostics.

Assessment:
Keep.


### test_order_parameter_amplitude_squared

Uses:
- psi = 3 + 4i
- psi = 1

Expected:
- |psi|² = 25
- |psi|² = 1

Proves:
- The amplitude-squared diagnostic is implemented correctly.

Classification:
Mathematical unit test.

Assessment:
Keep.


### test_order_parameter_phase

Uses:
- psi = 1
- psi = i

Expected phases:
- 0
- pi/2

Proves:
- The phase diagnostic correctly extracts the complex phase.

Classification:
Mathematical unit test.

Scientific significance:
Moderate because phase is physically important for superconducting current and gauge behavior.

Assessment:
Keep.


### test_supercurrent_magnitude

Uses:
- Jx = [[3,0],[0,1]]
- Jy = [[4,0],[0,0]]

Expected magnitude:
- [[5,0],[0,1]]

Proves:
- Current magnitude is calculated as sqrt(Jx² + Jy²).

Classification:
Mathematical diagnostic test.

Assessment:
Keep.


### test_covariant_gradient_zero_for_uniform_zero_field

Uses:
- Uniform psi = 1
- Ax = 0
- Ay = 0

Proves:
- The covariant gradient vanishes for a spatially uniform order parameter in zero vector potential.

Classification:
TDGL limiting-case physics test.

Scientific significance:
High. This is an important invariant.

Assessment:
Strong test. Keep.


### test_covariant_gradient_detects_vector_potential

Uses:
- Uniform psi = 1
- Ax = 1
- Ay = 0

Proves:
- A nonzero vector potential produces a nonzero covariant derivative even when psi itself is spatially uniform.

For the tested continuous form:
D_x psi = -i A_x psi

The resulting magnitude is expected to be 1.

Classification:
Gauge-covariant TDGL physics test.

Scientific significance:
High. It establishes that the vector potential actually enters the covariant derivative.

Assessment:
Strong test. Keep.


### test_free_energy_density_equilibrium

Uses:
- Reduced temperature = 0.5
- |psi| = sqrt(0.5)

The expected free-energy density is -0.125 for the implemented normalized expression.

Proves:
- The free-energy expression produces the expected equilibrium value for a superconducting state at reduced temperature 0.5.

Classification:
TDGL analytical physics unit test.

Scientific significance:
Moderate-to-high.

Important limitation:
This validates the implemented normalized free-energy expression, but does not independently establish that the chosen normalization corresponds to the intended physical GL free energy.

Assessment:
Keep.


### test_free_energy_density_normal_state

Uses:
- psi = 0
- reduced temperature = 1.5

Expected:
- free-energy density = 0

Proves:
- The normal state produces the expected value under the implemented free-energy convention.

Classification:
TDGL limiting-case test.

Assessment:
Keep.


### test_total_free_energy

Uses:
- 4 x 4 grid
- psi = 1
- reduced temperature = 0

Expected total free energy:
- 16 * (-1 + 0.5)
- -8

Proves:
- Total free energy is correctly summed from the local free-energy density.

Classification:
Numerical/physics diagnostic test.

Important limitation:
The test does not include gradient or magnetic contributions, so it validates only the portion represented by the current free-energy implementation.

Assessment:
Keep.


### test_amplitude_statistics

Uses a 2 x 2 array with amplitudes:
- 1
- 2
- 3
- 4

Expected:
- Mean = 2.5
- Minimum = 1
- Maximum = 4

Proves:
- Basic order-parameter amplitude statistics are calculated correctly.

Classification:
Diagnostic utility unit test.

Scientific significance:
Low physically.

Assessment:
Keep.


Overall assessment of test_tdgl_diagnostics.py:

This is a clean collection of inexpensive diagnostic tests.

The most scientifically important tests are:
- covariant gradient with zero field
- covariant gradient with nonzero vector potential
- free-energy equilibrium

Future improvements:
- Add known phase-gradient cases.
- Add gauge-invariant diagnostic tests.
- Add free-energy tests including gradient energy.
- Add tests confirming that physically equivalent gauge choices produce equivalent gauge-invariant diagnostics.


============================================================
## tests/tdgl/test_tdgl_equilibrium.py
============================================================

Purpose:
Tests initialization of uniform and equilibrium superconducting order-parameter states.

Tests:
- test_uniform_superconducting_state
- test_uniform_superconducting_state_with_amplitude
- test_equilibrium_state_at_zero_temperature
- test_equilibrium_state_at_half_tc
- test_equilibrium_state_near_tc
- test_equilibrium_state_above_tc


### test_uniform_superconducting_state

Creates a 10 x 10 uniform state.

Expected:
- Shape = 10 x 10
- dtype = complex
- psi = 1 + 0i everywhere

Proves:
- Uniform superconducting initialization works as intended.

Classification:
Initialization unit test.

Assessment:
Keep.


### test_uniform_superconducting_state_with_amplitude

Creates a uniform state with amplitude 0.5.

Proves:
- The initialization function accepts a requested amplitude.
- The resulting |psi| is 0.5 everywhere.

Classification:
Initialization unit test.

Assessment:
Keep.


### test_equilibrium_state_at_zero_temperature

Uses reduced temperature = 0.

Expected:
- |psi| = 1

Proves:
- The equilibrium initialization reproduces the normalized zero-temperature superconducting amplitude.

Classification:
TDGL analytical limiting-case test.

Assessment:
Keep.


### test_equilibrium_state_at_half_tc

Uses reduced temperature = 0.5.

Expected:
- |psi| = sqrt(0.5)

Proves:
- The equilibrium initialization follows the normalized GL relation:

|psi|² = 1 - T/Tc

for the tested superconducting regime.

Classification:
Analytical TDGL physics test.

Assessment:
Strong test.


### test_equilibrium_state_near_tc

Uses reduced temperature = 0.95.

Expected:
- |psi| = sqrt(0.05)

Proves:
- The equilibrium order parameter approaches zero as Tc is approached from below.

Classification:
TDGL limiting-case physics test.

Assessment:
Strong test.


### test_equilibrium_state_above_tc

Uses reduced temperature = 1.1.

Expected:
- psi = 0

Proves:
- The equilibrium initialization correctly produces the normal state above Tc.

Classification:
TDGL phase-transition limiting-case test.

Assessment:
Strong test.


Overall assessment:

This file provides good coverage of the normalized equilibrium amplitude relationship and the superconducting/normal transition.

It should remain as a foundational initialization test.


============================================================
## tests/tdgl/test_tdgl_gauge.py
============================================================

Purpose:
Tests gauge-link construction and gauge covariance of the covariant derivative.

Tests:
- test_gauge_links_have_unit_magnitude
- test_zero_vector_potential_reduces_to_gradient
- test_gauge_transformation_preserves_covariant_gradient


### test_gauge_links_have_unit_magnitude

Generates random Ax and Ay arrays and constructs gauge links using:
- spacing = 0.5

Proves:
- Gauge links have unit magnitude.

This is expected because the links are phase factors of the form exp(i theta).

Classification:
Gauge numerical unit test.

Scientific significance:
Moderate. Unit magnitude is an important structural property of a pure phase link variable.

Assessment:
Keep.


### test_zero_vector_potential_reduces_to_gradient

Uses:
- psi = X + iY
- Ax = 0
- Ay = 0
- dx = dy = 1

Expected:
- Dx psi = 1
- Dy psi = i
away from the final finite-difference boundary.

Proves:
- The gauge-covariant gradient reduces to the ordinary gradient when the vector potential is zero.

Classification:
TDGL operator limiting-case test.

Scientific significance:
High.

Assessment:
Strong test. Keep.


### test_gauge_transformation_preserves_covariant_gradient

Constructs:
- A spatially varying complex psi
- Constant Ax = 0.3
- Constant Ay = -0.2
- A nonlinear gauge function chi(x,y)

The transformed fields use:
- psi' = psi exp(i chi)
- A' = A + grad(chi)

The test then verifies:

D'psi' = exp(i chi) Dpsi

for the interior regions.

Proves:
- The implemented covariant derivative transforms correctly under the tested gauge transformation.

Classification:
Advanced TDGL gauge-invariance/covariance test.

Scientific significance:
Very high relative to the current TDGL infrastructure.

Important limitation:
The gauge-transformed vector potential is calculated using forward finite differences while the covariant derivative uses its own discretization. This means the test validates the particular discrete implementation rather than an independent continuum gauge-invariance proof.

Assessment:
One of the strongest TDGL tests currently present. Keep and protect carefully.


============================================================
## tests/tdgl/test_tdgl_initialization.py
============================================================

Purpose:
Tests the equilibrium superconducting-state initialization function.

Tests:
- test_equilibrium_initialization_at_half_tc
- test_equilibrium_initialization_above_tc
- test_equilibrium_initialization_preserves_phase


### test_equilibrium_initialization_at_half_tc

Uses:
- Shape = 10 x 10
- Reduced temperature = 0.5

Expected:
- |psi| = sqrt(0.5)

Proves:
- Equilibrium initialization follows the expected normalized temperature dependence.

Classification:
TDGL analytical initialization test.

Assessment:
Keep.


### test_equilibrium_initialization_above_tc

Uses:
- Reduced temperature = 1.1

Expected:
- psi = 0

Proves:
- Initialization correctly enters the normal state above Tc.

Classification:
Phase-transition limiting-case test.

Assessment:
Keep.


### test_equilibrium_initialization_preserves_phase

Uses:
- Reduced temperature = 0.5
- Phase = pi/2

Expected:
- amplitude = sqrt(0.5)
- real part approximately zero

Proves:
- The requested global phase is preserved while the equilibrium amplitude is determined independently.

Classification:
Complex-order-parameter initialization test.

Scientific significance:
Moderate. Phase is fundamental to the superconducting order parameter.

Assessment:
Keep.


Overall assessment:

This file overlaps somewhat with test_tdgl_equilibrium.py.

The duplication is not harmful because equilibrium.py tests the public initialization behavior while this file specifically emphasizes the initialization implementation and phase preservation.

Future cleanup could consolidate redundant tests if the suite becomes unnecessarily large.


============================================================
## tests/tdgl/test_tdgl_laplacian.py
============================================================

Purpose:
Tests the ordinary and gauge-covariant Laplacian.

Tests:
- test_covariant_laplacian_constant_state
- test_covariant_laplacian_uniform_vector_potential
- test_laplacian_constant_state_including_boundaries


### test_covariant_laplacian_constant_state

Uses:
- psi = 1
- Ax = 0
- Ay = 0
- dx = dy = 1

Proves:
- The covariant Laplacian of a constant state in zero vector potential is zero.

Classification:
TDGL numerical limiting-case test.

Scientific significance:
High.

Assessment:
Strong invariant. Keep.


### test_covariant_laplacian_uniform_vector_potential

Uses:
- psi = 1
- Ax = 1
- Ay = 0
- dx = dy = 1

Expected interior result:

2 cos(1) - 2

This corresponds to the discrete gauge-covariant second derivative for the selected uniform vector potential.

Proves:
- The covariant Laplacian responds correctly to a uniform vector potential under the implemented gauge-link discretization.

Classification:
Gauge-covariant numerical operator test.

Scientific significance:
High.

Assessment:
Keep.


### test_laplacian_constant_state_including_boundaries

Uses:
- psi = 1
- dx = dy = 1

Proves:
- The ordinary Laplacian returns zero for a constant field, including at the boundaries.

Classification:
Numerical operator invariant test.

Assessment:
Strong test. Keep.


Overall assessment:

This file is important because the Laplacian is a central component of the TDGL evolution equation.

The uniform-vector-potential case is particularly useful because it checks behavior beyond the trivial zero-field case.

Future improvements:
- Add known quadratic-function analytical tests.
- Add convergence tests with decreasing grid spacing.
- Add additional gauge-transformation tests for the Laplacian.


============================================================
## tests/tdgl/test_tdgl_model.py
============================================================

Purpose:
Tests the TDGL model parameters, equilibrium relations, order parameter construction, and superconducting current density.

Tests:
- test_tdgl_parameters_are_valid
- test_tdgl_equilibrium_amplitude_zero_temperature
- test_tdgl_equilibrium_amplitude_half_tc
- test_tdgl_equilibrium_amplitude_near_tc
- test_tdgl_equilibrium_amplitude_at_tc
- test_tdgl_equilibrium_amplitude_above_tc
- test_tdgl_equilibrium_amplitude_squared
- test_tdgl_equilibrium_order_parameter
- test_tdgl_equilibrium_order_parameter_phase
- test_supercurrent_uniform_zero_field
- test_supercurrent_from_vector_potential
- test_supercurrent_from_phase_gradient
- test_supercurrent_small_vector_potential_approaches_continuum


### Parameter and equilibrium tests

The equilibrium amplitude tests verify:

| Reduced temperature | Expected |psi| |
| 0 | 1 |
| 0.5 | sqrt(0.5) |
| 0.95 | sqrt(0.05) |
| 1.0 | 0 |
| 1.2 | 0 |

Proves:
- The TDGL model implements the normalized equilibrium relation.
- The superconducting order parameter disappears at and above Tc.
- The squared equilibrium amplitude follows max(1 - T/Tc, 0).

Classification:
Analytical TDGL model tests.

Scientific significance:
High for the normalized model.

Assessment:
Strong set of tests. Keep.


### test_tdgl_equilibrium_order_parameter

Uses:
- Reduced temperature = 0.75
- Phase = 0

Expected:
- |psi| = sqrt(0.25)
- imaginary part = 0

Proves:
- The model combines equilibrium amplitude and phase correctly.

Classification:
TDGL model unit test.

Assessment:
Keep.


### test_tdgl_equilibrium_order_parameter_phase

Uses:
- Reduced temperature = 0
- Phase = pi/2

Expected:
- real part = 0
- imaginary part = 1

Proves:
- The model correctly applies phase to the equilibrium order parameter.

Classification:
Complex-order-parameter unit test.

Assessment:
Keep.


### test_supercurrent_uniform_zero_field

Uses:
- psi = 1
- Ax = 0
- Ay = 0

Proves:
- A uniform order parameter in zero vector potential carries zero supercurrent.

Classification:
Superconducting physics limiting-case test.

Scientific significance:
High.

Assessment:
Strong invariant. Keep.


### test_supercurrent_from_vector_potential

Uses:
- psi = 1
- Ax = 1
- Ay = 0
- dx = dy = 1

Expected:
- Jx = -sin(1)
- Jy = 0

Proves:
- The implemented supercurrent responds to vector potential through the gauge-covariant phase relation.

Classification:
Gauge-coupled superconducting-current test.

Scientific significance:
High.

Assessment:
Strong test.


### test_supercurrent_from_phase_gradient

Uses:
- psi = exp(i k x)
- k = 0.1
- Ax = Ay = 0

Expected:
- Jx approximately k
- Jy approximately 0

Proves:
- A phase gradient generates supercurrent.

Classification:
Core superconducting physics test.

Scientific significance:
Very high.

Assessment:
Strong test. Keep.


### test_supercurrent_small_vector_potential_approaches_continuum

Uses:
- psi = 1
- Ax = 0.01
- Ay = 0
- dx = dy = 1

Expected:
- Jx approximately -0.01

This uses the small-angle approximation:

sin(A) approximately A

Proves:
- The discrete current expression approaches the expected continuum result for small vector potential.

Classification:
Continuum-limit numerical test.

Scientific significance:
High.

Assessment:
Strong test.


Overall assessment of test_tdgl_model.py:

This is one of the most important TDGL unit-test files.

It verifies:
- normalized equilibrium behavior
- complex order parameter construction
- phase dependence
- zero-current equilibrium
- current induced by vector potential
- current induced by phase gradients
- continuum small-field behavior

The supercurrent tests are particularly valuable because they establish the basic gauge-coupled superconducting transport behavior before the full solver is validated.

Future improvements:
- Test current under simultaneous phase gradient and vector potential.
- Test gauge covariance of the supercurrent itself.
- Test current scaling with |psi|².
- Eventually compare against known GL current-density relations using physical units.


============================================================
## tests/tdgl/test_tdgl_operators.py
============================================================

Purpose:
Tests the basic ordinary and gauge-covariant differential operators used by TDGL.

Tests:
- test_tdgl_gradient_constant
- test_tdgl_gradient_linear
- test_covariant_gradient_zero_field
- test_covariant_gradient_phase_field
- test_tdgl_laplacian_preserves_shape
- test_covariant_laplacian_uniform_zero_field
- test_covariant_laplacian_uniform_vector_potential
- test_covariant_laplacian_gauge_transformation


### test_tdgl_gradient_constant

Uses:
- Uniform complex psi
- dx = dy = 1

Proves:
- Ordinary gradient of a constant field is zero.

Classification:
Numerical operator invariant.

Assessment:
Keep.


### test_tdgl_gradient_linear

Uses:
- psi = x
- dx = dy = 1

Expected:
- dpsi/dx = 1
- dpsi/dy = 0

Proves:
- The ordinary gradient correctly differentiates a linear field.

Classification:
Analytical numerical-operator test.

Assessment:
Strong test.


### test_covariant_gradient_zero_field

Uses:
- psi = 1
- Ax = Ay = 0

Proves:
- Covariant gradient reduces to zero for the constant zero-field case.

Classification:
Gauge-covariant limiting case.

Assessment:
Keep.


### test_covariant_gradient_phase_field

Uses:
- psi = 1
- Ax = 1
- Ay = 0

Expected:
- Dx psi has imaginary component -1
- Dy psi = 0

Proves:
- The vector potential enters the covariant gradient with the correct sign and imaginary structure for the tested case.

Classification:
Gauge-covariant operator test.

Scientific significance:
High.

Assessment:
Keep.


### test_tdgl_laplacian_preserves_shape

Uses:
- 100 x 100 psi

Proves:
- The Laplacian preserves field shape.

Classification:
Numerical/data integrity test.

Assessment:
Keep.


### test_covariant_laplacian_uniform_zero_field

Uses:
- psi = 1
- Ax = Ay = 0

Proves:
- Covariant Laplacian is zero for a constant zero-field state.

Classification:
TDGL limiting-case test.

Assessment:
Keep.


### test_covariant_laplacian_uniform_vector_potential

Uses:
- psi = 1
- Ax = 1
- Ay = 0
- dx = dy = 1

Expected:

2 cos(1) - 2

Proves:
- The gauge-covariant Laplacian produces the expected discrete result for a uniform vector potential.

Classification:
Gauge-covariant numerical operator test.

Scientific significance:
High.

Assessment:
Strong test.


### test_covariant_laplacian_gauge_transformation

Uses:
- psi = 1
- A = 0
- chi = kx
- k = 0.3
- dx = dy = 0.2

Transforms:
- psi' = exp(i chi)
- A' = grad(chi)

Then checks:

D²psi' = exp(i chi) D²psi

in the interior.

Proves:
- The covariant Laplacian is gauge covariant under the tested transformation.

Classification:
Advanced gauge-covariance test.

Scientific significance:
Very high.

Important limitation:
As with the gradient gauge test, the gauge transformation and discrete operator use related finite-difference assumptions. This is strong evidence for internal consistency, but it is not independent physical validation.

Assessment:
One of the strongest operator tests in the suite. Keep.


Overall assessment:

test_tdgl_operators.py overlaps with test_tdgl_gauge.py and test_tdgl_laplacian.py.

This is not necessarily a problem. The files emphasize different levels:
- operators.py = general operator behavior
- gauge.py = explicit gauge transformations
- laplacian.py = focused Laplacian behavior

Future organization could consolidate overlapping tests, but this should not be done until the TDGL infrastructure is stable.


============================================================
## tests/tdgl/test_tdgl_scaling.py
============================================================

Purpose:
Tests the physical-to-dimensionless scaling system used by TDGL.

Representative NbN parameters:
- coherence length xi = 5e-9 m
- penetration depth lambda = 2e-7 m
- critical temperature Tc = 15.5 K
- time scale = 1e-12 s

Tests:
- validation of positive physical scales
- length conversion
- temperature conversion
- time conversion
- vector-potential scale
- magnetic-field scale
- current-density scale
- reversibility of conversions
- array shape preservation


### Validation tests

The scaling system rejects:
- xi <= 0
- lambda <= 0
- Tc <= 0
- time scale <= 0

Proves:
- Invalid physical scaling parameters are rejected.

Classification:
Software validation / physics-parameter validation.

Assessment:
Keep.


### test_length_conversion

Uses:
- physical length = 5 nm
- xi = 5 nm

Expected:
- dimensionless length = 1

Proves:
- xi is correctly used as the length scale.

Classification:
Physical normalization unit test.

Assessment:
Keep.


### test_length_conversion_is_reversible

Tests several lengths:
- 1 nm
- 5 nm
- 100 nm
- 5 um

Converts physical -> dimensionless -> physical.

Proves:
- Length conversion is reversible.

Classification:
Scaling numerical correctness test.

Assessment:
Strong test.


### test_temperature_conversion

Uses:
- T = 15.5 K
- Tc = 15.5 K

Expected:
- reduced temperature = 1

Proves:
- Temperature normalization is correctly defined by Tc.

Classification:
Scaling unit test.

Assessment:
Keep.


### test_temperature_conversion_is_reversible

Tests:
- 0 K
- 3 K
- 7.75 K
- 15.5 K
- 20 K

Proves:
- Temperature conversion is reversible across below-Tc, at-Tc, and above-Tc values.

Classification:
Scaling numerical correctness.

Assessment:
Strong test.


### test_time_conversion_is_reversible

Tests:
- 1e-15 s
- 1e-12 s
- 1e-9 s
- 1e-6 s

Proves:
- Time scaling is reversible across several orders of magnitude.

Classification:
Scaling unit test.

Assessment:
Keep.


### test_vector_potential_scale

Checks:

A_scale = Phi0 / (2 pi xi)

Proves:
- Vector-potential normalization uses the expected flux-quantum/coherence-length relationship.

Classification:
Analytical physical-scaling test.

Scientific significance:
High.

Assessment:
Strong test.


### test_magnetic_field_scale

Checks:

B_scale = Phi0 / (2 pi xi²)

Proves:
- Magnetic-field normalization is consistent with the selected TDGL length scale.

Classification:
Analytical physical-scaling test.

Assessment:
Strong test.


### test_current_density_scale

Checks:

J_scale = Phi0 /
          (2 pi mu0 lambda² xi)

Proves:
- Current-density scaling follows the selected GL normalization.

Classification:
Analytical physical-scaling test.

Scientific significance:
High.

Important limitation:
The test verifies the formula implemented by the scaling class, but the physical normalization itself should eventually be checked against the exact nondimensionalization used by the TDGL equations.

Assessment:
Keep.


### test_vector_potential_conversion_is_reversible

Tests multiple physical vector-potential values.

Proves:
- Vector-potential conversion is reversible.

Classification:
Scaling numerical correctness.

Assessment:
Keep.


### test_magnetic_field_conversion_is_reversible

Tests multiple magnetic-field values.

Proves:
- Magnetic-field conversion is reversible.

Classification:
Scaling numerical correctness.

Assessment:
Keep.


### test_current_density_conversion_is_reversible

Tests:
- 0
- 1e6
- 1e8
- 1e10 A/m²

Proves:
- Current-density conversion is reversible over several orders of magnitude.

Classification:
Scaling numerical correctness.

Assessment:
Keep.


### test_array_conversion_preserves_shape

Uses a 20 x 20 array.

Proves:
- Scaling operations work elementwise on arrays without changing their shape.

Classification:
Numerical/data integrity test.

Assessment:
Keep.


Overall assessment of test_tdgl_scaling.py:

This is an important infrastructure test file because incorrect normalization would contaminate every subsequent physical TDGL result.

The strongest aspects are:
- analytical scale formulas
- reversible conversions
- multi-order-of-magnitude testing
- invalid-parameter checks

Major future requirement:
The scaling formulas themselves need eventual validation against the exact nondimensionalization adopted by the project. Passing these tests only proves internal consistency with the implemented formulas.


============================================================
## tests/tdgl/test_tdgl_solver.py
============================================================

Purpose:
Integration-level tests of the TDGL time-step solver using the standard NbN simulation.

Shared setup:
- Loads configs/geometry/NbN_film.json
- Creates the mesh
- Loads NbN material
- Builds region map
- Builds material map
- Creates Fields
- Creates insulating TDGL boundaries on all four sides
- Creates a Simulation
- Uses TDGLModel and tdgl_step

Standard physical material:
- NbN
- Tc = 15.5 K
- coherence length = 5 nm
- penetration depth = 200 nm
- thickness = 100 nm

Tests:
- test_tdgl_low_temperature_stability
- test_tdgl_uniform_state_matches_analytical_solution
- test_tdgl_normal_state_above_tc
- test_tdgl_equilibrium_above_tc
- test_tdgl_uniform_state_approaches_equilibrium
- test_tdgl_solver_enforces_insulating_boundaries


### test_tdgl_low_temperature_stability

Initializes:
- Temperature = 3 K
- psi = 1
- dt = 0.001

Performs one TDGL step.

Proves:
- The resulting order parameter remains finite.
- The amplitude remains nonnegative.
- The change in amplitude is small.

Required change:
- maximum amplitude change < 1e-2

Classification:
TDGL solver stability test.

Scientific significance:
Moderate.

Important limitation:
This is primarily a sanity/stability check. It does not establish convergence, accuracy, or correct physical time evolution.

Assessment:
Keep.


### test_tdgl_uniform_state_matches_analytical_solution

Initializes:
- Temperature = 3 K
- psi = 1
- dt = 0.001
- 1000 steps
- total dimensionless time = 1

Uses the analytical solution of the spatially uniform TDGL equation:

u dpsi/dt = a psi - psi^3

where:

a = 1 - T/Tc

The numerical final amplitude is compared against the analytical solution with:
- absolute tolerance = 5e-4

Also verifies that the solution decreases from the initial amplitude.

Proves:
- The TDGL time integrator reproduces the analytical transient for the spatially uniform nonlinear TDGL equation.
- The relaxation direction is correct.

Classification:
TDGL solver analytical-validation test.

Scientific significance:
Very high.

This is currently one of the strongest solver tests in the suite.

Assessment:
Strong. Keep and protect.

Future improvement:
Repeat for several dt values to establish temporal convergence and determine the observed order of accuracy.


### test_tdgl_normal_state_above_tc

Uses:
- Temperature = 1.1 Tc
- psi = 1
- dt = 0.001
- 2000 steps

Proves:
- Above Tc, the initially superconducting order parameter decreases.

Classification:
TDGL phase-transition solver behavior test.

Scientific significance:
High qualitatively.

Important limitation:
The test only verifies decay, not the quantitative decay rate.

Assessment:
Keep.


### test_tdgl_equilibrium_above_tc

Directly checks the model equilibrium amplitude at reduced temperature 1.1.

Expected:
- 0

Proves:
- The model identifies the normal state above Tc.

Classification:
Analytical model test.

Important note:
This overlaps with test_tdgl_model.py.

Assessment:
Redundant but harmless. Could eventually be consolidated.


### test_tdgl_uniform_state_approaches_equilibrium

Uses:
- Temperature = 3 K
- psi = 1
- dt = 0.001
- 100,000 steps

Expected equilibrium:

|psi| = sqrt(max(1 - T/Tc, 0))

Proves:
- Long-time uniform TDGL evolution approaches the expected equilibrium amplitude.

Classification:
TDGL long-time integration/physics test.

Scientific significance:
Very high.

Important limitation:
The simulation performs 100,000 explicit steps, so this is relatively expensive. It also tests only a single physical condition.

Assessment:
Keep, but this may eventually be optimized or supplemented with a more systematic convergence test.


### test_tdgl_solver_enforces_insulating_boundaries

Creates a uniform psi = 1 state and deliberately perturbs:
- psi[:,1] = 0.8 + 0.2i

Performs one TDGL step.

Then evaluates the left-boundary covariant derivative.

Expected:
- D_x psi = 0

Proves:
- The TDGL solver applies the insulating boundary condition during evolution.

Classification:
TDGL solver/boundary integration test.

Scientific significance:
High.

Important limitation:
Only the left boundary is explicitly checked even though all four boundaries are configured as insulating.

Also, the test derives the dimensionless mesh spacing using:

mesh.dx / coherence_length

This makes the test dependent on the current scaling convention.

Assessment:
Keep.

Future improvement:
Check all four boundaries and use the same operator/discretization used internally by the solver to avoid testing a subtly different formulation.


Overall assessment of test_tdgl_solver.py:

This is the most important integration-level TDGL test file currently present.

It verifies:
- one-step stability
- analytical transient evolution
- above-Tc decay
- long-time equilibrium
- boundary enforcement

The analytical uniform-state transient test is particularly valuable because it is an actual quantitative comparison against an independent solution rather than merely checking whether a value moved in the expected direction.

Major future improvements:
- Temporal convergence study.
- Spatial convergence study.
- Multiple temperatures.
- Multiple initial amplitudes.
- Known analytical phase-gradient/current cases.
- Tests with nonzero vector potential.
- Tests with nonuniform temperature.
- Eventually vortex and magnetic-field benchmarks.


============================================================
## tests/tdgl/test_tdgl_thermal_coupling.py
============================================================

Purpose:
Tests whether a localized increase in temperature suppresses the superconducting order parameter.

Test:
- test_local_temperature_suppresses_order_parameter

Setup:
- NbN geometry
- NbN material
- Standard mesh
- Initial temperature = 3 K
- Uniform initial psi = 1
- Circular hotspot at the mesh center
- Radius = 10 cells
- Hotspot temperature = 0.95 Tc
- Insulating boundaries
- dt = 0.001
- 1000 TDGL steps

The hotspot is created using:

(X - center_x)^2 + (Y - center_y)^2 < radius^2

The expected local equilibrium amplitude is:

|psi| = sqrt(1 - T/Tc)

For the hotspot:

T = 0.95 Tc

so the expected local equilibrium amplitude is:

sqrt(0.05)

The outside region remains at:
- 3 K

The test compares:
- mean |psi| inside hotspot
- mean |psi| outside hotspot

Proves:
- Higher local temperature suppresses the superconducting order parameter relative to the colder surrounding region.

Classification:
TDGL/thermal coupling integration test.

Scientific significance:
High qualitatively.

Important diagnostic output:
The test also calculates:
- mean covariant-Laplacian magnitude inside hotspot
- mean nonlinear term magnitude inside hotspot

These are printed for diagnostic purposes but are not currently used as assertions.

Important limitation:
The test does NOT quantitatively assert that the hotspot reaches its analytical equilibrium amplitude.

It only asserts:

hotspot |psi| < outside |psi|

The printed suppression ratio is also not currently asserted.

Therefore the test establishes the expected direction of thermal suppression but not quantitative agreement with the local equilibrium prediction.

Additional limitation:
The hotspot temperature is imposed directly and is not dynamically generated by the thermal solver. Therefore this is a one-way thermal-to-TDGL coupling test, not yet a full electrothermal feedback test.

Additional limitation:
The hotspot radius is specified in mesh cells rather than physical units. Its physical size therefore depends on the current mesh resolution.

Assessment:
Keep.

This is an important bridge between the TDGL system and the existing thermal infrastructure.

Future improvements:
- Assert a quantitative suppression threshold.
- Compare hotspot amplitude against sqrt(0.05).
- Verify the spatial profile around the hotspot.
- Test several hotspot temperatures.
- Test hotspot removal and recovery.
- Couple to the actual thermal solver.
- Eventually test bidirectional feedback:
  temperature -> superconductivity -> current -> Joule heating -> temperature.


============================================================
# Overall Assessment of tests/tdgl/
============================================================

The TDGL test suite is now substantially more mature than a simple collection of sanity checks.

It covers five major levels:

1. Mathematical operators
   - gradient
   - Laplacian
   - covariant gradient
   - covariant Laplacian

2. Gauge structure
   - gauge links
   - gauge transformations
   - gauge covariance of derivatives
   - gauge covariance of the Laplacian

3. TDGL model physics
   - equilibrium order parameter
   - temperature dependence
   - phase
   - supercurrent
   - free energy

4. Numerical scaling
   - physical/dimensionless conversion
   - GL scales
   - vector potential
   - magnetic field
   - current density

5. Solver/integration behavior
   - stability
   - analytical uniform-state evolution
   - long-time equilibrium
   - above-Tc decay
   - insulating boundary enforcement
   - thermal suppression


============================================================
## Strongest tests currently present
============================================================

The most scientifically meaningful tests are:

1. test_gauge_transformation_preserves_covariant_gradient

Demonstrates gauge covariance of the discrete covariant gradient.

2. test_covariant_laplacian_gauge_transformation

Demonstrates gauge covariance of the discrete covariant Laplacian.

3. test_supercurrent_from_phase_gradient

Demonstrates that superconducting phase gradients generate current.

4. test_supercurrent_from_vector_potential

Demonstrates coupling between vector potential and superconducting current.

5. test_tdgl_uniform_state_matches_analytical_solution

Quantitatively compares numerical TDGL evolution against an analytical solution.

6. test_tdgl_uniform_state_approaches_equilibrium

Tests long-time convergence toward the expected GL equilibrium.

7. test_current_density_scale

Tests one of the fundamental physical normalization relationships.

8. test_local_temperature_suppresses_order_parameter

Demonstrates qualitative coupling between thermal state and superconductivity.


============================================================
## Current weaknesses
============================================================

The largest weakness is not lack of basic tests. The suite already has good coverage of the basic infrastructure.

The major missing component is systematic quantitative validation.

Many tests currently establish qualitative behavior:

- current is zero
- current is nonzero
- temperature suppresses psi
- psi decays above Tc
- equilibrium amplitude goes to zero
- boundaries remain insulating

The strongest solver test goes further and compares numerical evolution against an analytical solution.

The next stage should expand this style of testing.

Important future validation categories:

1. Temporal convergence

Run the same analytical problem at several dt values.

Determine whether the numerical error decreases at the expected rate.

2. Spatial convergence

Repeat operator and solver tests at multiple mesh resolutions.

Determine whether the discretization converges toward the continuum result.

3. Gauge invariance

Expand beyond the current covariance tests to verify that observable quantities such as:
- |psi|
- supercurrent
- magnetic field
- free energy

are invariant under equivalent gauge transformations.

4. Physical scaling validation

Verify the exact nondimensionalization used by the project against the governing dimensional TDGL equations.

5. Current-density validation

Test the full current relation with:
- phase gradient
- vector potential
- order-parameter amplitude
- simultaneous phase and vector potential

6. Magnetic-field/vortex validation

Eventually introduce benchmark configurations involving:
- applied magnetic field
- flux quantization
- vortices
- vortex motion

7. Thermal feedback

The current thermal-coupling test is one-way.

The eventual coupled model must establish:

temperature
    ->
order parameter
    ->
supercurrent/electrical response
    ->
Joule heating
    ->
temperature


============================================================
## Important current test-suite interpretation
============================================================

The TDGL infrastructure should currently be considered:

MATHEMATICALLY / NUMERICALLY TESTED:
- Basic order-parameter operations
- Equilibrium GL relationship
- Ordinary derivatives
- Covariant derivatives
- Covariant Laplacian
- Gauge-link magnitude
- Gauge covariance
- Basic supercurrent behavior
- Basic free-energy calculations
- Physical/dimensionless scaling conversions
- Basic TDGL solver evolution
- Basic thermal suppression

NOT YET FULLY PHYSICALLY VALIDATED:
- Complete dimensional TDGL formulation
- Exact nondimensionalization
- Spatial convergence
- Temporal convergence
- Magnetic-field dynamics
- Vortex physics
- Electromagnetic self-consistency
- Quantitative electrothermal feedback
- Experimental validation
- Research-grade accuracy under realistic mesh/time-step conditions


============================================================
## Important mesh consideration
============================================================

The standard NbN geometry currently uses approximately:
- 5 um x 5 um physical domain
- 100 x 100 mesh
- coherence length xi = 5 nm

This makes the numerical cell size approximately 50 nm, or roughly 10 coherence lengths.

This is very coarse for resolving spatial TDGL structure on the coherence-length scale.

Therefore:

Passing the current TDGL tests does NOT demonstrate that realistic spatial TDGL phenomena are resolved.

The current tests are nevertheless valuable because many of them use uniform states where spatial resolution is not the dominant issue.

This distinction should be preserved:

The current TDGL test suite demonstrates that the implemented equations and numerical infrastructure behave correctly in controlled cases.

It does not yet demonstrate that the production mesh is sufficiently resolved for realistic superconducting structures, vortices, or localized order-parameter variations.


============================================================
## Recommended next TDGL testing phase
============================================================

Do NOT immediately rewrite the existing tests.

First preserve the current suite as the baseline.

The next testing phase should add quantitative validation in approximately this order:

1. Temporal convergence of the uniform analytical TDGL solution.

2. Spatial convergence of the ordinary and covariant operators.

3. Gauge invariance of observable quantities.

4. Quantitative hotspot equilibrium comparison.

5. Supercurrent validation with combined phase gradient and vector potential.

6. More complete boundary-condition tests on all four sides.

7. Applied magnetic-field benchmark.

8. Vortex benchmark.

9. Full TDGL/thermal feedback validation.

The existing tests should remain as regression tests while these stronger validation tests are added.


============================================================
## Current status of tests/tdgl/
============================================================

Status:
TDGL infrastructure is substantially implemented and has a meaningful regression-test foundation.

Completed/tested:
- TDGL parameter validation
- Equilibrium order parameter
- Initialization
- Gauge links
- Covariant gradient
- Covariant Laplacian
- Gauge covariance
- Supercurrent
- Free energy diagnostics
- Physical scaling
- Insulating TDGL boundaries
- Basic TDGL time stepping
- Analytical uniform-state evolution
- Long-time equilibrium
- Above-Tc behavior
- One-way thermal suppression

Not yet completed:
- Systematic convergence validation
- Quantitative spatially varying benchmarks
- Vortex physics validation
- Magnetic-field validation
- Fully coupled electrothermal TDGL validation
- Experimental validation

Overall assessment:
The TDGL test suite is in good shape for the current development stage.

The implementation has moved beyond merely testing whether functions run. Several tests now verify meaningful mathematical properties of the TDGL formulation, especially gauge covariance and analytical time evolution.

The next major step should therefore be quantitative validation rather than simply adding more basic sanity tests.