# SHS Test Suite Context

## Purpose

Compact inventory of the SHS test suite. This document records what each test file tests, what it depends on, what it proves, and which hard-coded values may eventually warrant configuration.

This is an **audit document**, not the authoritative description of implementation. Source code remains authoritative.

---

# Test Inventory

## `tests/conftest.py`

**Purpose:**
Provides a shared pytest fixture for constructing the standard NbN hotspot simulation.

**Tests:**
None. This is test infrastructure.

**Provides / Proves:**
Provides `nbn_simulation` to tests that need a fully constructed SHS simulation.

**Dependencies:**

* `shs.config.builder.build_simulation`
* `configs/simulations/nbn_hotspot_test.json`

**Hard-coded controls:**

* `nbn_hotspot_test.json` path — **CONFIG/TEST FIXTURE**
* No physical values directly embedded in this file.

**Classification:**
Test infrastructure / fixture.

**Assessment:**
Keep. The fixture avoids repeatedly constructing the same standard simulation manually.

---

## `tests/test_boundaries.py`

**Purpose:**
Tests creation and basic management of boundary conditions and boundary sets.

**Tests:**

* `test_boundary_creation` — verifies a `BoundaryCondition` stores its side, type, and temperature correctly.
* `test_boundary_set` — verifies a `BoundarySet` can contain multiple boundary conditions and retrieve them correctly.
* `test_boundary_add` — verifies a boundary can be added to an initially empty `BoundarySet`.

**Proves:**

* Boundary conditions can be constructed.
* Boundary side/type information is preserved.
* Fixed-temperature boundaries can store temperature values.
* Insulating boundaries can be represented.
* `BoundarySet` can store and retrieve boundary conditions.

**Dependencies:**

* `shs.boundaries.BoundaryCondition`
* `BoundarySide`
* `BoundaryType`
* `BoundarySet`

**Hard-coded controls:**

* `4.2 K` — **INTERNAL TEST VALUE** representing a boundary temperature.
* Specific boundary sides/types — **INTERNAL TEST CASE DEFINITION**.

**Classification:**
Software correctness / data-model test.

**Scientific significance:**
Low–moderate. Establishes correct representation of boundary conditions but does not validate their physical application to a solver.

**Assessment:**
Likely keep. These are focused unit tests for the boundary-condition data structures.

---

## `tests/test_boundary_database.py`

**Purpose:**
Tests loading boundary conditions from the simulation configuration.

**Tests:**

* `test_load_boundary_json` — verifies the configured left boundary is loaded as fixed temperature at `4.2 K`.
* `test_insulating_boundary` — verifies the configured top boundary is loaded as insulating.

**Proves:**

* Boundary configuration can be loaded from the NbN simulation JSON.
* Boundary types and associated temperature values are transferred correctly from configuration into boundary objects.

**Dependencies:**

* `shs.boundaries.load_boundary_set`
* `BoundarySide`
* `BoundaryType`
* `configs/simulations/nbn_hotspot_test.json`

**Hard-coded controls:**

* Configuration file path — **CONFIG/TEST FIXTURE**
* `4.2 K` — **CONFIG-DERIVED EXPECTED VALUE**, duplicated in the test.
* Expected boundary types — **CONFIG-DERIVED EXPECTATIONS**.

**Classification:**
Configuration integration / software correctness.

**Scientific significance:**
Low. Tests configuration loading rather than physical boundary behavior.

**Assessment:**
Keep for now. Possible future cleanup if configuration tests become overly repetitive.

---

## `tests/test_builder.py`

**Purpose:**
Tests that the simulation builder correctly constructs the major components of an SHS simulation from configuration.

**Tests:**

* `test_build_simulation` — verifies geometry, mesh, region map, material map, boundaries, and fields are constructed.
* `test_builder_material` — verifies the constructed material is NbN.
* `test_builder_boundary` — verifies the configured left boundary temperature is `4.2 K`.
* `test_builder_fields` — verifies the initial temperature is `3.0 K` and field dimensions match the mesh dimensions.

**Proves:**

* `build_simulation()` can construct the core simulation object.
* Major simulation subsystems are connected during construction.
* Material configuration reaches the material map.
* Boundary configuration reaches the simulation.
* Initial temperature reaches the fields.
* Field array dimensions correspond to the mesh.

**Dependencies:**

* `shs.config.builder.build_simulation`
* `configs/simulations/nbn_hotspot_test.json`
* `BoundarySide`
* Simulation geometry/mesh/mapping/material/boundary/field systems indirectly through the builder.

**Hard-coded controls:**

* Configuration path — **CONFIG/TEST FIXTURE**
* `"NbN"` — **CONFIG-DERIVED EXPECTATION**
* `4.2 K` — **CONFIG-DERIVED EXPECTATION**
* `3.0 K` — **CONFIG-DERIVED EXPECTATION**
* `(mesh.ny, mesh.nx)` — **INTERNAL STRUCTURAL EXPECTATION**

**Classification:**
Integration / software correctness.

**Scientific significance:**
Moderate at most. Demonstrates that configuration is propagated into the simulation state, but does not validate the physics.

**Assessment:**
Keep. Important integration test for the simulation construction pipeline. Some assertions overlap with dedicated configuration/boundary tests and should be reviewed during the final audit.

---

## `tests/test_config.py`

**Purpose:**
Tests loading the primary simulation configuration and verifies selected configured parameters.

**Tests:**

* `test_simulation_config` — verifies temperature, current, and timestep loaded from the NbN hotspot configuration.

**Proves:**

* Simulation configuration can be loaded.
* Selected physical/simulation parameters are correctly parsed.

**Dependencies:**

* `shs.config.simulation.load_simulation`
* `configs/simulations/nbn_hotspot_test.json`

**Hard-coded controls:**

* Configuration path — **CONFIG/TEST FIXTURE**
* `temperature = 3.0 K` — **CONFIG-DERIVED EXPECTATION**
* `current = 0.001 A` — **CONFIG-DERIVED EXPECTATION**
* `dt = 1e-9` — **CONFIG-DERIVED EXPECTATION**

**Classification:**
Configuration/software correctness.

**Scientific significance:**
Low. Establishes configuration parsing, not whether the selected values are physically appropriate.

**Assessment:**
Keep, but review for overlap with builder/integration tests.

---

## `tests/test_contact_map.py`

**Purpose:**
Tests creation and structural validity of the contact masks generated from the NbN geometry and mesh.

**Tests:**

* `test_contact_map_exists` — verifies at least one contact mask exists.
* `test_contact_masks_are_boolean` — verifies masks use Boolean arrays.
* `test_contact_mask_shape` — verifies masks match mesh dimensions.
* `test_contact_contains_cells` — verifies each contact contains at least one mesh cell.
* `test_contact_map_shapes` — repeats the mesh-shape validation.
* `test_current_contacts_exist` — verifies `left_current` and `right_current` contacts exist.
* `test_contact_types` — verifies `left_current` is classified as a current contact.

**Proves:**

* A contact map can be constructed from geometry and mesh.
* Contact masks have the expected Boolean representation.
* Contact masks are spatially compatible with the mesh.
* Contacts actually contain mesh cells.
* Expected current contacts exist.
* Current contact classification is preserved.

**Dependencies:**

* `numpy`
* `shs.geometry.database.load_geometry`
* `shs.geometry.mesh.create_mesh`
* `shs.mapping.build_contact_map`
* `configs/geometry/NbN_film.json`

**Hard-coded controls:**

* Geometry configuration path — **CONFIG/TEST FIXTURE**
* `"left_current"` / `"right_current"` — **CONFIG-DERIVED / DOMAIN-SPECIFIC EXPECTATIONS**
* `"current"` — **CONFIG-DERIVED EXPECTATION**
* Boolean dtype and mesh dimensions — **INTERNAL STRUCTURAL EXPECTATIONS**

**Classification:**
Geometry/mapping correctness / integration.

**Scientific significance:**
Low–moderate. Establishes that physical contacts are correctly represented on the numerical mesh, but does not yet test whether electrical boundary conditions produce physically correct transport.

**Assessment:**
Keep the core structural tests. There is clear redundancy between `test_contact_mask_shape` and `test_contact_map_shapes`; this should be reviewed during the final audit.

---

# SHS Test Suite Context — Continued

## `tests/test_coupled_solver.py`

**Purpose:**
Tests that the coupled electrothermal solver can execute a short simulation and produce Joule heating.

**Tests:**

* `test_coupled_electrothermal` — constructs an NbN simulation, runs five coupled steps, and checks basic output validity.

**Proves:**

* Geometry, mesh, material, mappings, fields, and thermal model can be assembled manually into a `Simulation`.
* `run_coupled_simulation()` can execute multiple steps.
* The reported step count is correct.
* Temperature and heat-source fields retain mesh-compatible dimensions.
* The coupled run produces nonzero heat generation.

**Dependencies:**

* `shs.solvers.run_coupled_simulation`
* `shs.physics.Fields`
* Geometry database / mesh
* Region/material/contact mapping
* Material database
* `shs.physics.thermal.ThermalModel`
* `shs.config.simulation_state.Simulation`
* `configs/geometry/NbN_film.json`

**Hard-coded controls:**

* Initial temperature `3.0 K` — **PHYSICAL/TEST CONDITION**
* Bath temperature `3.0 K` — **PHYSICAL/TEST CONDITION**
* Thermal relaxation rate `1.0` — **SIMULATION CONTROL**
* `steps=5` — **TEST CONTROL**
* `dt=1e-6` — **SIMULATION CONTROL**
* `"NbN"` — **MATERIAL SELECTION**
* Geometry configuration path — **TEST FIXTURE**

**Classification:**
Coupled numerical/integration test with a basic physical-behavior assertion.

**Scientific significance:**
Moderate. Demonstrates that the coupled solver produces heating, but does not establish that the magnitude or evolution of the heating is physically correct.

**Assessment:**
Keep for now. The hard-coded timestep and thermal parameters should be reviewed during the configuration audit. The test is also likely related to the more sophisticated coupled tests elsewhere in the suite.

---

## `tests/test_electrical_solver.py`

**Purpose:**
Tests that the electrical solver produces correctly shaped electrical fields and a finite transport/heating response.

**Tests:**

* `test_electrical_transport` — runs `electrical_step()` and checks voltage, current-density, electric-field, and heat-source fields.

**Proves:**

* `electrical_step()` executes successfully.
* Electrical output fields have the expected mesh dimensions.
* The current configuration produces nonzero Joule heating.
* Voltage remains within the assumed `0–1 V` range.
* Current density is nonzero.

**Dependencies:**

* Geometry database / mesh
* `shs.physics.Fields`
* `shs.solvers.electrical_step`
* Region/material/contact mapping
* Material database
* `configs/geometry/NbN_film.json`

**Hard-coded controls:**

* Initial temperature `3.0 K` — **PHYSICAL/TEST CONDITION**
* Voltage bounds `0.0` and `1.0 V` — **TEST ASSUMPTION / POTENTIALLY QUESTIONABLE**
* `"NbN"` — **MATERIAL SELECTION**
* Geometry configuration path — **TEST FIXTURE**

**Classification:**
Electrical solver integration / numerical behavior test.

**Scientific significance:**
Moderate, but limited. It establishes that the solver produces a finite transport response and heating. It does **not** validate the actual current-voltage relationship, current conservation, or physical magnitude of the solution.

**Assessment:**
Keep for now, but this test requires later review. In particular, `voltage.max() <= 1.0` appears to be an arbitrary assumption rather than a meaningful physics benchmark. The duplicated `heat_source.max() > 0` assertion should also be removed or consolidated eventually.

---

## `tests/test_electromagnetics.py`

**Purpose:**
Tests the default configuration/state of the electromagnetic model.

**Tests:**

* `test_default_electromagnetic_model` — verifies default reference voltage, ground voltage, and self-field setting.

**Proves:**

* `ElectromagneticModel` can be constructed.
* Its default values are assigned as expected.
* Self-field calculation is disabled by default.

**Dependencies:**

* `shs.physics.ElectromagneticModel`

**Hard-coded controls:**

* `reference_voltage = 1.0 V` — **DEFAULT MODEL PARAMETER**
* `ground_voltage = 0.0 V` — **DEFAULT MODEL PARAMETER**
* `include_self_field = False` — **DEFAULT MODEL PARAMETER**

**Classification:**
Software correctness / model configuration.

**Scientific significance:**
Low. This does not validate electromagnetic physics.

**Assessment:**
Keep for now. Later review whether these defaults represent genuine model defaults or user-configurable simulation conditions.

**Important status note:**
The existence of this test does **not** mean self-consistent electromagnetic evolution is implemented or validated. It only tests the electromagnetic model's default state.

---

## `tests/test_fields.py`

**Purpose:**
Tests creation and initialization of the shared simulation field container.

**Tests:**

* `test_create_fields` — verifies basic field dimensions, initial temperature, complex `psi`, and initial order-parameter amplitude/phase.
* `test_electromagnetic_fields` — verifies selected electromagnetic/current fields have mesh-compatible dimensions.

**Proves:**

* `Fields.create()` constructs the expected arrays.
* Temperature and voltage fields have the expected dimensions.
* Initial temperature is correctly assigned.
* `psi` exists as a complex-valued field.
* Initial `psi` has unit amplitude and zero phase.
* Selected electric/current/magnetic fields exist with mesh-compatible dimensions.

**Dependencies:**

* `shs.geometry.load_geometry`
* `shs.geometry.create_mesh`
* `shs.physics.Fields`
* `numpy`
* `configs/geometry/NbN_film.json`

**Hard-coded controls:**

* `initial_temperature=3.0 K` — **PHYSICAL/TEST CONDITION**
* `(100, 100)` — **HARDCODED MESH ASSUMPTION**
* Initial `|psi| = 1` — **PHYSICS/INITIAL-CONDITION ASSUMPTION**
* Initial phase `0` — **PHYSICS/INITIAL-CONDITION ASSUMPTION**

**Classification:**
Data-model / initialization correctness with basic TDGL initialization validation.

**Scientific significance:**
Moderate for the TDGL infrastructure. The test verifies the intended initial order-parameter representation, but does not establish that `psi=1` is the physically correct equilibrium state for the supplied `3 K` temperature.

**Assessment:**
Keep the structural and complex-field tests. The hard-coded `(100,100)` should likely be replaced with mesh-derived dimensions, as the test already does elsewhere. The physical validity of the `psi=1` initialization should be reviewed separately.

---

## `tests/test_geometry.py`

**Purpose:**
Tests loading the NbN geometry and verifies that its major components and contact information exist.

**Tests:**

* `test_geometry_contacts` — verifies the first contact is `left_current` and has type `current`.
* `test_geometry_container` — verifies film, regions, contacts, width, and height exist.

**Proves:**

* Geometry configuration can be loaded.
* The geometry contains a film, regions, and contacts.
* Film dimensions are represented.
* At least the first contact has the expected name and type.

**Dependencies:**

* `shs.geometry.load_geometry`
* `shs.geometry.create_mesh` (imported but not used)
* `configs/geometry/NbN_film.json`

**Hard-coded controls:**

* `"left_current"` — **CONFIG/DOMAIN-SPECIFIC EXPECTATION**
* `"current"` — **CONFIG/DOMAIN-SPECIFIC EXPECTATION**
* No numerical physical values directly tested.

**Classification:**
Geometry/configuration correctness.

**Scientific significance:**
Low. Establishes that the expected geometry structure exists but does not validate geometric dimensions or spatial correctness.

**Assessment:**
Keep the useful geometry-container checks. Review the unused `create_mesh` import. The assumption that `geometry.contacts[0]` is specifically the left current contact may be fragile if contact ordering changes; a lookup by name may eventually be preferable.
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## Cross-Test Observations Added So Far

### Emerging Hard-Coded Controls

The following values have now appeared repeatedly:

| Value                         | Appears in                                                 | Preliminary classification                                          |
| ----------------------------- | ---------------------------------------------------------- | ------------------------------------------------------------------- |
| `3.0 K`                       | builder, config, coupled solver, electrical solver, fields | Physical/test condition; likely configuration-derived               |
| `4.2 K`                       | boundaries, boundary database, builder                     | Boundary/test condition; configuration-derived in integration tests |
| `0.001 A`                     | config                                                     | Simulation operating condition                                      |
| `1e-9`                        | config                                                     | Simulation timestep                                                 |
| `1e-6`                        | coupled solver                                             | Simulation timestep                                                 |
| `1.0` thermal relaxation rate | coupled solver                                             | Simulation control                                                  |
| `100 × 100`                   | fields                                                     | Mesh assumption; should generally be mesh-derived                   |
| `0–1 V`                       | electrical solver                                          | Potentially arbitrary test assumption                               |
| `psi = 1`                     | fields                                                     | Initial-condition assumption requiring physical-context review      |

These are **not yet recommendations to move them into configuration**. They are items for the later configuration audit.

### Emerging Redundancy

Already identified:

* `test_contact_mask_shape`
* `test_contact_map_shapes`

appear to duplicate the same assertion.

There is also increasing overlap between:

* `test_config`
* `test_boundary_database`
* `test_builder`
* `test_fields`

because they all test pieces of the same configuration → simulation construction pipeline. This may be perfectly reasonable, but we should evaluate the overlap once the complete suite is mapped.

### Important Testing Distinction Emerging

The suite is beginning to reveal three different levels of testing:

**1. Unit/software correctness**

Example:

> Does `BoundarySet.add()` actually store a boundary?

**2. Integration/numerical sanity**

Example:

> Can the electrical solver run and produce finite current and heating?

**3. Scientific validation**

Example:

> Does the numerical TDGL solution reproduce an analytical TDGL solution?

The final audit should explicitly preserve this distinction. A passing integration test should **not** be counted as equivalent to a physics validation test.

# SHS Test Suite Context — Continued

## `tests/test_material_map.py`

**Purpose:**
Tests construction and basic correctness of the material-property mapping from geometry regions to mesh cells.

**Tests:**

* `test_material_map_shape` — verifies material IDs match mesh dimensions.
* `test_material_dictionary` — verifies material ID `0` corresponds to NbN.
* `test_property_arrays` — verifies mapped thermal/electrical/critical-temperature properties match the source material.
* `test_all_cells_are_material_zero` — verifies every cell is assigned material ID `0`.
* `test_conductivity_map` — verifies electrical conductivity is the reciprocal of normal resistivity.

**Proves:**

* Material maps have correct spatial dimensions.
* Materials are associated with material IDs.
* Material properties are transferred onto the mesh.
* The current single-material NbN geometry assigns every cell to NbN.
* Normal electrical conductivity is calculated as `1 / normal_resistivity`.

**Dependencies:**

* `numpy`
* Geometry database
* Mesh generation
* Region mapping
* Material mapping
* Material database
* `configs/geometry/NbN_film.json`

**Hard-coded controls:**

* `"NbN"` — **MATERIAL SELECTION**
* Material ID `0` — **CURRENT MAPPING ASSUMPTION**
* `"NbN_film.json"` — **TEST FIXTURE**
* Property comparisons — **SOURCE-MATERIAL-DERIVED EXPECTATIONS**

**Classification:**
Mapping correctness / integration.

**Scientific significance:**
Moderate. Establishes that material properties are correctly represented numerically, but does not validate whether the material parameters themselves are physically correct.

**Assessment:**
Keep. Particularly useful as the project eventually expands to multiple materials/regions. `test_all_cells_are_material_zero` is specific to the current single-material geometry and should be reconsidered if the geometry becomes more complex.

---

## `tests/test_materials.py`

**Purpose:**
Tests basic retrieval and validity of the NbN material definition.

**Tests:**

* `test_load_nbn` — verifies NbN can be retrieved and has positive critical temperature and thickness.

**Proves:**

* The NbN material exists in the material database.
* Basic required material properties are populated with positive values.

**Dependencies:**

* `shs.materials.get_material`

**Hard-coded controls:**

* `"NbN"` — **MATERIAL SELECTION**
* `Tc > 0` and `thickness > 0` — **MINIMAL VALIDITY REQUIREMENTS**

**Classification:**
Software correctness / material database sanity test.

**Scientific significance:**
Low. It does not establish that the NbN parameters are experimentally accurate.

**Assessment:**
Keep, although the test is minimal. The following module-level code is not part of the test itself:

```text
nbn = get_material("NbN")
print(nbn)
```

This produces output merely by importing the test module and should likely be removed from the permanent automated suite.

---

## `tests/test_numerics.py`

**Purpose:**
Tests fundamental numerical differential operators and iterative solvers.

**Tests:**

* `test_gradient` — verifies the numerical gradient of a simple linear field.
* `test_laplacian` — verifies the Laplacian of a constant field is zero.
* `test_gauss_seidel` — verifies Gauss-Seidel converges on a simple boundary-value problem.
* `test_red_black_sor` — verifies Red-Black SOR converges on a simple boundary-value problem.

**Proves:**

* Gradient correctly handles a simple linear field.
* Laplacian correctly handles a constant field.
* Gauss-Seidel can converge on the constructed test problem.
* Red-Black SOR can converge on the constructed test problem.
* Boundary values are retained by the iterative solvers.

**Dependencies:**

* `numpy`
* `shs.numerics.gradient`
* `shs.numerics.laplacian`
* `shs.numerics.gauss_seidel`
* `shs.numerics.red_black_sor`

**Hard-coded controls:**

* Small test-grid sizes `(5,5)`, `(10,10)`, `(20,20)` — **TEST PROBLEM DEFINITION**
* `dx=1`, `dy=1` — **TEST PROBLEM DEFINITION**
* Boundary value `1` — **TEST PROBLEM DEFINITION**
* `max_iterations` is not explicitly specified for these tests — solver default is therefore being tested indirectly.

**Classification:**
Numerical correctness / unit testing.

**Scientific significance:**
Moderate foundational significance. These tests establish mathematical behavior of numerical infrastructure but are not superconductivity-specific validation.

**Assessment:**
Keep. These are valuable foundational tests. Later review whether additional analytical numerical tests are needed, especially for convergence accuracy and boundary-condition behavior.

---

## `tests/test_region_map.py`

**Purpose:**
Tests mapping of geometry regions onto the numerical mesh.

**Tests:**

* `test_region_map_shape` — verifies region IDs match mesh dimensions.
* `test_region_map_all_film` — verifies all current cells belong to region `0`.
* `test_region_dictionary` — verifies region ID `0` is named `"film"`.
* `test_region_type` — verifies region `0` is classified as `"superconductor"`.

**Proves:**

* Region maps have the expected dimensions.
* The current NbN film geometry is mapped to the mesh.
* Region IDs have associated names.
* The film is classified as a superconductor.

**Dependencies:**

* `numpy`
* Geometry database
* Mesh generation
* `build_region_map`
* `configs/geometry/NbN_film.json`

**Hard-coded controls:**

* Region ID `0` — **CURRENT GEOMETRY ASSUMPTION**
* `"film"` — **CONFIG/DOMAIN-SPECIFIC EXPECTATION**
* `"superconductor"` — **CONFIG/DOMAIN-SPECIFIC EXPECTATION**
* Geometry path — **TEST FIXTURE**

**Classification:**
Geometry/mapping correctness.

**Scientific significance:**
Moderate. Correct region classification is important to later material and superconducting physics, but this test does not validate physical behavior.

**Assessment:**
Keep the structural tests. `test_region_map_all_film` is specifically tied to the current simple geometry and should be reconsidered once multiple regions are introduced.

---

## `tests/test_simulation_boundaries.py`

**Purpose:**
Tests loading boundary-condition information through the simulation configuration system.

**Tests:**

* `test_simulation_loads_boundaries` — verifies the left boundary exists and is fixed-temperature.
* `test_boundary_temperature` — verifies the configured left boundary temperature is `4.2 K`.

**Proves:**

* Simulation configuration contains boundary information.
* Boundary type is parsed correctly.
* Boundary temperature is parsed correctly.

**Dependencies:**

* `shs.config.simulation.load_simulation`
* `configs/simulations/nbn_hotspot_test.json`

**Hard-coded controls:**

* Configuration path — **TEST FIXTURE**
* `"left"` — **CONFIG/DOMAIN-SPECIFIC EXPECTATION**
* `"fixed_temperature"` — **CONFIG-DERIVED EXPECTATION**
* `4.2 K` — **CONFIG-DERIVED EXPECTATION**

**Classification:**
Configuration correctness.

**Scientific significance:**
Low. It tests configuration parsing rather than application of the boundary condition to the physical solver.

**Assessment:**
Keep for now, but likely overlaps substantially with `test_boundary_database.py` and parts of `test_builder.py`. Final audit should determine whether all three layers provide enough distinct coverage to justify retaining them.

---

## `tests/test_solver_comparison.py`

**Purpose:**
Compares Gauss-Seidel and Red-Black SOR on the same numerical boundary-value problem.

**Tests:**

* `create_test_problem` — constructs a standardized `100 × 100` solver problem with fixed left/right voltage boundaries.
* `test_solver_comparison` — runs both solvers, checks convergence/residuals, compares their resulting fields, and prints performance diagnostics.

**Proves:**

* Both Gauss-Seidel and Red-Black SOR converge on the same test problem.
* Both achieve residual below `1e-7`.
* Their resulting fields agree to within `1e-4`.
* Runtime and iteration counts can be observed for comparison.

**Does NOT currently prove:**

* That Red-Black SOR is faster by an asserted threshold.
* That either solver is physically correct for the electrical model.
* That the solution is the analytical solution to this particular boundary-value problem.

**Dependencies:**

* `time`
* `numpy`
* `shs.numerics.gauss_seidel`
* `shs.numerics.red_black_sor`

**Hard-coded controls:**

* Grid size `100 × 100` — **BENCHMARK PROBLEM DEFINITION**
* `dx=1`, `dy=1` — **BENCHMARK PROBLEM DEFINITION**
* Left boundary `1.0` — **BENCHMARK PROBLEM DEFINITION**
* Right boundary `0.0` — **BENCHMARK PROBLEM DEFINITION**
* Gauss-Seidel `max_iterations=50000` — **SOLVER CONTROL**
* Residual threshold `1e-7` — **NUMERICAL VALIDATION CRITERION**
* Solver agreement threshold `1e-4` — **NUMERICAL VALIDATION CRITERION**

**Classification:**
Numerical validation / performance diagnostic.

**Scientific significance:**
Moderate for numerical infrastructure; low for superconducting physics.

**Assessment:**
Important test/benchmark, but likely better separated from the normal automated test suite eventually. Its `print()` statements and runtime measurements make it partly a diagnostic/benchmark. The test currently validates convergence and agreement, while the performance information is only printed rather than asserted.

---

# Cross-Test Observations — Updated

## Hard-Coded Values Appearing Across the Suite

| Value / assumption            | Current interpretation                                                            |
| ----------------------------- | --------------------------------------------------------------------------------- |
| `3.0 K`                       | Repeated physical/test condition; likely configuration-derived                    |
| `4.2 K`                       | Repeated boundary/test condition; configuration-derived                           |
| `0.001 A`                     | Configured operating condition                                                    |
| `1e-9`                        | Configured timestep                                                               |
| `1e-6`                        | Coupled-solver timestep                                                           |
| `1.0` thermal relaxation rate | Coupled-solver control                                                            |
| `100 × 100`                   | Appears both as current mesh assumption and numerical benchmark size              |
| `0–1 V`                       | Potentially arbitrary electrical test assumption                                  |
| `psi = 1`                     | Initial-condition assumption                                                      |
| `1e-7` residual               | Numerical acceptance criterion; likely **INTERNAL TEST VALUE**, not configuration |
| `1e-4` solver difference      | Numerical comparison criterion; likely **INTERNAL TEST VALUE**                    |
| `50000` iterations            | Benchmark/solver control                                                          |

The distinction between **configuration candidates** and **legitimate test constants** is becoming increasingly important. We should not automatically move numerical tolerances or mathematical test-problem definitions into the simulation configuration.

---

## New Redundancy Candidates

### Region mapping

`test_region_map.py` is structurally clean, but its assumptions overlap with material mapping:

```text
region 0 → film
film → superconductor
all cells → region 0
```

This is appropriate for the current geometry, but several tests will need reconsideration once SHS supports multi-region geometries.

### Material mapping

`test_material_map.py` has several tests that validate different layers of the same mapping operation. Most appear useful, but:

* `test_material_map_shape`
* `test_all_cells_are_material_zero`

are partly consequences of the current simple geometry.

They may need to evolve rather than simply be deleted when more complex geometries are introduced.

### Boundary configuration

There is now significant overlap among:

```text
test_boundaries.py
test_boundary_database.py
test_simulation_boundaries.py
test_builder.py
```

They operate at different architectural layers, however, so **do not remove any yet**.

The final audit should determine whether each layer deserves independent coverage.

---

## Test Quality Issue Identified

`tests/test_materials.py` contains module-level executable code:

```text
nbn = get_material("NbN")
print(nbn)
```

This is not an assertion and does not contribute to automated correctness. It is likely leftover diagnostic code.

**Candidate final action:** remove from the permanent pytest test file or move to a diagnostic script.

---

## Benchmark vs Automated Test

`test_solver_comparison.py` is becoming a clear example of a test that serves **two purposes**:

1. Automated numerical correctness:

   * both solvers converge
   * residuals are sufficiently small
   * solutions agree

2. Performance investigation:

   * iteration counts
   * runtime measurements
   * printed comparison

These purposes may eventually be separated into:

```text
tests/
    test_iterative_solvers.py

diagnostics/
    compare_solvers.py
```

No change should be made yet; this is an audit observation only.

---

## Current Test Categories Identified

The suite now contains examples of:

**Data/model unit tests**

* boundaries
* materials
* fields

**Configuration tests**

* config
* simulation boundaries
* boundary database

**Geometry/mapping tests**

* geometry
* region map
* material map
* contact map

**Numerical tests**

* numerical operators
* Gauss-Seidel
* Red-Black SOR

**Solver/integration tests**

* electrical solver
* coupled solver
* builder

**Benchmark/diagnostic tests**

* solver comparison

**TDGL tests**

* Not yet inventoried in this audit batch.

# SHS Test Suite Context — Continued

## `tests/test_superconducting_transport.py`

**Purpose:**
Tests the relationship between the superconducting current density and the complex order parameter.

**Tests:**

* `test_superconducting_current_zero_when_psi_zero`
* `test_superconducting_current_preserves_current_at_unit_psi`
* `test_superconducting_current_scales_with_order_parameter`

**Proves:**

* Superconducting current vanishes when `psi = 0`.
* At `|psi| = 1`, the supplied current is unchanged.
* For `|psi|² = 0.5`, superconducting current is reduced to 50%.

This establishes that the current implementation applies the expected scaling:

$$
J_s \propto |\psi|^2
$$

for the tested cases.

**Dependencies:**

* `numpy`
* `shs.physics.superconducting_transport.superconducting_current_density`

**Hard-coded controls:**

* `10 × 10` arrays — **TEST SIZE**
* `psi = 0`, `1`, and `sqrt(0.5)` — **ANALYTICAL TEST CASES**
* Current values `1`, `2`, and `-3` — **TEST INPUTS**
* Expected scaling `0.5` — **ANALYTICAL EXPECTATION**

**Classification:**
Physics unit test.

**Scientific significance:**
Moderate. It verifies the implemented mathematical relationship but does not validate the complete physical superconducting-current equation, gauge coupling, units, or electromagnetic consistency.

**Assessment:**
Keep. These are clean limiting-case tests. Later, once gauge-covariant superconducting transport is implemented, these tests will likely need to evolve.

---

## `tests/test_thermal.py`

**Purpose:**
Tests the analytical calculation of thermal diffusivity.

**Test:**

* `test_thermal_diffusivity`

Uses:

$$
\alpha = \frac{k}{C}
$$

with `k = 10` and `C = 2×10^6`, expecting:

$$
\alpha = 5×10^{-6}.
$$

**Proves:**

* The implemented thermal diffusivity calculation gives the expected result for the supplied parameters.

**Dependencies:**

* `ThermalModel`

**Hard-coded controls:**

* Thermal conductivity `10`
* Heat capacity `2e6`
* Expected value `5e-6`

These are **TEST INPUTS**, not necessarily configuration candidates.

**Classification:**
Physics unit test / analytical calculation.

**Scientific significance:**
Low-to-moderate. It verifies the mathematical implementation of a material property calculation.

**Assessment:**
Keep. Very simple and inexpensive test.

---

## `tests/test_thermal_solver.py`

**Purpose:**
Tests thermal evolution under diffusion, heating, and thermal relaxation.

**Shared setup:**

* Builds the NbN simulation from `nbn_hotspot_test.json`.
* Uses a bath temperature of `3 K`.

### `test_heat_diffusion`

Raises the center cell to `10 K` and performs one thermal step with no bath relaxation.

**Proves:**

* A localized temperature perturbation does not increase at the center during the tested diffusion step.
* Temperature field retains the correct mesh shape.

**Hard-coded controls:**

* Center cell `[50,50]`
* Initial hotspot temperature `10 K`
* Bath temperature `3 K`
* Relaxation rate `0`
* `dt = 1e-9`

**Classification:**
Thermal physics / solver behavior.

**Important limitation:**
The test only checks that the center temperature is `<= 10 K`. It does not quantitatively verify the diffusion equation or conservation behavior.

---

### `test_boundary_stability`

Places a heat source at the center and checks the corner temperature remains `3 K`.

**Proves:**

* The tested boundary/corner remains unchanged during the thermal step.

**Hard-coded controls:**

* Heat source `10`
* Location `[50,50]`
* Bath temperature `3 K`
* Relaxation rate `0`
* `dt = 1e-6`

**Classification:**
Thermal boundary-condition test.

**Important limitation:**
The test assumes `[0,0]` represents a boundary that should remain fixed, but it does not explicitly establish which boundary condition is responsible. This should be reviewed.

---

### `test_hotspot_heating`

Applies a very large localized heat source:

$$
Q = 10^{12}
$$

and verifies that the center temperature rises above `3 K`.

**Proves:**

* Positive localized heat input produces a positive temperature response.

**Hard-coded controls:**

* Heat source `1e12`
* Center `[50,50]`
* `dt = 1e-9`
* Bath temperature `3 K`
* Relaxation rate `0`

**Classification:**
Thermal physics sanity test.

**Scientific significance:**
Moderate qualitatively, but the chosen source magnitude is an artificial test stimulus rather than experimental validation.

---

### `test_uniform_temperature_remains_constant`

Runs the thermal solver from a uniform `3 K` state with no relaxation or explicit heating.

**Proves:**

* A uniform temperature field remains uniform under the tested diffusion operation.

**Classification:**
Thermal numerical correctness / limiting case.

**Assessment:**
Strong test. This is an important invariant.

---

### `test_bath_cooling`

Initializes the entire system to `10 K` and couples it to a `3 K` bath.

**Proves:**

* Thermal relaxation lowers the temperature toward the bath.

**Hard-coded controls:**

* Initial temperature `10 K`
* Bath temperature `3 K`
* Relaxation rate `100`
* `dt = 1e-6`

**Classification:**
Thermal physics behavior.

**Limitation:**
Only verifies the direction of temperature change, not the quantitative relaxation rate.

---

### `test_bath_equilibrium`

Starts at the configured `3 K` temperature and applies bath coupling.

**Proves:**

* A system already at bath temperature remains at bath temperature.

**Classification:**
Thermal equilibrium/invariant test.

**Assessment:**
Keep. Important limiting case.

---

### Overall assessment

`test_thermal_solver.py` is an important physics test file and should remain.

The biggest future improvement is **quantitative validation**. Several tests currently establish only:

```text
temperature went up
temperature went down
temperature did not change
```

rather than comparing the result against an analytical thermal solution.

This should be addressed later during the physics-validation stage, not during this initial audit.

---

## `tests/test_validator.py`

**Purpose:**
Tests simulation validation and detection of inconsistent field/property dimensions.

**Tests:**

* `test_invalid_material_map`
* `test_simulation_validation`
* `test_invalid_field_shape`

**Proves:**

* A correctly constructed simulation passes validation.
* An incorrectly shaped material-property array raises `ValueError`.
* An incorrectly shaped temperature field raises `ValueError`.

**Dependencies:**

* `build_simulation`
* `pytest`
* `nbn_hotspot_test.json`

**Hard-coded controls:**

* Configuration path — **TEST FIXTURE**
* Removing one element from arrays — **DELIBERATELY INVALID TEST INPUT**

**Classification:**
Software validation / defensive correctness.

**Scientific significance:**
Low directly, but important for preventing malformed simulation states.

**Assessment:**
Keep. These are legitimate structural validation tests.

---

## `tests/test_hotspot.py`

**Purpose:**
Tests generation of a localized Gaussian optical heat source.

**Setup:**

* Loads the NbN geometry and mesh.
* Creates a `GaussianHotspot`.

Parameters:

$$
x_0 = 2.5\,\mu m
$$

$$
y_0 = 2.5\,\mu m
$$

$$
A = 10^8
$$

$$
\sigma = 0.5\,\mu m
$$

**Proves:**

* The Gaussian hotspot generator produces an array matching the mesh dimensions.
* The generated heat source contains positive values.

**Dependencies:**

* Geometry loading
* Mesh generation
* `GaussianHotspot`

**Hard-coded controls:**

* Center `(2.5 μm, 2.5 μm)` — **PHYSICAL TEST CONDITION**
* Amplitude `1e8` — **TEST CONDITION**
* Width `0.5 μm` — **TEST CONDITION**
* Expected shape `100 × 100` — **CURRENT MESH ASSUMPTION**

**Classification:**
Optical-physics unit/integration test.

**Scientific significance:**
Moderate. It verifies that localized optical heating infrastructure exists, but does not establish that the Gaussian has the correct normalization, spatial width, power interpretation, or physical units.

**Assessment:**
Keep. Later strengthen with tests of:

* hotspot center location,
* symmetry,
* peak position,
* width,
* normalization/integrated power.
# SHS Test Suite Context — Continued

## `tests/test_superconducting_transport.py`

**Purpose:**
Tests the relationship between the superconducting current density and the complex order parameter.

**Tests:**

* `test_superconducting_current_zero_when_psi_zero`
* `test_superconducting_current_preserves_current_at_unit_psi`
* `test_superconducting_current_scales_with_order_parameter`

**Proves:**

* Superconducting current vanishes when `psi = 0`.
* At `|psi| = 1`, the supplied current is unchanged.
* For `|psi|² = 0.5`, superconducting current is reduced to 50%.

This establishes that the current implementation applies the expected scaling:

$$
J_s \propto |\psi|^2
$$

for the tested cases.

**Dependencies:**

* `numpy`
* `shs.physics.superconducting_transport.superconducting_current_density`

**Hard-coded controls:**

* `10 × 10` arrays — **TEST SIZE**
* `psi = 0`, `1`, and `sqrt(0.5)` — **ANALYTICAL TEST CASES**
* Current values `1`, `2`, and `-3` — **TEST INPUTS**
* Expected scaling `0.5` — **ANALYTICAL EXPECTATION**

**Classification:**
Physics unit test.

**Scientific significance:**
Moderate. It verifies the implemented mathematical relationship but does not validate the complete physical superconducting-current equation, gauge coupling, units, or electromagnetic consistency.

**Assessment:**
Keep. These are clean limiting-case tests. Later, once gauge-covariant superconducting transport is implemented, these tests will likely need to evolve.

---

## `tests/test_thermal.py`

**Purpose:**
Tests the analytical calculation of thermal diffusivity.

**Test:**

* `test_thermal_diffusivity`

Uses:

$$
\alpha = \frac{k}{C}
$$

with `k = 10` and `C = 2×10^6`, expecting:

$$
\alpha = 5×10^{-6}.
$$

**Proves:**

* The implemented thermal diffusivity calculation gives the expected result for the supplied parameters.

**Dependencies:**

* `ThermalModel`

**Hard-coded controls:**

* Thermal conductivity `10`
* Heat capacity `2e6`
* Expected value `5e-6`

These are **TEST INPUTS**, not necessarily configuration candidates.

**Classification:**
Physics unit test / analytical calculation.

**Scientific significance:**
Low-to-moderate. It verifies the mathematical implementation of a material property calculation.

**Assessment:**
Keep. Very simple and inexpensive test.

---

## `tests/test_thermal_solver.py`

**Purpose:**
Tests thermal evolution under diffusion, heating, and thermal relaxation.

**Shared setup:**

* Builds the NbN simulation from `nbn_hotspot_test.json`.
* Uses a bath temperature of `3 K`.

### `test_heat_diffusion`

Raises the center cell to `10 K` and performs one thermal step with no bath relaxation.

**Proves:**

* A localized temperature perturbation does not increase at the center during the tested diffusion step.
* Temperature field retains the correct mesh shape.

**Hard-coded controls:**

* Center cell `[50,50]`
* Initial hotspot temperature `10 K`
* Bath temperature `3 K`
* Relaxation rate `0`
* `dt = 1e-9`

**Classification:**
Thermal physics / solver behavior.

**Important limitation:**
The test only checks that the center temperature is `<= 10 K`. It does not quantitatively verify the diffusion equation or conservation behavior.

---

### `test_boundary_stability`

Places a heat source at the center and checks the corner temperature remains `3 K`.

**Proves:**

* The tested boundary/corner remains unchanged during the thermal step.

**Hard-coded controls:**

* Heat source `10`
* Location `[50,50]`
* Bath temperature `3 K`
* Relaxation rate `0`
* `dt = 1e-6`

**Classification:**
Thermal boundary-condition test.

**Important limitation:**
The test assumes `[0,0]` represents a boundary that should remain fixed, but it does not explicitly establish which boundary condition is responsible. This should be reviewed.

---

### `test_hotspot_heating`

Applies a very large localized heat source:

$$
Q = 10^{12}
$$

and verifies that the center temperature rises above `3 K`.

**Proves:**

* Positive localized heat input produces a positive temperature response.

**Hard-coded controls:**

* Heat source `1e12`
* Center `[50,50]`
* `dt = 1e-9`
* Bath temperature `3 K`
* Relaxation rate `0`

**Classification:**
Thermal physics sanity test.

**Scientific significance:**
Moderate qualitatively, but the chosen source magnitude is an artificial test stimulus rather than experimental validation.

---

### `test_uniform_temperature_remains_constant`

Runs the thermal solver from a uniform `3 K` state with no relaxation or explicit heating.

**Proves:**

* A uniform temperature field remains uniform under the tested diffusion operation.

**Classification:**
Thermal numerical correctness / limiting case.

**Assessment:**
Strong test. This is an important invariant.

---

### `test_bath_cooling`

Initializes the entire system to `10 K` and couples it to a `3 K` bath.

**Proves:**

* Thermal relaxation lowers the temperature toward the bath.

**Hard-coded controls:**

* Initial temperature `10 K`
* Bath temperature `3 K`
* Relaxation rate `100`
* `dt = 1e-6`

**Classification:**
Thermal physics behavior.

**Limitation:**
Only verifies the direction of temperature change, not the quantitative relaxation rate.

---

### `test_bath_equilibrium`

Starts at the configured `3 K` temperature and applies bath coupling.

**Proves:**

* A system already at bath temperature remains at bath temperature.

**Classification:**
Thermal equilibrium/invariant test.

**Assessment:**
Keep. Important limiting case.

---

### Overall assessment

`test_thermal_solver.py` is an important physics test file and should remain.

The biggest future improvement is **quantitative validation**. Several tests currently establish only:

```text
temperature went up
temperature went down
temperature did not change
```

rather than comparing the result against an analytical thermal solution.

This should be addressed later during the physics-validation stage, not during this initial audit.

---

## `tests/test_validator.py`

**Purpose:**
Tests simulation validation and detection of inconsistent field/property dimensions.

**Tests:**

* `test_invalid_material_map`
* `test_simulation_validation`
* `test_invalid_field_shape`

**Proves:**

* A correctly constructed simulation passes validation.
* An incorrectly shaped material-property array raises `ValueError`.
* An incorrectly shaped temperature field raises `ValueError`.

**Dependencies:**

* `build_simulation`
* `pytest`
* `nbn_hotspot_test.json`

**Hard-coded controls:**

* Configuration path — **TEST FIXTURE**
* Removing one element from arrays — **DELIBERATELY INVALID TEST INPUT**

**Classification:**
Software validation / defensive correctness.

**Scientific significance:**
Low directly, but important for preventing malformed simulation states.

**Assessment:**
Keep. These are legitimate structural validation tests.

---

## `tests/test_hotspot.py`

**Purpose:**
Tests generation of a localized Gaussian optical heat source.

**Setup:**

* Loads the NbN geometry and mesh.
* Creates a `GaussianHotspot`.

Parameters:

$$
x_0 = 2.5\,\mu m
$$

$$
y_0 = 2.5\,\mu m
$$

$$
A = 10^8
$$

$$
\sigma = 0.5\,\mu m
$$

**Proves:**

* The Gaussian hotspot generator produces an array matching the mesh dimensions.
* The generated heat source contains positive values.

**Dependencies:**

* Geometry loading
* Mesh generation
* `GaussianHotspot`

**Hard-coded controls:**

* Center `(2.5 μm, 2.5 μm)` — **PHYSICAL TEST CONDITION**
* Amplitude `1e8` — **TEST CONDITION**
* Width `0.5 μm` — **TEST CONDITION**
* Expected shape `100 × 100` — **CURRENT MESH ASSUMPTION**

**Classification:**
Optical-physics unit/integration test.

**Scientific significance:**
Moderate. It verifies that localized optical heating infrastructure exists, but does not establish that the Gaussian has the correct normalization, spatial width, power interpretation, or physical units.

**Assessment:**
Keep. Later strengthen with tests of:

* hotspot center location,
* symmetry,
* peak position,
* width,
* normalization/integrated power.


# SHS Test Suite Summary

## 1. Executive Summary

The SHS test suite effectively validates **software correctness, data structures, and system integration** for the Superconducting Hotspot Simulator framework. It successfully ensures that simulation builds, configuration data flows properly through the pipeline, numerical solvers run to completion without errors, and internal state fields maintain correct dimensions. 

However, the suite currently functions primarily as a **sanity and integration check** rather than a **scientific validation framework**. It lacks rigorous physics benchmarks (such as comparing solver outputs against analytical TDGL solutions, validating physical conservation laws, or verifying current-voltage relationships).

---

## 2. Successes & Strengths

### Robust Pipeline Integration
* **Construction & Mapping:** Successfully proves that geometries, region maps, material properties, and contact masks are correctly assigned and discretized onto the numerical mesh (`test_builder`, `test_material_map`, `test_region_map`, `test_contact_map`).
* **Configuration Flow:** Confirms that parameters loaded from JSON files (e.g., $3.0\text{ K}$ initial temperature, $4.2\text{ K}$ boundaries, $0.001\text{ A}$ current) properly reach the underlying objects.
* **Complex Data State:** Verifies correct initialization of complex-valued Order Parameter ($\psi$) fields and electromagnetic fields (`test_fields`).

### Solid Numerical Foundations
* **Basic Differential Operators:** Confirms core mathematical tools behave correctly on simple fields (gradient of linear fields, zero Laplacian for constant fields).
* **Iterative Solvers:** Validates that both Gauss-Seidel and Red-Black SOR converge below a $10^{-7}$ residual threshold and produce consistent field values within $10^{-4}$ of each other (`test_numerics`, `test_solver_comparison`).

### Execution Sanity
* **Solver Execution:** Demonstrates that the coupled electrothermal solver and standalone electrical solver execute multi-step runs without crashing and generate non-zero Joule heating (`test_coupled_solver`, `test_electrical_solver`).

---

## 3. Limitations & Gaps

### Lack of Physical Validation
* **Arbitrary Physical Assumptions:** Several tests check for non-zero responses rather than physical accuracy. For example, `voltage.max() <= 1.0 V` in `test_electrical_solver` acts as an arbitrary bound rather than a physics-derived check.
* **Unvalidated Initial Conditions:** Initializing $\psi = 1$ at $T = 3\text{ K}$ (`test_fields`) is verified as a software initial state, but is not validated as a physical equilibrium state.
* **No Conservation Checks:** Solvers are not yet evaluated for physical current conservation, thermal evolution accuracy, or quantitative agreement with known physical models.

### Test Architecture & Redundancies
* **Single-Region Hard-coding:** Many mapping tests assume a single-material, single-region geometry (e.g., assuming all cells belong to region `0` or material `0`). These will break or require refactoring once multi-region geometries are introduced.
* **Layer Overlap:** Significant test redundancy exists across boundary loading and builder tests (`test_boundaries`, `test_boundary_database`, `test_simulation_boundaries`, and `test_builder`).
* **Duplicate Assertions:** Redundant assertions exist within individual files (e.g., `test_contact_mask_shape` vs. `test_contact_map_shapes`).
* **Test Pollution:** `test_materials.py` contains top-level execution code (`print(nbn)`) outside of pytest functions, which pollutes test output.

---

## 4. Key Hard-Coded Control Observations

The audit identified two distinct categories of hard-coded values across the suite:

1. **Configuration Candidates (To be consolidated/audited):**
   * **Temperatures:** Initial $3.0\text{ K}$, Boundary $4.2\text{ K}$
   * **Operational:** Current $0.001\text{ A}$, Timesteps ($10^{-9}\text{ s}$, $10^{-6}\text{ s}$), Thermal relaxation rate ($1.0$)
   * **Grid Assumptions:** Fixed $100 \times 100$ mesh assumptions (should be dynamically mesh-derived in unit tests)
2. **Legitimate Test Constants (Keep inside tests):**
   * Small test-grid bounds ($(5,5)$, $(10,10)$ in `test_numerics`)
   * Numerical convergence tolerances ($10^{-7}$ residual, $10^{-4}$ solver agreement)

---

## 5. Next Steps & Recommendations

1. **Keep & Refactor Integration Layer:** Retain current tests for pipeline validation, but remove duplicate shape/boundary assertions.
2. **Clean Up Test Infrastructure:** Remove top-level `print` statements in `test_materials.py` and move benchmark/diagnostic code (`test_solver_comparison.py`) out of the automated unit test path.
3. **Build a Physics Validation Layer:** Add a dedicated test suite focused on comparing solver outputs against analytical solutions and physical conservation principles.