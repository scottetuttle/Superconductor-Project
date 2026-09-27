> September 14, 2026 update: see [solver corrections and validation](SolverCorrections_2026-09-14.md) for the current implementation and 181-test result. The text below records the earlier development stage.

# SHS Current Development Plan

**Project:** Superconducting Hotspot Simulator (SHS)

**Current Phase:** TDGL Foundation → Gauge-Covariant TDGL

---

# PASS 1 — TDGL Formulation and Validation

**Status: ✅ Substantially Complete**

Completed:

* dimensionless TDGL formulation
* TDGL parameter system
* GL equilibrium model
* complex order parameter
* temperature-dependent equilibrium
* equilibrium phase
* TDGL spatial operators
* initial covariant operators
* TDGL timestep evolution
* uniform analytical transient validation
* normal-state validation
* long-time equilibrium validation
* initial temperature-dependent hotspot response

The remaining work in this pass is primarily strengthening validation rather than building the initial formulation.

---

# PASS 2 — Gauge-Covariant Discretization

**Status: 🚧 Current Development Priority**

The current implementation contains covariant differential operators, but the complete gauge structure has not yet been established.

## Objectives

### 2.1 Review Existing Covariant Operators

Verify:

* discretization
* sign conventions
* grid spacing
* vector-potential scaling
* boundary behavior

---

### 2.2 Link Variables

Introduce link-variable methods where appropriate.

The goal is to preserve gauge consistency directly within the discrete representation.

---

### 2.3 Gauge Transformations

Develop tests demonstrating that physically equivalent gauge choices produce equivalent physical behavior.

Validation should include:

* transformed order parameter
* transformed vector potential
* invariant physical observables

---

### 2.4 Boundary Treatment

Develop physically consistent TDGL boundary conditions.

The implementation must distinguish between:

* mathematical boundary conditions
* superconducting transport boundaries
* insulating boundaries
* electrical contacts

---

# PASS 3 — Superconducting Dynamics

**Status: ⏳ Planned**

Once gauge consistency is established:

* derive supercurrent from the TDGL state
* validate phase-gradient behavior
* validate current-induced order-parameter suppression
* establish transport-current boundary conditions
* investigate phase evolution
* investigate phase-slip behavior

Validation should begin with analytically controlled systems.

---

# PASS 4 — Electromagnetic Self-Consistency

**Status: ⏳ Planned**

Introduce:

* vector potential evolution/solution
* magnetic field calculation
* self-fields
* screening currents
* electromagnetic boundary conditions
* TDGL/electromagnetic coupling

The objective is to move from externally specified `A` toward self-consistent electromagnetic behavior.

---

# PASS 5 — Vortex Physics

**Status: ⏳ Planned**

Develop and validate:

* vortex nucleation
* vortex identification
* vortex position tracking
* vortex motion
* vortex interactions
* pinning
* depinning
* flux flow

Initial validation should use controlled magnetic-field configurations.

---

# PASS 6 — Electrothermal TDGL

**Status: 🚧 Early Foundation**

The first temperature-to-order-parameter connection has already been demonstrated.

Current behavior:

[
T(\mathbf r)\rightarrow\psi(\mathbf r).
]

The eventual goal is:

[
T
\rightarrow
\psi
\rightarrow
J
\rightarrow
Q_J
\rightarrow
T.
]

Development should proceed only after current and superconducting transport are sufficiently established.

Objectives:

* superconducting current calculation
* superconducting/resistive transport
* current redistribution
* Joule heating from the superconducting state
* thermal feedback
* hotspot formation
* hotspot recovery
* optical heating
* comparison with NbN experimental behavior

---

# PASS 7 — Optical Excitation

**Status: ⏳ Planned**

Introduce:

* Gaussian beam profiles
* spatially localized heating
* temporal pulses
* moving excitation
* absorption models
* experimentally motivated optical parameters

This pass will ultimately provide the optical-control mechanism central to the SHS research vision.

---

# Validation Strategy

Every major pass should proceed from simple to complex.

Preferred progression:

```text
Analytical solution
        ↓
Numerical test
        ↓
Controlled physical system
        ↓
Benchmark problem
        ↓
Coupled system
        ↓
Experimental comparison
```

No major physical capability should be considered complete merely because the code executes successfully.

---

# Immediate Next Steps

The immediate development sequence is:

1. Inspect and validate the existing covariant operators.
2. Determine the appropriate link-variable formulation.
3. Establish gauge-transformation tests.
4. Establish physically meaningful TDGL boundary conditions.
5. Develop supercurrent calculations.
6. Validate phase-gradient and current behavior.
7. Only then advance toward self-consistent electromagnetics.

The objective is to make the TDGL foundation scientifically robust before adding additional multiphysics complexity.
