# SHS AI Development Context

**Project:** Superconducting Hotspot Simulator (SHS)
**Version:** 2.1
**Purpose:** Development Guidelines for AI Contributors

---

# Purpose of This Document

This document provides the context required for an AI assistant to contribute effectively to the SHS project.

It describes the project's scientific objectives, architectural philosophy, current physics capabilities, numerical approach, development state, and preferred development practices.

It is not intended to describe every implementation detail of the codebase. The source code and tests remain the authoritative representation of implementation.

An AI should read and understand this document before proposing architectural changes or writing new code.

---

# Project Overview

The Superconducting Hotspot Simulator (SHS) is a research-oriented multiphysics simulation framework for superconducting materials and devices.

The long-term objective is to model the interaction between:

* superconductivity
* thermal transport
* electrical transport
* electromagnetics
* optical excitation
* vortex dynamics
* nonequilibrium effects

using a modular scientific software architecture.

The project is intended to develop toward research-grade simulation rather than educational demonstrations or purely engineering approximations.

---

# Long-Term Scientific Goal

SHS is intended to become a unified computational platform capable of investigating superconducting phenomena that are difficult to study analytically or experimentally.

Major research interests include:

* electrothermal hotspot dynamics
* current redistribution
* nonequilibrium superconductivity
* Time-Dependent Ginzburg-Landau (TDGL) physics
* vortex dynamics
* Josephson devices
* SQUIDs
* superconducting electronics
* optical control of superconductivity
* phase-slip phenomena
* Higgs-mode and amplitude dynamics
* advanced superconducting device concepts

Every architectural decision should support these long-term objectives.

---

# Central Design Philosophy

The architecture of SHS is intentionally modular.

Different physical systems should remain independent whenever possible.

The preferred dependency chain is:

Configuration

↓

Geometry

↓

Mesh

↓

Mappings

↓

Simulation State

↓

Physics Models

↓

Numerical Solvers

↓

Analysis

↓

Visualization

No module should bypass this structure without a strong scientific or architectural reason.

---

# Responsibility of Each Layer

## Configuration

Stores simulation parameters.

Configuration files should not contain implementation logic.

---

## Geometry

Defines physical structure only.

Geometry should not perform physics calculations.

Geometry should not know about temperatures, currents, superconducting states, or numerical solvers.

---

## Mesh

Represents the numerical discretization of geometry.

Mesh generation remains independent from physical models.

---

## Mapping

Maps physical properties onto the numerical mesh.

Current mapping systems include:

* RegionMap
* MaterialMap
* ContactMap

Future mapping systems should extend this concept rather than replacing it.

---

## Fields

Fields represent evolving simulation quantities.

Current and planned examples include:

* temperature
* voltage
* electric field
* current density
* magnetic field
* vector potential
* superconducting order parameter `psi`

Material properties do not belong inside Fields.

---

## Physics Modules

Physics modules define governing equations.

Current and planned modules include:

* thermal
* electrical
* superconductivity / TDGL
* electromagnetics
* optics
* vortices

Physics modules should define physical behavior rather than embedding solver-specific numerical algorithms.

---

## Numerical Solvers

Solvers evolve the simulation state.

Reusable numerical infrastructure is preferred.

Numerical methods should remain as independent from physical models as practical.

---

# Current Scientific State

SHS has completed its foundational simulation infrastructure and its initial coupled electrothermal system.

The project has now entered active superconducting physics development.

The first major TDGL foundation has been implemented and validated.

The current TDGL framework is based on a dimensionless formulation of the form:

[
u\frac{\partial\psi}{\partial t}
================================

D^2\psi
+
\left(1-\frac{T}{T_c}\right)\psi
--------------------------------

|\psi|^2\psi
]

where the gauge-covariant derivative is represented conceptually as:

[
D=\nabla-i\mathbf A.
]

The order parameter is complex-valued.

The current implementation supports temperature-dependent equilibrium superconductivity and time evolution of the order parameter.

---

# Validated TDGL Capabilities

The following capabilities have been implemented and tested.

## Complex Order Parameter

SHS now has a complex superconducting order parameter field:

```text
psi
```

with amplitude and phase represented naturally by its complex value.

---

## Equilibrium Superconducting State

For reduced temperature

[
t=\frac{T}{T_c},
]

the homogeneous equilibrium amplitude is:

[
|\psi|_{\mathrm{eq}}
====================

\sqrt{1-t}
]

for

[
t<1,
]

and:

[
|\psi|_{\mathrm{eq}}=0
]

for:

[
t\geq1.
]

This behavior has been explicitly tested at zero temperature, intermediate temperature, near-(T_c) temperature, (T_c), and above (T_c).

---

## Equilibrium Order-Parameter Phase

The TDGL model can construct equilibrium complex order parameters with specified phase.

The amplitude and phase are treated independently through the complex representation.

---

## Spatial Differential Operators

The TDGL numerical infrastructure includes:

* gradient
* Laplacian
* covariant gradient
* covariant Laplacian

Basic analytical cases have been tested.

---

## Gauge-Covariant Operators

The current implementation recognizes the vector potential through:

[
D=\nabla-i\mathbf A.
]

Tests verify expected behavior for constant fields and constant vector potential.

This represents the initial gauge-aware numerical foundation.

It does not yet constitute a complete gauge-invariant electromagnetic implementation.

---

## TDGL Time Evolution

The TDGL solver can evolve the superconducting order parameter in time.

The implementation has been tested for:

* finite, stable evolution
* superconducting-state relaxation
* decay above (T_c)
* convergence toward equilibrium

---

## Analytical Uniform TDGL Validation

For the spatially uniform case, the numerical solution has been compared against the analytical solution of:

[
u\frac{d\psi}{dt}
=================

a\psi-\psi^3,
]

with:

[
a=1-\frac{T}{T_c}.
]

The numerical evolution reproduces the analytical transient within the tested tolerance.

This is an important validation of the TDGL time integration rather than merely a software correctness test.

---

# Temperature-Dependent Superconductivity

The TDGL model responds to the local temperature field.

A localized region raised to approximately:

[
0.95T_c
]

produces local suppression of the superconducting order parameter relative to the surrounding colder material.

This demonstrates the first spatial connection between the thermal field and superconducting order parameter.

The current test establishes qualitative hotspot suppression.

It does not yet establish a fully self-consistent electrothermal TDGL system.

---

# Current Foundational Capabilities

The underlying SHS infrastructure includes:

* repository architecture
* configuration system
* geometry
* mesh generation
* material database
* region mapping
* material mapping
* contact mapping
* boundary-condition infrastructure
* simulation construction
* simulation validation
* shared field storage
* thermal physics
* thermal solver
* electrical transport
* electrical PDE solver
* Joule heating
* coupled electrothermal solver
* reusable differential operators
* iterative solver framework
* Red-Black SOR
* TDGL parameters
* GL equilibrium model
* complex order parameter
* gauge-covariant operators
* TDGL time integration
* initial temperature-dependent TDGL coupling

---

# Current Scientific Limitations

The following should NOT be considered complete:

* full gauge-invariant discretization
* link-variable formulation
* gauge-invariance validation
* physically complete TDGL boundary conditions
* self-consistent electromagnetic evolution
* magnetic screening
* current-induced superconducting dynamics
* vortex nucleation
* vortex motion
* vortex pinning
* phase-slip dynamics
* Josephson physics
* complete electrothermal-TDGL feedback
* optical excitation
* experimental NbN validation

These are future development stages.

---

# Numerical Philosophy

Numerical methods should remain independent from physical models.

Reusable numerical infrastructure is strongly preferred.

The preferred structure is:

Numerics

↓

Physics

↓

Solver

↓

Simulation

Avoid implementing specialized numerical routines inside individual physics modules unless scientifically necessary.

Correctness takes priority over performance.

Optimization should proceed approximately in this order:

1. Correct implementation.
2. Physical validation.
3. Profiling.
4. Algorithmic improvement.
5. Vectorization.
6. Parallelization.
7. GPU acceleration.

---

# Testing Philosophy

Every significant physical or numerical implementation should have corresponding tests.

Whenever practical:

* test mathematical identities
* test limiting cases
* test analytical solutions
* test equilibrium states
* test stability
* test convergence
* test physical transitions
* compare against published benchmarks

Passing software tests is necessary but not sufficient for scientific validation.

A successful test should be interpreted according to exactly what physical or numerical property it establishes.

---

# Preferred Development Style

Development should remain incremental.

Before modifying existing code:

1. Inspect the current implementation.
2. Understand the existing architecture.
3. Identify the relevant physics.
4. Extend existing systems where appropriate.
5. Add tests.
6. Validate the result.
7. Update project documentation at the milestone.

Avoid unnecessary rewrites.

Avoid introducing parallel implementations of existing functionality.

---

# Documentation and Continuity

SHS uses several levels of project documentation.

## AI_CONTEXT.md

Describes the project's architecture, scientific philosophy, current capabilities, limitations, and preferred AI development behavior.

It should change relatively infrequently.

---

## CurrentStatus.md

Describes the actual current implementation state.

It should be updated after significant milestones.

---

## CurrentPlan.md

Describes the immediate development sequence.

It should be updated when a development pass is completed or priorities change.

---

## DevelopmentLog.md

Records completed development milestones chronologically.

This is the project's continuity record.

DevelopmentLog should be updated frequently enough that an interrupted development session can be reconstructed without relying on conversation history.

This document is particularly important for maintaining continuity when external development assistance becomes temporarily unavailable.

---

# Preferred AI Behavior

Before proposing new code:

* understand the current architecture
* inspect relevant existing files
* distinguish implemented capabilities from planned capabilities
* identify scientific assumptions
* ask for missing information rather than inventing implementation details

When generating code:

* provide complete implementations whenever practical
* preserve project architecture
* reuse existing infrastructure
* maintain compatibility with existing tests
* avoid unnecessary abstraction
* separate physical equations from numerical algorithms

When interpreting tests:

* state exactly what the test proves
* do not claim stronger validation than the test establishes
* distinguish software correctness from scientific validation

The AI should prioritize scientific correctness over simply making tests pass.

---

# Final Guiding Principle

Every contribution should move SHS closer to becoming a modular, validated, research-grade superconductivity simulation platform.

Code should not merely solve today's task.

It should strengthen the foundation for the physics implemented tomorrow.

The eventual objective is not merely a functioning program.

The objective is a scientifically credible computational laboratory for exploring superconducting physics.
