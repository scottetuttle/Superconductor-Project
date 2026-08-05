# SHS AI Development Context

Project: Superconducting Hotspot Simulator (SHS)

Version: 2.0

Purpose: Development Guidelines for AI Contributors

---

# Purpose of This Document

This document provides the context required for an AI assistant to contribute effectively to the SHS project.

It is **not** intended to describe every implementation detail of the codebase. Instead, it explains the project's philosophy, architecture, scientific objectives, and preferred development practices.

An AI should read and understand this document before proposing architectural changes or writing new code.

---

# Project Overview

The Superconducting Hotspot Simulator (SHS) is a research-grade multiphysics simulation framework for superconducting materials and devices.

The long-term objective is to model the interaction between:

* superconductivity
* thermal transport
* electrical transport
* electromagnetics
* optical excitation
* vortex dynamics

using a modular scientific software architecture.

The project is intended for scientific research rather than educational demonstrations or engineering approximations.

---

# Long-Term Scientific Goal

SHS is ultimately intended to become a unified research platform capable of investigating complex superconducting phenomena that are difficult to study analytically or experimentally.

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

Visualization

No module should bypass this structure.

---

# Responsibility of Each Layer

## Configuration

Stores simulation parameters.

Configuration files should never contain implementation logic.

---

## Geometry

Defines physical structure only.

Geometry should never perform physics calculations.

Geometry should not know about materials, temperatures, currents, or numerical solvers.

---

## Mesh

Represents the numerical discretization of geometry.

Mesh generation should remain independent from physics.

---

## Mapping

Maps physical properties onto the numerical mesh.

Current mapping layers include:

* RegionMap
* MaterialMap
* ContactMap

Future mapping systems should extend this concept rather than replacing it.

---

## Fields

Fields represent evolving simulation quantities.

Examples:

* temperature
* voltage
* electric field
* current density
* magnetic field
* vector potential
* superconducting order parameter

Material properties do **not** belong inside Fields.

---

## Physics Modules

Physics modules define governing equations.

Examples:

* thermal
* electrical
* superconductivity
* electromagnetics
* optics
* vortices

Physics modules should not contain numerical algorithms specific to a particular solver.

---

## Numerical Solvers

Solvers evolve the simulation state.

They should use reusable numerical infrastructure whenever possible.

Avoid embedding numerical methods directly inside physics modules.

---

# Current Development Status

The project has completed its foundational infrastructure.

Completed systems include:

* repository architecture
* configuration system
* geometry
* mesh generation
* material database
* region mapping
* material mapping
* contact mapping
* boundary conditions
* simulation builder
* validation layer
* field storage
* thermal solver
* electrical PDE solver
* coupled electrothermal solver
* reusable numerical operators
* iterative solver framework
* Red-Black SOR implementation

The project is now transitioning toward full superconducting physics.

---

# Current Development Priority

The highest development priority is implementing a research-grade Time-Dependent Ginzburg-Landau (TDGL) framework.

The intention is **not** to build a simplified superconductivity model that will later be discarded.

Instead, new infrastructure should directly support the final TDGL implementation whenever practical.

Examples include:

* complex-valued order parameter fields
* gauge-covariant numerical operators
* Ginzburg-Landau material parameters
* dimensionless formulations
* stable nonlinear PDE solvers

Future work should move SHS closer to a complete TDGL simulator.

---

# Numerical Philosophy

Numerical methods should remain independent from physical models.

Reusable numerical infrastructure is strongly preferred.

The same numerical solver should be reusable by multiple physics modules.

Preferred organization:

Numerics

↓

Physics

↓

Solver

↓

Simulation

Avoid implementing specialized numerical routines inside individual physics modules unless absolutely necessary.

---

# Preferred Development Style

Development should always be incremental.

Before modifying existing code:

1. Inspect the current implementation.
2. Understand existing architecture.
3. Extend existing systems.
4. Avoid unnecessary rewrites.

Large architectural changes should only be proposed when they provide significant long-term benefit.

---

# Code Generation Rules

When writing code:

* Preserve existing architecture.
* Avoid duplicate systems.
* Reuse existing infrastructure.
* Maintain compatibility with current tests whenever possible.
* Follow the project's existing coding style.
* Keep implementations modular.
* Separate data, physics, and numerics.

Do not introduce parallel implementations of existing functionality.

---

# Testing Philosophy

Every significant implementation should include corresponding tests.

Whenever practical:

* create tests before major refactoring
* validate numerical correctness
* compare against analytical solutions
* compare against published benchmark problems

Passing tests are considered mandatory before advancing to new development stages.

---

# Performance Philosophy

Correctness takes priority over speed.

Optimization should proceed in the following order:

1. Correct implementation.
2. Validation.
3. Profiling.
4. Algorithmic improvement.
5. Vectorization.
6. Parallelization.
7. GPU acceleration.

Avoid premature optimization.

However, when a substantially better numerical algorithm exists (for example, Red-Black SOR instead of classical Gauss-Seidel), prefer the superior algorithm rather than optimizing an inferior one.

---

# Documentation Philosophy

Whenever architecture changes significantly:

* update CurrentStatus.md
* update Roadmap.md if development priorities change
* update Architecture.md if responsibilities change
* preserve MasterPlan.md unless long-term scientific goals change

Documentation should remain synchronized with the implementation.

---

# Preferred AI Behavior

Before proposing new code:

* understand the current architecture
* inspect relevant files if they are available
* ask for missing files rather than making assumptions
* explain architectural implications before suggesting major changes

When generating code:

* provide complete implementations rather than fragments whenever practical
* avoid placeholder code
* preserve formatting consistency
* maintain compatibility with the existing project structure

Avoid introducing unnecessary abstraction or complexity.

---

# Scientific Philosophy

The objective of SHS is scientific realism rather than the fastest possible implementation.

Whenever a tradeoff exists between convenience and long-term scientific capability, favor the design that better supports future research.

Architectural decisions should be evaluated based on whether they make the eventual TDGL-centered multiphysics simulator more robust, extensible, and scientifically credible.

---

# Final Guiding Principle

Every contribution should move SHS closer to becoming a modular, validated, research-grade superconductivity simulation platform.

Code should not merely solve today's task—it should strengthen the foundation for the physics that will be implemented tomorrow.
