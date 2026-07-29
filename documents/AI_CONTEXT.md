# SHS AI Development Context

## Project

Superconducting Hotspot Simulator (SHS) is a modular scientific simulation framework for studying superconducting devices controlled through localized heating.

The project focuses on coupling:

- thermal physics
- superconductivity
- optical heating
- electrical transport
- electromagnetic effects
- vortex dynamics

The goal is to build a research-grade simulation platform.

---

## AI Rules

When assisting with SHS:

- Inspect existing files before suggesting changes.
- Do not assume file contents.
- Make incremental changes.
- Preserve existing tests.
- Do not implement future features prematurely.
- Separate configuration, physics, solvers, and visualization.
- Use existing architecture instead of creating parallel systems.
- Give exact implementation steps for the current task.

---

## Architecture Principles

SHS follows:

Configuration
→ Data Objects
→ Geometry/Mesh
→ Fields
→ Physics Models
→ Solvers
→ Visualization

Physics modules should not parse configuration files.

Geometry describes physical layout only.

Materials describe physical properties.

Solvers update fields.

---

## Current Architecture Summary

Main modules:

geometry/
- device layout
- regions
- contacts
- mesh generation

materials/
- material definitions
- database

physics/
- thermal
- superconductivity
- electromagnetics
- vortices

solvers/
- numerical evolution

optics/
- laser and hotspot models

---

## Current Development Stage

The foundation layer is being built.

Completed:

- materials system
- geometry system
- JSON configuration
- mesh generation
- fields
- thermal solver

Current priority:

Complete geometry → mesh → material mapping before implementing advanced superconducting physics.

---

## Documentation

For details consult:

MasterPlan.md:
Long-term vision

Architecture.md:
Software design

CurrentStatus.md:
Current implementation state

Roadmap.md:
Development sequence

References.md:
Scientific sources