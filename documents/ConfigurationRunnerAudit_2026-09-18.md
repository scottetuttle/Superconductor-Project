# Configuration and runner audit — 2026-09-18

## Scope

This audit covers every JSON file under `configs/`, `tools/config/`, and
`src/shs/config/defaults.json`, plus each Python simulation or diagnostics
runner under `tools/`.

## Configuration authority

The reusable simulation JSON owns the mesh, material, physical timestep,
base duration, initial drive mode, and initial current. A TDGL diagnostics
file now controls only experiment behavior: sampling, a dimensionless duration
multiplier, case sweeps, diagnostics, and output rendering. A
`diagnostics.target_current_A` sweep explicitly selects current drive for that
case. The resolved values and the full effective simulation configuration are
written to each case `summary.json`.

This removed the former duplicate `run.dt`, `run.steps`,
`transport.drive_mode`, and `transport.target_current_A` fields from TDGL tool
configs. The old arrangement silently ignored one of the two declarations. In
particular, `ring_fluxoid_exploration.json` contained a dead 50 nA transport
value while its actual sweep ran at 5 nA.

`run.duration_multiplier` is relative to the simulation's duration, so longer
diagnostic experiments no longer redefine either physical time parameter.
The number of steps is resolved as
`round(simulation.duration * duration_multiplier / simulation.dt)`.

## Result-file safety

All result-producing simulation runners now reserve their output directory
atomically. If the requested directory exists, the runner creates `_1`, `_2`,
and so on. Existing results are never opened for replacement. This applies to:

- full TDGL diagnostics and Josephson plots;
- reduced-vortex diagnostics, including calls routed through the TDGL runner;
- SIS junction and vortex/microwave sweeps;
- Josephson current-reversal protocols;
- coupled convergence and adaptive coupled benchmarks.

`build_results_dashboard.py` intentionally rewrites only
`benchmark_results/index.html`. That file is a generated index of current
results rather than a simulation result, so replacement is expected.

## Documentation contract

Every configuration now starts with `_documentation`, containing:

- `purpose`: what the file is intended to model or test;
- `relies_on`: its upstream simulation/material/geometry dependency;
- `run_command`: the direct command, or an explanation that the file is loaded;
- `authority`: the precedence rule;
- `fields`: an explanation for every configuration leaf, including fields in
  array entries.

Run `python tools/document_configurations.py` after adding or removing a field.
`tests/test_configuration_audit.py` fails when JSON is invalid, a field lacks
documentation, a reusable physics config cannot load, a TDGL experiment
duplicates timestep/drive controls, or output reservation regresses.

## Other findings

`configs/materials/Aluminum.json` and `configs/materials/YBCO.json` were empty,
invalid JSON files with no references. They were removed rather than retained
as misleading selectable materials. Material loading now discards the
documentation object before constructing the material dataclass.

All remaining material, geometry, simulation, and tool references resolve to
existing files. The reduced-vortex routing bug was also fixed: when a numbered
directory is selected, the TDGL wrapper now reports and returns that actual
directory instead of the unnumbered requested path.

