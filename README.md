# Superconducting Hotspot Simulator (SHS)

SHS is a Python framework for coupled thermal, electrical, and simplified
Time-Dependent Ginzburg–Landau simulations on rectangular superconducting films.

## Human documentation

- Start with [`documents/UserGuide.md`](documents/UserGuide.md) for the simulator workflow, confidence levels, and result checklist.
- Read [`documents/PhysicsGuide.md`](documents/PhysicsGuide.md) for the physical meaning and limitations of every implemented subsystem.
- Read [`documents/MathematicsAndSoftwareGuide.md`](documents/MathematicsAndSoftwareGuide.md) for equations, numerical methods, coupling, performance, and source-code locations.

Start with [AI context](documents/AI_CONTEXT.md), then the
[September 2026 solver corrections and validation](documents/SolverCorrections_2026-09-14.md).
Older status documents describe earlier development stages.

## Visual experiment editor

Create rectangular devices, place circular or rectangular holes, seed signed vortices, and draw timed laser paths in the local browser editor:

```powershell
python tools/serve_config_editor.py
```

Then open `http://127.0.0.1:8765`. Export geometry to `configs/geometry`, simulation settings to `configs/simulations`, and diagnostics to `tools/config`. See [the editor guide](tools/config_editor/README.md) for details and physical-validation limits.

## Installation

From the repository root, with Python 3.10 or newer:

```powershell
python -m pip install -e ".[test,benchmark]"
```

## Small configured run

Run from the repository root because material and geometry paths are currently
resolved relative to `configs/`.

```python
from shs.config.builder import build_simulation
from shs.solvers.coupled_solver import run_coupled_simulation

simulation = build_simulation("configs/simulations/nbn_hotspot_test.json")
result = run_coupled_simulation(simulation, steps=2, dt=1e-14)
assert result.converged
print(result.elapsed_time, result.total_coupling_iterations)
```

Public timesteps are seconds. Shared vector potential and current are SI;
TDGL normalizes internally. External volumetric heating is assigned to
`simulation.fields.external_heat_source` in W/m^3. `heat_source` is the total
external plus Joule heating. The full duration in the example JSON is much more
expensive than this small run.

## Tests

```powershell
python -m pytest -q -p no:cacheprovider
```

The September 14 correction pass completed with **181 passing tests**.
See the validation report for the original failures and interpretation limits.

## Structure

- `configs/`: simulation, material, and geometry JSON.
- `src/shs/`: construction, mesh mappings, shared fields, physics, numerics, solvers.
- `tests/`: infrastructure, numerical, TDGL, coupled, and regression checks.
- `tools/`: experimental convergence and adaptation benchmarks.
- `documents/`: architecture, plans, scientific validation, and development history.
- `benchmark_results*/`: saved experiment outputs.

Generalized TDGL, vortices, perforated films, and an early thin-film magnetic
screening loop are implemented. Josephson weak links, electron–phonon
nonequilibrium, strong-screening qualification, and quantitative experimental
calibration remain future work. Passing tests does not establish that the
current model predicts experimental hotspot or optical-manipulation thresholds.
