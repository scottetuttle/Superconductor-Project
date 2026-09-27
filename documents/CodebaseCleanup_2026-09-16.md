# Codebase cleanup — 2026-09-16

## Purpose

Reduce repository noise without removing current experiments, scientific evidence, or user-authored configurations. Generated benchmark output is now ignored so future runs do not make source-control reviews unreadable.

## Removed

- 1,970 generated files (229.15 MiB) from superseded diagnostic, convergence, ring, laser-protocol, and validation runs.
- Empty top-level scaffolding directories and empty result directories.
- The obsolete `coupled_convergence_benchmark_v2.py` and `coupled_convergence_benchmark_v3.py` scripts. The maintained comparison tools are `coupled_convergence_benchmark.py` and `benchmark_adaptive_coupled_updated.py`.
- Empty or unused placeholder modules in optics, vortex physics, coupled solving, TDGL helpers, and visualization.
- The empty legacy `requirements.txt`; dependency declarations live in `pyproject.toml`.
- The superseded adaptive-benchmark working-notes file.

## Preserved evidence

The retained result set contains the user-edited laser-bounce and simple-tweezer experiments, reduced optical-tweezer comparison, latest ring exploration, latest TDGL diagnostics and laser protocol, phase-animation validation, direct-current validation, Meissner validation, and three gauge-equivalence cases.

`benchmark_results/` and `benchmark_results_v2/` are ignored because they are generated artifacts. Existing local results remain available for the dashboard and can be regenerated from the corresponding configurations.

## Maintenance rule

Keep source configurations, compact summaries needed as scientific records, and one representative validation run per capability. Treat plots, GIFs, frame directories, and repeated parameter sweeps as regenerable output. Do not remove a user experiment merely because a newer numbered run exists unless its configuration and scientific purpose are known to be duplicated.

