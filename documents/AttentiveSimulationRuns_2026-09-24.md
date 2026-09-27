# Attentive simulation runs

Long optical-transport runs can publish a lightweight spatial preview without
changing the simulated physics. The `live_preview` block controls this mode.
When enabled, the runner renders a four-panel view every configured number of
accepted steps. The panels show electron temperature, order-parameter
amplitude, phase, and current-density magnitude. All detected signed vortices,
the tracked vortex, laser position, and target are overlaid.

The run directory contains:

- `live/latest_state.png`: atomically replaced current preview.
- `live/latest_status.json`: current step, physical time, thermal state,
  progress, laser state, and timing summary.
- `live/progress.html`: automatically refreshing browser page.
- `live/preview_*.png`: bounded recent history. Old previews beyond
  `history_frames` are removed.

The preview is observational. Setting `live_preview.enabled` to `false` skips
all preview plotting and retains the uninterrupted long-run path. Raw
`field_frames`, restart checkpoints, and diagnostic CSV output remain separate.

To request a clean stop from PowerShell, create a marker in the printed run
directory:

```powershell
New-Item benchmark_results/<run-name>/STOP_REQUESTED -ItemType File
```

The runner consumes this marker after the current accepted step, writes its
normal diagnostics and checkpoint, and exits with `user_stop_requested`.

## Runtime accounting

With `profiling.enabled`, `performance.json` accumulates wall time for the
coupled solver, vortex tracking, diagnostic construction and flushes, field
snapshots, previews, and checkpoints. Each section reports total seconds, call
count, mean seconds per call, and fraction of invocation wall time. This is a
low-overhead coarse profiler intended to identify whether larger meshes are
limited by field evolution, compression and disk output, diagnostics, or
visualization. It does not yet split the coupled solver into TDGL, electrical,
two-temperature, and magnetic-screening kernels; that finer instrumentation
should be added after real runs establish that the coupled solver dominates.

Preview cadence should remain much coarser than diagnostic sampling. A value
of 200 steps is the initial interactive default. If `live_preview` occupies a
material fraction of wall time in `performance.json`, increase
`every_steps`, reduce retained history, or disable it.
