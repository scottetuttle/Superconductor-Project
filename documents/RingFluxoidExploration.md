# Ring fluxoid exploration

This editable demonstration combines the validated perforated geometry with a trapped `+1` phase winding, `5 nA` left-to-right current, a `0.4 T` perpendicular field, self-consistent magnetic screening, and a `0.8 K` off-center Gaussian hotspot. It is intended to make current crowding, phase winding, local order-parameter suppression, heat flow, and fluxoid balance visible in one compact experiment.

Run it from the repository root:

```powershell
python tools/tdgl_diagnostics.py --config tools/config/ring_fluxoid_exploration.json
python tools/build_results_dashboard.py
```

The editable diagnostic configuration is `tools/config/ring_fluxoid_exploration.json`. Its `_documentation` section explains the controls. The base solver configuration is `configs/simulations/nbn_ring_exploration.json`.

Useful changes include:

- Put several values in `sweeps.diagnostics.applied_Bz_T` or `sweeps.diagnostics.target_current_A` to compare cases.
- Change the seeded winding `charge` from `1` to `-1`, or set `vortices` to an empty list.
- Move the hotspot with `x_fraction` and `y_fraction` or change its strength with `peak_delta_K`.
- Increase `run.steps` from `12` to `50-200` to observe slower relaxation.
- Set `output.gif.enabled` to `false` for faster parameter sweeps.

The reference run completed 12 steps and converged. It recovered `5.000000000006688 nA`, reached `14.6834 K`, retained topological winding `1`, and produced a London fluxoid of `0.9999226 Phi0`. The relative current-continuity defect was `1.34e-15`. The run is physically interesting because the normal transport divides around the hole while the trapped winding creates a circulating supercurrent and the hotspot locally weakens one arm.

This remains a fluxoid-ring experiment. It does not contain two Josephson junctions and must not be interpreted as a quantitative SQUID calculation.
