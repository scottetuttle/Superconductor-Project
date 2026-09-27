# TDGL normalization, transport correction, and model strategy

## Implemented convention

The default TDGL convention now follows the generalized-TDGL normalization documented by pyTDGL:

- `tau0 = mu0 sigma lambda^2`
- `Bc2 = Phi0/(2 pi xi^2)`
- `A0 = xi Bc2`
- `J0 = 4 xi Bc2/(mu0 lambda^2)`
- `V0 = xi J0/sigma`
- `epsilon(T) = clip(Tc/T - 1, -1, 1)`

The source equations and units are documented at https://github.com/loganbvh/py-tdgl/blob/main/docs/background.rst.

All scale conversions are selected by `tdgl.normalization`. `pytdgl` is the packaged default. `legacy_gl` preserves the earlier SHS convention so historical comparisons remain reproducible, but mixing scales from the two conventions is prevented by constructing one scale set from the selected model.

## Normal-current correction

The standard model now uses the complete measured normal-state conductivity:

`Jn = sigma_n E`.

The previous `sigma_eff=(1-|psi|^2)sigma_n` expression incorrectly forced quasiparticle conductivity to zero when `|psi|=1`. It remains available only as `electrical.normal_conductivity_model="condensate_depletion"` for legacy phenomenological comparisons. The normal current and supercurrent remain separately recorded, and only normal current contributes to Joule heating.

## Reference regime

`configs/simulations/nbn_tdgl_transport.json` now uses:

- 14 K with `Tc=15.5 K`, giving `T/Tc≈0.903`;
- a 2 nm film, giving `d/xi=0.4` and `d/lambda=0.01`;
- approximately 0.10 coherence lengths per grid spacing;
- the matched pyTDGL normalization and constant normal conductivity.

The diagnostics record explicit warnings when temperature, thickness, or resolution exceed configured validity thresholds. This follows the stated generalized-TDGL assumptions that the dirty-film model is most defensible near `Tc` and that a 2D film should be thin relative to both `xi` and `lambda`.

## Multiple-model architecture

The next architecture should use a small registry of explicit physics models rather than conditionals scattered through solvers. A model profile should declare:

- equations and normalization;
- required material parameters;
- validity conditions;
- supported boundary conditions;
- available numerical methods;
- estimated computational cost; and
- diagnostics required before accepting a result.

Useful initial profiles would be:

1. `normal_ohmic`: scalar electrical and thermal solve with no order parameter.
2. `tdgl_applied_field`: generalized TDGL with prescribed applied field and no screening.
3. `tdgl_screened`: generalized TDGL with induced-vector-potential iterations.
4. `legacy_phenomenology`: explicitly reproduces earlier SHS behavior for comparison.

A deterministic policy should select the cheapest valid profile using geometry, `T/Tc`, field strength, requested observables, estimated screening importance, and recent diagnostic residuals. During a run it can promote the model when anomalies appear, such as growing fluxoid error, vortex nucleation, loss of current uniformity, unexpected heating, or repeated timestep rejection. It should use hysteresis before demoting complexity so the solver does not oscillate between models.

Machine learning is not the right first selector. There is not yet a trusted labeled dataset, and an unconstrained classifier could choose an invalid physical model while appearing fast. ML may later help predict runtime, propose a timestep, rank preconditioners, or flag unusual trajectories. Validity rules, conservation checks, and error estimators should remain authoritative. Every automated choice should record its inputs, selected model, reason, and diagnostics so a human can reconstruct the decision.

## Remaining normalization work

- Implement true current-flux terminal boundary conditions rather than relying on voltage control.
- Add self-consistent induced magnetic screening and fluxoid validation.
- Establish temperature-dependent `xi(T)` and `lambda(T)` if the selected material convention requires them.
- Compare dimensionless and SI observables against an external generalized-TDGL reference case.
- Calibrate NbN material inputs from a single internally consistent experimental dataset.
