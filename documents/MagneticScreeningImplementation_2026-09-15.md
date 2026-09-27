# Magnetic screening implementation checkpoint

## Implemented equation

The first self-consistent screening path uses the thin-film vector-potential relation

`A_induced(r) = mu0/(4*pi) integral K(r') / |r-r'| dA'`,

where sheet current is `K=J*d`, followed by

`A_total = A_applied + A_induced`.

This matches the screening formulation described by pyTDGL. Because total vector potential changes the TDGL covariant derivatives, which change current and therefore induced vector potential, the update runs inside the existing coupled fixed-point iteration.

Primary references:

- https://py-tdgl.readthedocs.io/en/latest/notebooks/screening.html
- https://github.com/loganbvh/py-tdgl/blob/main/docs/background.rst
- https://www.osti.gov/pages/servlets/purl/1985292

## Numerical implementation

The Cartesian kernel depends only on observer-source displacement, so two FFT convolutions calculate `A_induced,x` and `A_induced,y`. Source values are sheet-current components. The observation and source are collocated in SHS, unlike pyTDGL's edge/site layout, so the singular same-cell kernel is replaced by the area-average `2/a` for a circular cell of equal area, `a=sqrt(dx*dy/pi)`.

Each coupled iteration calculates current from the present total vector potential, computes a target induced vector potential, and applies

`A_induced <- A_induced + alpha*(A_target-A_induced)`.

The physical step is accepted only when the ordinary coupled residual and configured relative screening residual pass. Existing transactional rollback protects the accepted time level when screening fails.

## Current limits

- The same-cell regularization is physically motivated but not yet verified by mesh refinement against an edge/site reference discretization.
- Damped Picard iteration is implemented; pyTDGL's heavy-ball acceleration is not yet implemented.
- There is no fluxoid contour diagnostic yet.
- Gauge-equivalent screened trajectories, Meissner screening, equilibrium vortex fields, and flux quantization are not yet validated.
- The current vector-potential update is magnetoquasistatic and does not include displacement current or electromagnetic radiation.
- Source strides above one are approximate and should not be used for quantitative results.

Screening therefore remains disabled by default. The implementation is ready for the validation and acceleration work needed before using it for SQUID or quantitative vortex predictions.

## Qualification results

The physical fluxoid diagnostic evaluates

`Phi_f = integral A.dl + integral mu0*lambda^2*Js/|psi|^2 dl`

and compares it with `winding*Phi0`. This independent comparison exposed an exact factor-of-four error in conversion of the covariant supercurrent under the pyTDGL normalization. The conversion now accounts for the factor four already present in pyTDGL's declared `J0` scale.

A synthetic `+1` vortex on the 100x100 reference grid now gives a London fluxoid within `3e-4` relative error of `Phi0`. Three complete screened runs using symmetric, `landau_x`, and `landau_y` gauges each converged in 21 iterations and agreed in reported current, mean order parameter, self-field, continuity, and fluxoid. Their London fluxoid was `0.9998091048 Phi0`, while phase winding gave exactly `1 Phi0`.

For a uniform synthetic sheet current, the induced vector potential at the film center changed by `0.895%` between 51x51 and 101x101 grids and by `0.299%` between 101x101 and 151x151 grids. This supports convergence of the same-cell regularization but is not a formal convergence-order result.

The zero-vortex 50 nm NbN reference film at 0.01 T produced only weak Meissner screening: the center field changed from `0.010000000 T` without feedback to `0.009999647 T` with feedback. This is physically expected because the effective thin-film penetration length is much larger than the lateral device size. A larger or more strongly screening geometry is still needed for a quantitative Meissner benchmark.
