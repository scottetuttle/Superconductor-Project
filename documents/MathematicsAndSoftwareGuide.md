# Mathematics and software guide

This document connects the physical equations to the algorithms and source modules used by SHS. SI units are used at subsystem interfaces. TDGL converts to a matched dimensionless system internally.

## Data flow and source map

| System | Main source |
|---|---|
| Configuration and construction | `src/shs/config/`, `src/shs/config/builder.py` |
| Geometry and active mask | `src/shs/geometry/`, `src/shs/mapping/` |
| Shared fields | `src/shs/physics/fields.py` |
| TDGL model and scaling | `src/shs/tdgl/model.py`, `scaling.py`, `operators.py` |
| TDGL time integration | `src/shs/solvers/tdgl_solver.py` |
| Electrical continuity | `src/shs/solvers/electrical_solver.py` |
| Thermal evolution | `src/shs/solvers/thermal_solver.py` |
| Magnetic response | `src/shs/physics/electromagnetics.py`, `src/shs/solvers/magnetic_solver.py` |
| Coupled fixed point | `src/shs/solvers/coupled_solver.py` |
| Laser and pinning | `src/shs/optics/moving_laser.py`, `src/shs/physics/pinning.py` |
| Fluxoid and diagnostics | `src/shs/physics/fluxoid.py`, `diagnostics.py` |
| Diagnostic runner | `tools/tdgl_diagnostics.py` |
| Lumped SIS/RCSJ physics | `src/shs/physics/josephson_rcsj.py`, `tools/sis_junction_diagnostics.py` |

## TDGL equation

The base normalized equation is

```math
u(\partial_t+i\phi)\psi=(\nabla-i\mathbf A)^2\psi
 + \epsilon(T)\psi-|\psi|^2\psi.
```

The preferred `pytdgl` temperature coefficient is

```math
\epsilon(T)=\operatorname{clip}(T_c/T-1,-1,1),
```

and configured pinning suppression is subtracted from `ε`. The generalized Kramer–Watts–Tobin parameter `gamma` modifies the local amplitude update. With `gamma = 0`, the implementation reduces to the ordinary dissipative TDGL update with a scalar-potential phase rotation.

The normalization uses

```math
\mathbf x'=\mathbf x/\xi,\qquad t'=t/\tau_0,
```

with `τ0 = μ0 σ λ²` in the preferred convention. Electromagnetic scales include

```math
A_0=\frac{\Phi_0}{2\pi\xi},\qquad
B_0=\frac{\Phi_0}{2\pi\xi^2}.
```

The scale object performs every SI-to-dimensionless conversion. Mixing values from `legacy_gl` and `pytdgl` produces inconsistent physics and is prohibited by configuration validation.

### Spatial discretization

The covariant Laplacian uses link variables

```math
U_x=\exp(-iA_x\Delta x),\qquad U_y=\exp(-iA_y\Delta y).
```

A forward x contribution is formed from `Ux ψ(i+1) - ψ(i)`, with the conjugate link used in the reverse direction. This parallel-transports neighboring complex values before subtraction and preserves discrete gauge covariance. Links that cross an inactive cell are omitted, imposing a zero-normal-supercurrent boundary around holes.

Insulating outer boundaries apply the same covariant no-flow condition. Normal contacts impose `ψ = 0` on either a full side or a configured contact mask.

### Time integration

The current production TDGL integrator is explicit. It calculates a conservative stability limit from covariant diffusion, the temperature coefficient, the cubic slope, and the generalized-gamma factor. A requested physical step is divided into internal normalized substeps as necessary. The scalar potential is applied through a temporal link, and the generalized amplitude update solves its local quadratic relation stably.

This method is simple and valuable as a reference, but its diffusion limit scales roughly as `Δx²`. Fine meshes and long physical durations therefore become expensive. The planned improvement is an IMEX method that treats the covariant Laplacian implicitly while retaining error-controlled nonlinear evolution.

## Supercurrent

In normalized form the supercurrent is derived from the gauge-covariant phase gradient, schematically

```math
\mathbf j_s=\operatorname{Im}\{\psi^*(\nabla-i\mathbf A)\psi\}.
```

The TDGL scale converts this to SI current density. Current on inactive links is zero. The `pytdgl` current-density scale includes its required factor of four and is regression tested.

## Electrical continuity

Normal current is

```math
\mathbf J_n=\sigma_n\mathbf E=-\sigma_n\nabla V,
```

and total current is `J = Js + Jn`. Charge conservation gives

```math
\nabla\cdot(\mathbf J_s-\sigma_n\nabla V)=0,
```

or a variable-coefficient Poisson equation for `V`. Face conductivities are constructed on active-active links. No coefficient crosses a hole or insulating boundary.

Voltage mode applies Dirichlet contact potentials. Current mode converts requested terminal current to balanced Neumann flux densities using trapezoidal contact-area integration. Because a pure-Neumann system has an arbitrary constant voltage, the sparse solve adds a mean/reference-voltage constraint.

The preferred backend assembles a sparse matrix and solves it directly with SciPy. Red-black SOR remains available for comparison. Diagnostics recompute finite-volume divergence and integrated source, sink, and cross-section currents rather than trusting the linear-solver residual alone.

Joule heating is

```math
q_J=\mathbf J_n\cdot\mathbf E,
```

because the superconducting current is nondissipative in this model.

## Thermal equation

The one-temperature equation is

```math
C_v\frac{\partial T}{\partial t}
=\nabla\cdot(k\nabla T)+q_J+q_{ext}+q_{laser}
-C_v\Gamma(T-T_{bath}).
```

The solver uses a conservative finite-volume flux balance. Face conductance vanishes across inactive geometry. Fixed-temperature and insulating outer boundaries are supported. The explicit step is subdivided according to a mesh-dependent diffusion stability bound and configured maximum thermal substep.

The conservative form supports meaningful total-energy checks. It does not compensate for missing two-temperature or microscopic energy channels.

## Laser deposition and pinning

For beam center `rL(t)`, width `σL`, absorbed power `P(t)`, and local thickness `d`, the volumetric source is

```math
q_{laser}(\mathbf r,t)=
\frac{P(t)}{2\pi\sigma_L^2d}
\exp\left[-\frac{|\mathbf r-\mathbf r_L(t)|^2}{2\sigma_L^2}\right].
```

Position and power fraction are linearly interpolated between timed waypoints. The coordinate grid and repeated evaluations at the same time are cached. The Gaussian is normalized on an infinite plane; clipped power is not placed back into the device.

Static pinning modifies the TDGL coefficient by

```math
\epsilon_{eff}=\epsilon(T)-
\min\left(s_{max},\sum_k s_k e^{-|\mathbf r-\mathbf r_k|^2/(2\sigma_k^2)}\right).
```

The pinning map is cached until its configuration or mesh changes.

## Magnetic field and screening

Applied uniform perpendicular fields are represented by symmetric or Landau-gauge vector potentials. Magnetic field is recovered as

```math
B_z=\partial_x A_y-\partial_y A_x.
```

For self-field calculations, volume current is converted to sheet current using film thickness. A regularized thin-film Green function is evaluated efficiently using FFT convolution. The solver separates applied, induced, and total vector potential:

```math
\mathbf A=\mathbf A_{applied}+\mathbf A_{induced}.
```

Each screening update mixes the newly calculated induced potential with the previous value using a configured step size. The coupled step accepts only after the relative induced-potential residual and ordinary field residual pass their tolerances. Source stride can reduce cost by sampling current sources, with corresponding loss of resolution.

## Fluxoid

For a closed contour `C`, SHS compares phase winding with the London expression

```math
\Phi_L=\oint_C\mathbf A\cdot d\mathbf l
+\mu_0\oint_C\frac{\lambda^2\mathbf J_s}{|\psi|^2}\cdot d\mathbf l.
```

The topological value is `n Φ0`, where `n` is the wrapped phase winding. The rectangular contour must remain active and maintain `|psi|` above an amplitude floor. A contour crossing a core or normal region is rejected because phase and division by condensate density become ill-defined.

## Reduced vortex solver

The reduced model uses overdamped motion

```math
\eta\mathbf v_i=\mathbf F_{thermal}+\mathbf F_{pin}
+\mathbf F_{interaction}+\mathbf F_{Lorentz}.
```

It uses Bardeen–Stephen drag, a GL thermal-force estimate, screened pair interactions involving the modified Bessel function `K1`, Gaussian pinning energy, and `F_L = q Φ0 J × z`. Explicit position updates include an optional displacement cap. A first-order lag models delayed motion of the thermal center and peak temperature.

This solver is computationally inexpensive because it evolves point coordinates rather than a complex field. Its force formulas and material constants must be calibrated before quantitative prediction.

## Coupled algorithm and rollback

For one physical timestep, the coupled solver starts from the last accepted fields. Each fixed-point pass evaluates the laser at the timestep midpoint, advances TDGL from the accepted order parameter, solves electrical continuity, updates magnetic screening, and advances temperature from the accepted temperature. It measures scaled changes in shared fields, applies relaxation, and uses a convergence controller to identify convergence, slow progress, oscillation, or stalling.

All coupling iterations represent alternative approximations to the same new time level. They must not repeatedly advance physical time. The accepted state replaces the old state only after convergence. Failure preserves the prior fields, which is required for adaptive retry with a smaller timestep.

## Performance characteristics

The dominant costs are explicit TDGL substeps, repeated sparse electrical solves, magnetic FFTs/fixed-point iterations, and image generation. Current optimizations include NumPy vectorization, sparse linear algebra, cached laser coordinates, cached pinning maps, active-link masks, optional magnetic source stride, configurable output sampling, and transactional rollback without rebuilding the complete simulation.

The largest remaining performance improvement is semi-implicit TDGL integration. Sparse topology/factorization reuse and multirate thermal/magnetic updates should follow, guarded by physical error monitors.

## Verification strategy

Software unit tests cover configuration, geometry, operators, boundaries, and solver behavior. Physics qualification additionally checks gauge equivalence, current conservation, heat conservation, Meissner response, vortex winding, fluxoid consistency, hole isolation, and timestep behavior. Diagnostic tools record the effective configuration and raw time series.

A new numerical method is ready for scientific use only after it agrees with analytic limits, converges under mesh and timestep refinement, preserves gauge equivalence and conservation, and matches an independent implementation or experiment in its claimed regime.
