# SHS Experiment Designer

Start the local editor from the repository root:

```powershell
python tools/serve_config_editor.py
```

Open `http://127.0.0.1:8765` in a browser. Begin from one of the runnable defaults:

- Basic film
- Current-driven vortex
- Optical tweezer
- Optical junction transport
- SNS / resolved weak-link junction
- Perforated ring / SQUID
- Two-temperature hotspot

Loading a default establishes useful dimensions, resolution, physics controls,
and initial objects. You can then draw or edit holes, vortices, and laser paths,
change magnetic field or current, adjust junction width and suppression, or
change feedback-laser power, spot size, target, and follow distance. Controls
that do not apply to the selected experiment are hidden.

The editor downloads four independent JSON files:

- Geometry JSON goes in `configs/geometry`.
- Simulation JSON goes in `configs/simulations`.
- Experiment JSON goes in `tools/config`.
- Project JSON can be imported back into the editor and is not consumed by the simulator.

Keep the three generated filenames together. The exported simulation refers to `<name>_geometry.json`, and the diagnostics file refers to `configs/simulations/<name>_simulation.json`.

The Basic Film, Current Vortex, Junction, Ring, and Hotspot defaults use:

```powershell
python tools/tdgl_diagnostics.py --config tools/config/<name>_experiment.json
```

The Optical Tweezer and Optical Junction defaults use:

```powershell
python tools/optical_transport_overnight.py --config tools/config/<name>_experiment.json
```

For feedback experiments, the first vortex in the object list is controlled by
the laser. Additional vortices evolve freely. The junction editor exports the
implemented phenomenological TDGL weak link; the SNS label is an experimental
starting point and does not imply microscopic Usadel or tunnelling physics.

The editor checks structural concerns such as holes touching boundaries, under-resolved holes or laser spots, invalid laser timing, and missing paths or vortices. Passing these checks does not establish physical validity. Review material parameters, TDGL regime warnings, timestep convergence, thermal timescales, boundary conditions, and magnetic-domain size before interpreting results.

## Physics used by exported experiments

The editor writes ordinary SHS geometry, simulation, and experiment files. It
does not use a separate simplified solver. Runs therefore use the current
coupled implementation, including gauge-covariant TDGL operators, the pyTDGL
normalization, Heun time integration, scalar potential, condensate-depletion
normal conductivity, two-temperature electron-phonon dynamics, and optional
self-consistent magnetic screening.

The editor exposes only the controls needed to assemble common experiments. It
currently writes fixed expert defaults for TDGL parameters, electron-phonon
coefficients, electrical solver tolerances, screening iteration, thermal
boundaries, and output diagnostics. Edit the exported simulation JSON when
those values need material calibration or a convergence study. A generated
configuration is structurally valid input; it is not by itself evidence of
quantitative physical accuracy.

Each diagnostics run records its effective simulation configuration in the
result summary. Compare that record with the exported JSON files when checking
which values reached the solver.
