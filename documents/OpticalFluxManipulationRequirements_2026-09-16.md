# Optical manipulation of single flux quanta: simulator requirements

## Reported mechanism

Veshchunov et al. demonstrated manipulation of Abrikosov vortices in a 90 nm niobium film using a focused laser. The light acts primarily as a moving heat source. It lowers the local superconducting condensation energy and vortex line energy, creating an attractive thermal potential for a vortex. Successful motion requires that this attraction overcome material pinning while viscous vortex drag remains small enough for the vortex to follow the moving temperature profile.

This is analogous to optical tweezers operationally, but it is not a direct radiation-pressure trap. The relevant object is the temperature-dependent free-energy landscape of the superconductor. In the strong-heating regime the spot can instead create a normal region. Flux can be confined by currents and a geometrical barrier at the moving superconducting/normal boundary; cooling and collapse of that region can redistribute or release vortices.

Useful sources:

- I. S. Veshchunov et al., *Optical manipulation of single flux quanta*, Nature Communications 7, 12801 (2016), https://doi.org/10.1038/ncomms12801
- M. S. Scheurer et al., *Optical vortex manipulation for topological quantum computation*, Physical Review B 104, 104501 (2021), https://doi.org/10.1103/PhysRevB.104.104501
- *On-Demand Optical Generation of Single Flux Quanta*, Nano Letters (2020), https://doi.org/10.1021/acs.nanolett.0c02166

The 2021 force-balance treatment uses an overdamped vortex equation with viscous drag, a stationary pinning force, and a moving thermal force. Its Bardeen-Stephen drag coefficient is based on normal resistivity, while the thermal force follows the temperature dependence of vortex energy. This is a useful reduced model for cross-checking full TDGL results.

## Why the current reference run did not drag the vortex

The current route moves approximately 15 nm in 150 fs, a speed of order 100 km/s. The 2016 experiment discusses an estimated upper scale near 1 km/s from a roughly 1 micrometre hotspot and roughly 1 ns thermal response, while demonstrated imaging-compatible and scanner-limited speeds were much slower. The current route therefore outruns formation and translation of a realistic thermal well.

The current Gaussian width is 3 nm, smaller than the 5 nm configured coherence length and only several mesh cells wide. A vortex responds to the free energy of its core and circulating currents over coherence- and penetration-depth scales. A sub-coherence-length temperature perturbation is poorly suited to binding the complete vortex structure and is especially sensitive to grid resolution.

The present model deposits energy into one effective temperature with constant heat capacity and conductivity. It has no electron-phonon delay, substrate spreading, interface resistance, optical penetration, or wavelength-dependent absorption. Consequently the commanded laser coordinate is not yet a calibrated prediction of the actual electronic-temperature minimum that acts on TDGL.

The reference run starts near 14 K for a configured 15.5 K critical temperature and reaches only about 14.11 K. This produces a gradient but does not establish whether a sufficiently deep moving vortex-energy minimum forms. The result is also only 150 fs long, so it mainly probes ultrafast order-parameter response rather than experimentally comparable transport and repinning.

## Required simulation capabilities

### 1. Two-stage comparison model

Implement both of these and require them to agree in their shared validity regime:

1. An inexpensive overdamped point-vortex model containing Bardeen-Stephen viscosity, a temperature-dependent vortex energy, vortex-vortex forces, edge forces, Lorentz force, and the same configured pinning landscape.
2. Full generalized TDGL with the spatial temperature field, scalar potential and charge-conserving normal current, magnetic screening where needed, and identical geometry and material inputs.

The point-vortex model can sweep millions of paths and identify plausible capture windows. Full TDGL should qualify the boundaries of those windows and handle core deformation, annihilation, nucleation, normal regions, and topology changes.

Implementation status: the reduced backend is now available through `tools/tdgl_diagnostics.py` using `physics_backend: reduced_vortex`. Its reference configuration performs capture, 1 m/s transport over 14 micrometres, hold, and laser ramp-down. The screened trial follows within 103 nm and finishes at the destination without invoking its numerical motion limiter. This demonstrates internal force-balance capability, not experimental validation; pinning energies and the prescribed temperature profile remain illustrative.

TDGL integration status: laser waypoints now accept a linearly interpolated `power_fraction`, so the same capture, transport, hold, and release structure can be applied as absorbed heat rather than as an imposed vortex force. The reduced-model force must never be added to the TDGL right-hand side because TDGL already generates thermal, pinning, Lorentz, and vortex-interaction motion from its free-energy and electromagnetic fields.

### 2. Electron-phonon and substrate thermal dynamics

Add electron temperature `Te` and phonon temperature `Tph` with separate heat capacities and conductivities. Deposit optical and Joule power into the appropriate subsystem, couple them through a material-dependent electron-phonon transfer law, and couple phonons to the substrate through an interface/escape term. TDGL must use `Te`, because the condensate responds to the electronic distribution.

Material files must provide measured or cited functions or fitted coefficients for electronic heat capacity, phonon heat capacity, electronic and phonon thermal conductivity, electron-phonon coupling, interface conductance, optical absorption, and substrate properties. Unknown coefficients must not receive plausible-looking universal defaults.

### 3. Realistic optical source

Use incident power, wavelength, reflectivity/absorption, film thickness, penetration depth, and beam profile to calculate absorbed energy. The thermal source may move continuously, but temperature should lag and broaden according to the coupled heat equations. Beam widths should first be tested at several coherence lengths and resolved by multiple cells.

### 4. Calibrated pinning and boundaries

The new Gaussian GL-coefficient defects are useful phenomenology, but their strength must be related to a defect radius and pinning energy or critical depinning current. Include edge/image forces and distinguish point defects, holes, thickness variations, and local `Tc` variations. Otherwise successful motion in a perfectly clean film will overestimate experimental reliability.

### 5. Stochastic activation

Add fluctuation terms consistent with the TDGL normalization and discretization. Near a depinning threshold, success is probabilistic. Report capture and delivery probability across random seeds, not a single deterministic trajectory.

### 6. Capture, following, and release protocol

A valid manipulation trial needs separate phases:

1. Equilibrate a vortex at a known defect with the laser off.
2. Turn on a stationary laser offset from the vortex and test capture.
3. Ramp to constant translation speed and require bounded vortex-laser lag over many thermal times.
4. Stop over a destination defect, ramp the laser down, and require the vortex to remain there.
5. Continue after cooling to detect delayed escape, annihilation, or renucleation.

Record minimum and RMS separation, lag along and transverse to the path, depinning time, maximum temperature, normal-region area, vortex charge, flux per vortex, energy balance, destination retention time, and failures involving unwanted vortices.

## Validation ladder

1. Reproduce a stationary radial temperature profile and its relaxation from an analytical Gaussian heat-diffusion case.
2. Verify the temperature dependence of isolated-vortex free energy and compare the numerical gradient with the reduced thermal-force expression.
3. Measure free vortex mobility under a weak transport current and compare it with Bardeen-Stephen flux-flow drag.
4. Measure the critical force/current for each configured pin and compare it with the reduced pinning model.
5. Map capture versus laser offset and absorbed power for a stationary beam.
6. Map successful transport versus power, width, and velocity, including the no-capture, following, overheating, nucleation, and failed-release regimes.
7. Reproduce a published niobium experiment using its film, substrate, bath temperature, spot size, absorbed power, field, and scan speed before transferring conclusions to NbN nanoscale devices.

## Immediate experiment correction

Before the two-temperature solver is available, the current one-temperature model can perform a qualitative test by increasing the beam width to at least one to several coherence lengths, extending the route into the picosecond-to-nanosecond range, resolving the thermal time with adaptive or implicit stepping, equilibrating the initial pinned vortex, and adding capture/hold/release segments. These runs can expose numerical problems and reject extreme settings, but they cannot yet establish experimental feasibility.
