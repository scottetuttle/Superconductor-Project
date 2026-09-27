# Multiply connected TDGL validation checkpoint

## Purpose

This checkpoint adds the minimum geometry and numerical machinery needed to represent a superconducting ring as a true multiply connected domain. It is validation-first: every new behavior is paired with a known geometric, mathematical, or physical result before the ring is used for SQUID or strong-screening claims.

## Implemented capability

- Geometry JSON accepts named circular or rectangular `holes` in SI metres.
- `RegionMap.active_mask` distinguishes superconducting nodes from insulating cutouts.
- The gauge-link covariant Laplacian omits every link that crosses an internal boundary. This imposes zero gauge-covariant normal derivative, and therefore zero normal supercurrent, at an insulating hole.
- The gauge-covariant current operator also removes links crossing the hole. TDGL keeps `psi`, supercurrent, and magnetic source current zero in inactive nodes.
- Fluxoid diagnostics can generate a rectangular superconducting contour around a named hole with a physical clearance.
- `configs/geometry/NbN_ring.json` and `configs/simulations/nbn_ring_validation.json` provide a reproducible reference ring.

## Reference results

Automated tests in `tests/test_multiply_connected_validation.py` establish:

1. The rasterized circular-hole area agrees with `pi*r^2` within 3% on the 101 by 101 reference mesh.
2. A constant order parameter has zero covariant Laplacian everywhere in the active domain, including the internal insulating boundary.
3. Arbitrarily large values placed inside the inactive hole cannot affect any active-domain Laplacian value.
4. The masked covariant operator transforms correctly under a discrete gauge transformation to absolute tolerance `1e-12`.
5. A unit phase winding around the hole gives winding number 1 and a London fluxoid within `3e-4` relative error of one superconducting flux quantum.
6. For a uniform applied field, the discrete vector-potential contour integral agrees with `B` times the enclosed rectangular area to relative tolerance `2e-14`.
7. A real TDGL step preserves exactly zero order parameter and supercurrent inside the hole.

The complete repository suite passes: `236 passed` in `101.22 s`. The pytest cache warning is an environment permission issue and does not affect test execution.

## Validity boundary

The active-domain mask now also removes electrical and thermal face coefficients that cross a hole. The pure-Neumann current solver assembles only active nodes, contacts on inactive nodes are rejected, inactive cells generate no normal current or Joule heat, and thermal evolution excludes inactive cells from diffusion, sources, relaxation, and stability estimates.

`tests/test_perforated_transport_validation.py` verifies symmetric current diversion around the hole, exact zero inactive current, requested-current recovery through every physically integrated cross-section, zero thermal leakage, symmetry of heat flow around the obstacle, and active-domain energy conservation. The quantitative electrical benchmark uses the finite-volume current-controlled path. The legacy voltage-driven Dirichlet path conserves its nodal flux but retains rectangle-rule outer-boundary weighting; it is suitable for symmetry and qualitative voltage checks but should be converted to the same finite-volume convention before being used as an independent quantitative current benchmark.

The next magnetic milestone is a strong-screening equilibrium ring compared with an independent TDGL implementation across mesh, timestep, and screening-tolerance refinements.
