# Schwarzschild geodesics

The relativity module is separate from the Newtonian `ForceModel` system. It propagates
test-particle and photon geodesics in a fixed Schwarzschild spacetime; it does not treat
general relativity as an extra three-dimensional force.

## Units and state

The module uses geometric units, `G = c = 1`:

- the mass parameter is `M = G m / c²` in km;
- Schwarzschild `t` and `r`, and the affine parameter, are lengths in km;
- `theta` and `phi` are radians;
- the tangent is `dx^mu / d(lambda)`.

`SchwarzschildSpacetime::from_mass_kg()` performs the SI-mass conversion. Public states
must be outside the Schwarzschild horizon and away from the polar coordinate singularity.

## Equations and diagnostics

`schwarzschild_geodesic_derivative()` evaluates

`d²x^mu/dlambda² + Gamma^mu_ab (dx^a/dlambda)(dx^b/dlambda) = 0`

using the analytic Schwarzschild-coordinate connection. Adaptive DOPRI5 propagation
monitors the metric norm, Killing energy, and axial angular momentum. A horizon guard,
minimum useful step, rejection count, and maximum step count make termination explicit.

## Run

```bash
./solar blackhole circular 10 10 1
./solar blackhole photon 10 1
```

The first command propagates one timelike circular orbit at `r=10M` around a ten-solar-mass
black hole. The second propagates one null circular orbit at the `r=3M` photon sphere.

## Validation boundary

`tests/test_relativity.cpp` checks:

- timelike normalization and a circular orbit at `r=10M`;
- a null circular orbit at the `r=3M` photon sphere;
- radial infall direction and conservation diagnostics;
- termination at the exterior horizon guard;
- conversion from kilograms to geometric mass length;
- rejection of states at the Schwarzschild horizon.

These are analytic/regression checks, not an independent implementation comparison.
Kerr geodesics, weak-field deflection, periapsis advance, ISCO stability scans, and
cross-validation against a trusted relativity code remain future work.
