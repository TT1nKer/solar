# Validation: Solar relativity Phase 2 local and timelike states

## Claim

Solar initializes future-directed local photons and subluminal unit-mass
timelike states through observer tetrads, and evaluates/constructs analytic
equatorial Kerr circular timelike orbits with existence and stability kept
separate.

## Model boundary

Photon frequency is normalized to one for the chosen observer. Timelike local
energy is the Lorentz factor and is not rescaled. Circular-orbit
`OrbitSense` is relative to the black-hole spin; for negative spin the
coordinate signs of `Omega` and `Lz` reverse. Circular observers are
equatorial BL-exterior constructs. Unstable timelike circular orbits may
exist; the photon radius is rejected as non-timelike.

## Reference

- Project v3 sections 7.5, 8.1, 8.2, and 11 define local state
  initialization, future direction, Hamiltonian targets, and Kerr special
  orbits.
- Maintained independent special-orbit/geodesic reference:
  <https://bhptoolkit.org/KerrGeodesics/>.
- Sign and unit conventions:
  <https://bhptoolkit.org/conventions.html>.

## Command

```bash
make tests/relativity/test_local_initialization \
     tests/relativity/test_kerr_orbits
./tests/relativity/test_local_initialization
./tests/relativity/test_kerr_orbits
```

Both executables were included in the full Release and ASan/UBSan runs.
An 80-decimal independent probe evaluated the `chi=0.5` special radii.
Verified code commit:
`0f76d411aaf96e8abf31a2a88567c3bf3fedd876`.

## Inputs

- Minkowski static observer; photon local direction `(2,0,0)`.
- Minkowski timelike local velocity `(0.6,0,0)`, so `gamma=1.25`.
- Backward null integration with negative affine step.
- Schwarzschild circular orbits at `r=3,5,6`.
- Kerr `M=1`, `chi=+/-0.5` special radii and circular observer at
  `r=8`, `theta=pi/2`.
- Zero/non-finite photon directions, luminal/superluminal velocities,
  invalid tetrads, non-equatorial circular observers, and unknown orbit sense.
- Observer/timelike initialization frames immediately inside and outside the
  strict `1e-12` tetrad gate, plus special-orbit mass-scaling overflow.

## Expected

- Photon measured frequency `1` and null normalized constraint below `1e-14`.
- Timelike measured energy `1.25`, Hamiltonian `-0.5`, and normalized
  constraint below `1e-14`.
- Schwarzschild `r_ISCO=6M`, `r_ph=3M`, `r_mb=4M`.
- Kerr `chi=0.5` radii match analytic literals within `3e-15`.
- Lowering the circular observer velocity reproduces analytic `E` and `Lz`
  within `2e-13`.

## Actual

```text
local-state assertions:             29 passed, 0 failed
special-orbit assertions:           40 passed, 0 failed
photon frequency / constraint:       1 / 0
timelike frequency / constraint:     1.25 / 0
backward-traced photon frequency:    1
chi=.5 prograde/retrograde ISCO:     4.2330025295308262
                                       7.5545847145123579
chi=.5 prograde/retrograde photon:   2.3472963553338606
                                       3.5320888862379562
circular lowered-energy error:       1.1102230246251565e-16
circular lowered-Lz error:           0
circular observer tetrad error:      4.4408920985006262e-16
```

## Error

Against the 80-decimal oracle, the two ISCO errors were
`5.23e-16` and `5.42e-16`; the photon-radius errors were
`9.77e-17` and `1.30e-16`. Local Minkowski constraints and frequency errors
were exactly zero in double precision.

## Result

PASSED on Darwin 23.6.0 arm64 with Apple Clang 16.0.0. Local physical
initialization and analytic equatorial circular-orbit support are accepted for
Phase 2.

## Limitations

No generic circular-orbit finder, inclined/spherical timelike orbit solver,
ISCO perturbation evolution, self-force, backreaction, separated Kerr solver,
or long-term symplectic integration is present. Exact zero errors in the
Minkowski fixture do not imply zero error in curved spacetime.

## Fastest falsification

Run the two focused executables. Removing photon-direction normalization,
using a coordinate-time sign for future direction, accepting speed
`>=1`, treating `r=3M` as timelike, swapping prograde/retrograde signs for
negative spin, or changing the circular normalization must fail.
