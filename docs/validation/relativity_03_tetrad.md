# Validation: Solar relativity Phase 2 observers and tetrads

## Claim

Solar constructs normalized static, arbitrary, look-at, Kerr ZAMO, and
equatorial circular observer frames. Coordinate/local-vector round trips,
orientation, zero-angular-momentum behavior, and explicit failure domains
meet the Phase 2 gates.

## Model boundary

The tetrads use signature `(-,+,+,+)` and store contravariant basis vectors
`e_(a)^mu`. Kerr observers remain in the regular Boyer–Lindquist exterior.
A static observer is rejected where its coordinate worldline is not timelike;
the BL axis and invalid near-horizon points are rejected rather than repaired.
This does not provide Kerr–Schild horizon crossing or an accelerated camera
worldline integrator.

## Reference

- Project v3 sections 7.5 and 10 define local tetrad expansion, observer
  normalization, and ZAMO semantics.
- The backward-ray/ZAMO convention is cross-checked against
  <https://arxiv.org/abs/1605.08293>.
- Metric/sign conventions are cross-checked against
  <https://bhptoolkit.org/conventions.html>.

## Command

```bash
make tests/relativity/test_observers
./tests/relativity/test_observers
make clean && make && make test
```

The same observer test was also run in the documented full ASan/UBSan pass.
Verified code commit:
`0f76d411aaf96e8abf31a2a88567c3bf3fedd876`.

## Inputs

- Minkowski static observer at the origin.
- Minkowski observer with `u=(1.25,0.75,0,0)`, look direction `-z`, and
  up reference `+y`.
- Arbitrary boosted frame built from Cartesian spatial seeds.
- Kerr ZAMO with `M=1`, `chi=0.7`, `r=8`, `theta=1.1`, `phi=0.3`.
- Far-field ZAMO at `r=10^6`, plus ergosphere, horizon-margin, degenerate-seed,
  non-unit-velocity, and polar-axis failures.
- Explicit frame-norm errors near the v3 gate: approximately `5e-13`
  (accepted) and `5e-11` (rejected).

## Expected

- Construction accepts only frames whose maximum over all 16 components of
  `g(e_(a),e_(b))-eta_(a)(b)` is strictly below the v3 `1e-12` gate.
- The ordinary exact-fixture assertions remain within `2e-15`.
- Local/coordinate round trips agree within `2e-15`.
- ZAMO covariant `p_phi` is zero within `2e-17`.
- Invalid, non-timelike, or degenerate constructions fail explicitly.

## Actual

```text
static tetrad error:             0
boosted look-at tetrad error:    0
arbitrary boosted tetrad error:  0
Kerr ZAMO tetrad error:          3.3306690738754696e-16
observer assertions:             44 passed, 0 failed
```

All round-trip, handedness, zero-`p_phi`, far-field, and failure-path
assertions passed.

## Error

The largest reported ordinary-frame orthonormality error was the Kerr ZAMO
value `3.3306690738754696e-16`, about `3.33e-4` of its `1e-12` gate.
No component was clamped or projected to obtain this result.

## Result

PASSED on Darwin 23.6.0 arm64 with Apple Clang 16.0.0. The observer/tetrad
layer is accepted for Phase 2 exterior initialization.

## Limitations

Only deterministic double-precision cases were tested locally. ZAMO and
circular observers are BL-exterior constructs; no claim is made at the axis,
inside the chart-valid boundary, or through the horizon. Long accelerated
observer transport and Fermi–Walker transport remain unimplemented.

## Fastest falsification

Run `./tests/relativity/test_observers`. Negating the ZAMO shift sign,
accepting a static observer inside the ergosphere, changing the tetrad
handedness, or skipping Gram–Schmidt degeneracy checks must fail this
executable.
