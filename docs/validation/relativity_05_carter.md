# Validation: Solar relativity Phase 2 Kerr constants

## Claim

Solar evaluates `E=-p_t`, `Lz=p_phi`, and the Kerr Carter quantity from a
canonical Boyer–Lindquist state, and can opt into Carter drift monitoring
without coupling the generic geodesic integrator to Kerr.

## Model boundary

The evaluator is valid only for finite states accepted by the subextremal Kerr
BL metric and away from the coordinate axis. It supports null and unit-mass
timelike geodesics. It is a diagnostic invariant evaluator, not a separated
Kerr or Mino-time solver.

## Reference

- Project v3 sections 8.2 and 10 define the Carter expression and the
  max-one drift denominator.
- Black Hole Perturbation Toolkit conventions:
  <https://bhptoolkit.org/conventions.html>.
- Independent maintained Kerr geodesic reference:
  <https://bhptoolkit.org/KerrGeodesics/>.
- CEKG constant-of-motion documentation/source:
  <https://bhptoolkit.org/GremlinEq/doc/a00208.html> and
  <https://bhptoolkit.org/GremlinEq/doc/a00128_source.html>.

## Command

```bash
make tests/relativity/test_kerr_constants \
     tests/relativity/test_geodesics_kerr
./tests/relativity/test_kerr_constants
./tests/relativity/test_geodesics_kerr
```

Both executables were included in the full Release and ASan/UBSan runs.
Verified code commit:
`8a6a2533972685f31dea8eae5f361f948546f285`.

## Inputs

- Literal state: `M=1`, `chi=0.5`, `r=8`, `theta=pi/3`,
  `p=(-1,0.25,3,2)`.
- Generic monitored null ray: `M=1`, `chi=0.7`, initial
  `(r,theta)=(8,1.1)` from a ZAMO local direction `(0.3,0.4,0.5)`.
- Synthetic small invariant, throwing evaluator, non-finite evaluator,
  invalid metric point, unknown geodesic kind, and BL-axis cases.

## Expected

- Literal null `Q=10.270833333333333`; literal timelike
  `Q=10.333333333333333`, each within `2e-14`.
- Generic-ray Hamiltonian, Carter relative drift, and Carter absolute drift
  remain below `1e-10`.
- Small invariants use `max(1,abs(initial))`; evaluator failures terminate
  explicitly before accepting a step.

## Actual

```text
literal/evaluator assertions:     12 passed, 0 failed
Kerr integration assertions:      25 passed, 0 failed
generic initial H:                -1.2490009027033011e-16
generic max Hamiltonian error:     6.1423419174274421e-14
generic max Carter relative error: 1.2883872154118593e-13
generic max Carter absolute error: 3.4887648325820919e-12
accepted/rejected steps:           14 / 0
```

## Error

The generic Carter relative drift used `0.129%` of the `1e-10` gate; the
absolute drift used `3.49%`. E and Lz remained exact in the stationary,
axisymmetric metric test. The literal formulas were within their
double-precision tolerances.

## Result

PASSED on Darwin 23.6.0 arm64 with Apple Clang 16.0.0. Carter monitoring is
accepted as an opt-in Kerr-specific diagnostic.

## Limitations

The generic trajectory is short and not a bound long-duration orbit.
Near-extremal, near-axis, near-horizon, high-momentum, and long secular-drift
sweeps remain open. This phase does not implement separated radial/polar
potentials, Mino time, frequencies, or action-angle variables.

## Fastest falsification

Run `./tests/relativity/test_kerr_constants` and
`./tests/relativity/test_geodesics_kerr`. Changing the `cos^2(theta)` term,
using contravariant instead of covariant momentum, dividing relative drift by
the raw initial value, or swallowing evaluator failures must fail.
