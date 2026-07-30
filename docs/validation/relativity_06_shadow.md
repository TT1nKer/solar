# Validation: Solar relativity Phase 2 Kerr shadow

## Claim

Solar produces the asymptotic Bardeen Kerr critical curve and independently
recovers its two equatorial horizontal edges with CPU backward ray tracing
from a distant ZAMO.

## Model boundary

The analytic curve is the infinite-distance vacuum Kerr benchmark. Numerical
observers are at finite `r=1000M` and `2000M`; their local directions are
solved so the conserved impact parameter is exactly `alpha=-Lz/E`, rather
than approximated by a raw local angle. The two-radius comparison remains an
exterior far-distance convergence check. Backward tracing keeps
future-directed photon momentum and uses a negative affine step. In
Boyer–Lindquist coordinates, `r=r_+ + 1e-3M` is an explicit capture proxy; it
is not a claimed physical horizon crossing. `InvalidMetricPoint` is never
classified as capture.

## Reference

- Project v3 section 13 specifies `xi(r_p)`, `eta(r_p)`, celestial
  `(alpha,beta)`, the Schwarzschild limit, and finite-distance policy.
- Kerr shadow/ZAMO backward-ray reference:
  <https://arxiv.org/abs/1605.08293>.
- Kerr geodesic cross-check package:
  <https://bhptoolkit.org/KerrGeodesics/>.

No external implementation code was copied.

## Command

```bash
make tests/relativity/test_kerr_shadow \
     tests/relativity/test_kerr_shadow_raytrace
./tests/relativity/test_kerr_shadow
./tests/relativity/test_kerr_shadow_raytrace
```

Both executables were included in the full Release and ASan/UBSan runs.
An 80-decimal `mpmath` probe independently evaluated the special radii and
Bardeen endpoints. Verified code commit:
`0f76d411aaf96e8abf31a2a88567c3bf3fedd876`.

## Inputs

- Analytic: Schwarzschild `M=2`; Kerr `M=1`, `chi=+/-0.5`,
  `i=pi/2`, 65 upper-branch samples; inclined `i=pi/3`;
  near-equatorial and near-axis `i=1e-3,1e-14,1e-15`; `M=3` scaling;
  `chi=1e-8` limiting branch; invalid inclinations, sample counts, and
  overflowing Kerr/Schwarzschild mass scaling.
- Numerical: `M=1`, `chi=0.5`, equatorial ZAMOs at `r=1000,2000`;
  separate Schwarzschild `chi=0`, `r=1000`; screen brackets `[-8,0]` and
  `[0,8]`; binary tolerance `2e-7`; negative initial affine step `-0.5`;
  `max_step=2`, `max_affine=4*r_observer`; inner/escape events at
  `r_+ + 1e-3` and `1.1*r_observer`.

## Expected

- Schwarzschild `alpha^2+beta^2=27M^2`.
- Kerr analytic edges
  `[-4.096266658713869, 6.138155724715452]`.
- Kerr sampled screen-distance p95 below `2e-4M` and maximum below `1e-3M`.
- Corresponding edges at `r=1000M` and `2000M` differ by less than `1e-3M`.
- Schwarzschild critical-radius relative error below `1e-6`.
- Screen inputs satisfy `alpha=-Lz/E` within `1e-12`.
- Maximum Hamiltonian and Carter relative errors below `1e-10`.
- Every initialized photon has positive observer frequency.

## Actual

```text
analytic assertions:                 1340 passed, 0 failed
ray-trace assertions:                   9 passed, 0 failed
analytic left/right:                 -4.0962666587138692
                                       6.1381557247154532
Kerr left at r=1000 / r=2000:        -4.0962666869163513
                                      -4.0962666869163513
Kerr right at r=1000 / r=2000:        6.1381557583808899
                                       6.1381557583808899
Kerr sampled p95 / max error:          3.3665437548791033e-08
                                       3.3665437548791033e-08
Schwarzschild left/right:             -5.19615238904953
                                       5.19615238904953
Schwarzschild target:                  5.196152422706632
Schwarzschild maximum relative error:  6.4773122988548604e-09
maximum screen-mapping error:          1.7763568394002505e-15
maximum Hamiltonian error:             9.6027321368194377e-11
maximum Carter relative error:         5.9376738197937386e-29
```

The independent 80-decimal endpoint oracle differed from the double outputs
by `9.89e-16` (left) and `2.00e-15` (right).

## Error

The Kerr maximum sampled distance was `3.36654e-8M`, using `0.0168%` of the
stricter `2e-4M` p95 gate. The Schwarzschild maximum relative root error was
`6.47731e-9`, using `0.648%` of its gate. The Hamiltonian maximum used
`96.0%` of its strict gate. Horizontal-screen-sign, `xi`-denominator-sign,
near-axis residual-division, and invalid-metric-as-capture mutations each
caused the focused tests to fail.

## Result

PASSED on Darwin 23.6.0 arm64 with Apple Clang 16.0.0. The analytic curve and
CPU exterior backward-ray benchmark are accepted for Phase 2.

## Limitations

This is not an image renderer, finite-distance analytic shadow, radiative
transfer model, disk model, or horizon-crossing proof. The numerical check
samples only the two equatorial horizontal edges at two distant radii and one
moderate spin. Near-extremal spins, numerical off-equatorial boundaries, full
two-dimensional boundary convergence, and Kerr–Schild capture remain open.

## Fastest falsification

Run both shadow executables above. Using
`local_phi=-alpha/r_observer` violates the conserved screen mapping; changing
the `xi` denominator sign breaks both analytic edges; restoring
radius-uniform sampling makes the `i=1e-14` interior alpha wrong by about
`9.69e-2`; and classifying `InvalidMetricPoint` as capture breaks the explicit
termination-contract assertion.
