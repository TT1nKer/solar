# Validation: Solar relativity Phase 2 Kerr shadow

## Claim

Solar produces the asymptotic Bardeen Kerr critical curve and independently
recovers its two equatorial horizontal edges with CPU backward ray tracing
from a distant ZAMO.

## Model boundary

The analytic curve is the infinite-distance vacuum Kerr benchmark. The
numerical observer is at finite `r=1000M`, so agreement is a convergence test,
not exact equality. Backward tracing keeps future-directed photon momentum
and uses a negative affine step. In Boyer–Lindquist coordinates,
`r=r_+ + 1e-3M` is an explicit capture proxy; it is not a claimed physical
horizon crossing. `InvalidMetricPoint` is never classified as capture.

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
`8a6a2533972685f31dea8eae5f361f948546f285`.

## Inputs

- Analytic: Schwarzschild `M=2`; Kerr `M=1`, `chi=+/-0.5`,
  `i=pi/2`, 65 upper-branch samples; `M=3` scaling; `chi=1e-8`
  limiting branch; invalid inclinations and sample counts.
- Numerical: `M=1`, `chi=0.5`, equatorial ZAMO at `r=1000`;
  screen brackets `[-8,0]` and `[0,8]`; negative initial affine step `-0.5`;
  `max_step=2`, `max_affine=4000`; inner/escape events at
  `r_+ + 1e-3` and `1100`.

## Expected

- Schwarzschild `alpha^2+beta^2=27M^2`.
- Kerr analytic edges
  `[-4.096266658713869, 6.138155724715452]`.
- Numerical edge error below `3e-2`.
- Maximum Hamiltonian and Carter relative errors below `1e-10`.
- Every initialized photon has positive observer frequency.

## Actual

```text
analytic assertions:                 1320 passed, 0 failed
ray-trace assertions:                   7 passed, 0 failed
analytic left/right:                 -4.0962666587138692
                                       6.1381557247154523
numerical left/right at r=1000:      -4.09228515625
                                       6.13232421875
left/right finite-distance errors:    0.0039815024638691909
                                       0.0058315059654523438
maximum Hamiltonian error:            8.329275879277567e-11
maximum Carter relative error:        2.0809373184981404e-29
```

The independent 80-decimal endpoint oracle differed from the double outputs
by `9.89e-16` (left) and `2.00e-15` (right).

## Error

The larger finite-distance edge mismatch was `0.00583151`, or `19.4%` of the
allowed `0.03`. The Hamiltonian maximum used `83.3%` of its strict gate.
Horizontal-screen-sign, `xi`-denominator-sign, and invalid-metric-as-capture
mutations each caused the focused tests to fail.

## Result

PASSED on Darwin 23.6.0 arm64 with Apple Clang 16.0.0. The analytic curve and
CPU exterior backward-ray benchmark are accepted for Phase 2.

## Limitations

This is not an image renderer, finite-distance analytic shadow, radiative
transfer model, disk model, or horizon-crossing proof. The numerical check
samples only the two equatorial horizontal edges at one distant radius and
one moderate spin. Near-extremal spins, arbitrary inclinations, full
two-dimensional boundary convergence, and Kerr–Schild capture remain open.

## Fastest falsification

Run both shadow executables above. Reversing
`local_phi=-alpha/r_observer` mirrors the recovered edges to approximately
`[-6.132,4.092]`; changing the `xi` denominator sign breaks both analytic
edges; classifying `InvalidMetricPoint` as capture breaks the explicit
termination-contract assertion.
