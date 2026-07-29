# Validation: Solar relativity Phase 0B metrics

## Claim

Solar now has a C++17 fixed-background metric layer for Minkowski Cartesian,
Schwarzschild Boyer–Lindquist, and subextremal Kerr Boyer–Lindquist
spacetimes. The covariant matrices, explicit inverse matrices, analytic
covariant derivatives, inverse-derivative identity, characteristic surfaces,
domain rejection, and real CLI entry point passed the Phase 0B gates.

## Model boundary

This phase evaluates analytic, stationary spacetime metrics in geometrized
units with signature `(-,+,+,+)`. It does not solve the Einstein equations,
evolve matter or spacetime, model backreaction, combine multiple exact
spacetimes, integrate a geodesic, cross a horizon, construct an observer, or
render a black hole. Boyer–Lindquist evaluation is intentionally restricted to
the regular exterior away from the polar coordinate singularity.

The model is therefore a validated foundation for later test-particle and
photon work, not yet a black-hole simulator.

## Architecture

- L0: `metric.h` defines chart identity and the fixed-background metric
  contract over the Phase 0A tensor wrappers.
- L1: the three metric implementations own analytic formulas, inverse
  formulas, derivatives, surfaces, and representable-domain checks.
- L3: `relativity_metric.cpp` converts untrusted argv strings to internal
  values, selects a metric, validates the point, and emits human or JSON
  diagnostics.
- `cli/main.cpp` contains only the existing top-level command diplomacy plus
  one relativity dispatch branch.

Dependency direction remains `CLI -> Metric implementations -> Metric/math
contracts`. No physics rule was added to the CLI.

## Formula audit

The implementation was checked line by line against v3:

- `Sigma = r^2 + a^2 cos(theta)^2`;
- `Delta = r^2 - 2 M r + a^2`;
- `A = (r^2 + a^2)^2 - a^2 Delta sin(theta)^2`;
- all five nonzero Kerr covariant components and five inverse components;
- all listed `r` and `theta` covariant derivatives;
- exact zero `t` and `phi` derivatives from stationarity and axisymmetry;
- inverse derivatives computed only as
  `-g_inverse * partial(g) * g_inverse`;
- stable but algebraically equivalent horizon and stationary-limit formulas;
- strict exterior, axis, finite-value, `Sigma`, `Delta`, `A`, and
  representability checks, without clamping.

Schwarzschild uses the independent `f = 1 - 2M/r` diagonal form. Kerr at
`chi=0` is compared component by component with Schwarzschild. Minkowski uses
the exact diagonal `(-1,1,1,1)` and exact zero derivatives.

## Verified source state and platform

```text
Verified code commit:
a234b160890d127cc546bbecc4d2addb97410a62

macOS 14.8.7 arm64 (Build 23J520)
Apple clang version 16.0.0 (clang-1600.0.26.6)
```

## Commands run

Release configuration:

```bash
make clean
make
make test
./solar relativity metric --metric kerr-bl --M 1 --spin 0.9 \
  --x 0,10,1.5707963267948966,0 --json
git diff --check
```

The CLI JSON was piped through Python's standard `json` parser; parsing
succeeded.

Sanitizer configuration:

```bash
make clean
make CXXFLAGS='-std=c++17 -O1 -g -Wall -Wextra -Iinclude \
  -fsanitize=address,undefined -fno-omit-frame-pointer' \
  solar tests/relativity/test_dual4 tests/relativity/test_kerr_bl \
  tests/relativity/test_math tests/relativity/test_metric_cli \
  tests/relativity/test_metric_derivatives tests/relativity/test_metrics \
  tests/relativity/test_types tests/relativity/test_units
ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 \
UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1 \
  ./tests/relativity/test_dual4
```

The same sanitizer environment was applied individually to all eight
relativity test executables listed below.

After sanitizer testing, a second clean release build and full test run
restored and revalidated the normal `-O2` artifacts.

## Test results

| Executable | Passed | Failed |
|---|---:|---:|
| `test_dual4` | 25 | 0 |
| `test_kerr_bl` | 49 | 0 |
| `test_math` | 15 | 0 |
| `test_metric_cli` | 18 | 0 |
| `test_metric_derivatives` | 3 | 0 |
| `test_metrics` | 70 | 0 |
| `test_types` | 8 | 0 |
| `test_units` | 14 | 0 |
| Relativity subtotal | 202 | 0 |
| Fixture-independent legacy assertions | 56 | 0 |
| Total | 258 | 0 |

The sanitizer run passed all 202 relativity assertions without an
AddressSanitizer or UndefinedBehaviorSanitizer report.

The optional `data/de440.asc` fixture was absent, so its eight assertions were
visibly skipped and are not counted. `test_horizons` remains a print-only
legacy comparison. Seven pre-existing unused-variable warnings remain in
`cr3bp.cpp` and `integrator.cpp`; new Phase 0B sources compile without
warnings.

## Numerical evidence

At `M=1`, `chi=0.9`, `r=10`, `theta=pi/2`, the CLI and literal tests found:

```text
g_tt       = -0.80000000000000004
g_tphi     = -0.17999999999999999
g_rr       =  1.2374706100730106
g_thetatheta = 100
g_phiphi   = 100.97200000000001

g^tt       = -1.2494988244029206
g^tphi     = -0.0022274470981314192
g^rr       =  0.80810000000000004
g^thetatheta = 0.01
g^phiphi   =  0.0098997648805840867
```

Measured maxima:

```text
Kerr sampled inverse identity error:       2.2204460492503131e-16
Kerr analytic derivative vs Dual4:         2.86872e-16
Kerr inverse derivative vs five-point FD:  4.38019e-11
Schwarzschild inverse derivative vs FD:    3.70074e-12
```

All are below their gates of `5e-13` for ordinary inverse identity and `1e-8`
for normalized derivative comparison.

Characteristic-surface literals:

```text
r_plus  = 1.4358898943540672
r_minus = 0.5641101056459328
r_ergo(pi/2) = 2
```

Large finite masses also exercise overflow-resistant surface formulas.
Unrepresentable derived lengths, matrix scales, horizon/axis points,
non-finite inputs, non-positive mass, and `abs(chi)>=1` are rejected.

The derivative test was mutation-checked: negating a production
`partial_r g_tt` sign caused two of three derivative assertions to fail with
normalized errors around `0.443` and `0.669`; restoring the expression returned
the test to green.

## Most likely bugs

1. Near-extremal `abs(chi)` very close to one and points close to the configured
   BL margin are more ill-conditioned than the sampled ordinary exterior
   points; the code rejects exact extremality but has no arbitrary-precision
   reference sweep.
2. The fixed axis floor `abs(sin(theta)) > 1e-12` and epsilon-scaled `Delta`
   floor are numerical policy choices. A later adaptive integrator may need a
   tolerance tied to its state/error scale without weakening the no-clamp
   boundary.
3. Dual4 and five-point checks sample several deterministic points rather than
   a randomized or interval-based domain. A correlated transcription error
   could survive if it also preserves the literals, inverse identity, and
   `chi=0` limit.
4. Extreme but still representable magnitudes can lose substantial relative
   precision before an overflow check rejects them; CPU `double` remains the
   authoritative backend and no multiprecision oracle was run.
5. Only Apple Clang 16 on macOS arm64 was verified. GCC/Linux portability
   remains unverified locally.
6. The real CLI is tested for strict input and JSON structure, but no stable
   versioned output schema has yet been declared for downstream applications.

## Fastest way to falsify

1. Run `./tests/relativity/test_kerr_bl`; changing the sign or factor of
   `g_tphi`, an inverse component, or a horizon formula must fail literals or
   the inverse identity.
2. Run `./tests/relativity/test_metric_derivatives`; an analytic derivative
   mutation must exceed the `1e-8` five-point gate or the tighter Dual4 check.
3. Change Kerr spin to zero and run `test_kerr_bl`; every sampled Kerr matrix
   must still equal independently implemented Schwarzschild.
4. Run the documented Kerr CLI command and parse it as JSON; `metric` must be
   `kerr-bl` and `inverse_error` must stay below `5e-13`.
5. Run the same CLI at `r=1.4358898943540672`, at `theta=0`, with `spin=1`, or
   with `nan`; every command must exit nonzero.
6. Rebuild the relativity tests with AddressSanitizer and
   UndefinedBehaviorSanitizer; any diagnostic invalidates this gate.

## Not completed

- Hamiltonian or separated geodesic integration, dense output, events, and
  constraint projection.
- Observer tetrads, local ray initialization, conserved quantities, or
  timelike proper-time evolution.
- Kerr–Schild coordinates and reliable horizon crossing.
- Radiative transfer, reference rendering, WASM/GPU, Solar Local Patch, UI, or
  movie pipeline.
- DE440's eight external-fixture assertions and GCC/Linux verification.

## Result

Phase 0B passes on the verified source state. The only next allowed
implementation phase is Phase 1, the general Hamiltonian geodesic foundation.
