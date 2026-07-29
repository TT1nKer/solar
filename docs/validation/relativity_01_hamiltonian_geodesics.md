# Validation: Solar relativity Phase 1 Hamiltonian geodesics

## Claim

Solar now has a general eight-dimensional canonical Hamiltonian geodesic
foundation for fixed-background metrics. The implementation provides the v3
Hamilton equations, a fixed-size adaptive Dormand–Prince 5(4) solver with
component tolerances and fourth-order dense output, directed bracketed events,
explicit termination diagnostics, and independently checked Minkowski,
Schwarzschild, and Kerr Boyer–Lindquist trajectories.

This is the CPU reference foundation for later observer and black-hole ray
work. It is not yet a complete black-hole simulator or renderer.

## Model boundary

The state is ordered
`(x^0,x^1,x^2,x^3,p_0,p_1,p_2,p_3)`. Coordinates are contravariant,
canonical momenta are covariant, the signature is `(-,+,+,+)`, and units are
geometrized with `G=c=1`. The background metric is fixed and the test particle
or photon does not backreact.

Phase 1 accepts already constructed coordinate-basis canonical states. v3
section 7.5 local initialization is intentionally deferred because `Tetrad`,
observer normalization, measured frequency, and local round-trip contracts
belong to Phase 2. This phase therefore does not claim a general
observer-independent future-direction test.

Kerr and Schwarzschild integration in this phase uses the regular
Boyer–Lindquist exterior only. An invalid BL trial is rejected and shrunk; it
is never called a physical horizon crossing. There is no Kerr–Schild chart,
interior evolution, singularity treatment, matter, radiation, rendering,
constraint projection, Carter implementation, or separated Kerr solver.

## Architecture

- L0: `geodesic_types.h` owns event direction, every v3 termination reason
  including `InteriorCutoff`, diagnostics, and event payloads.
- L1: `dopri5.h` owns the fixed-size solver contract; internal detail headers
  separately own the seven-stage kernel and continuous-extension matrix.
- L1: `hamiltonian.cpp` owns only canonical contractions and state packing.
- L1: `event_root.cpp` owns one safeguarded secant/bisection bracket.
- Internal L1 adapters classify one metric-aware DOPRI trial and arbitrate
  events without depending on the L2 integration flow.
- L2: `geodesic_integrator.cpp` owns one accepted/rejected-step state machine,
  limits, diagnostics, invariant opt-in policy, and truthful termination.
- L3 compatibility: the existing dynamic-vector
  `dopri5_generic_step` adapts to the same stage engine with its established
  scalar tolerance and maximum norm. The specialized N-body solver remains
  unchanged.

Dependency direction is
`Geodesic flow -> trial/event adapters -> Hamilton/DOPRI/event algorithms ->
contracts`. No CLI or UI business rule enters the physics layers.

## Formula audit

The implementation was checked line by line against v3:

```text
H       = 1/2 g^mu,nu p_mu p_nu
dx^mu   = g^mu,nu p_nu
dp_mu   = -1/2 partial_mu(g^alpha,beta) p_alpha p_beta
```

The normalized constraint is exactly:

```text
abs(H-H0) /
(1 + abs(H0) + 1/2 sum_mu,nu abs(g^mu,nu p_mu p_nu))
```

The null target is `0`; the unit-mass timelike target is `-1/2`.
The RHS uses all four metric derivative dimensions. It does not force
`p_t` or `p_phi` to be conserved. Existing stationary and axisymmetric
metrics produce exact zero derivatives naturally, while a test-only
time-dependent metric produces `dp_t=+0.05`.

DOPRI5 uses the standard seven-stage 5(4) tableau, the fifth-order trial,
embedded error coefficients, v3 per-component scale, a numerically stable RMS
norm, sign-preserving controller, and Shampine fourth-order continuous
extension. Event roots require a sign bracket and use safeguarded
secant/bisection with a 100-iteration cap. No default momentum projection is
present.

## Verified source state and platform

```text
Verified code commit:
57b65256094e6d6b41c32ffd2906310b0ed4b134

Darwin 23.6.0 arm64
Apple clang version 16.0.0 (clang-1600.0.26.6)
```

GCC/Linux and other architectures were not locally verified.

## Commands run

Release configuration:

```bash
make clean
make
make test
git diff --check
```

`make test` invoked all discovered test executables, including all 15
`tests/relativity/test_*` binaries.

Sanitizer configuration:

```bash
make clean
make CXXFLAGS='-std=c++17 -O1 -g -Wall -Wextra -Iinclude \
  -fsanitize=address,undefined -fno-omit-frame-pointer' \
  $(find tests/relativity -name 'test_*.cpp' -type f | sort | \
    sed 's/\.cpp$//')
```

Every relativity executable was then run individually with:

```bash
ASAN_OPTIONS=detect_leaks=0:halt_on_error=1
UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1
```

A second `make clean`, normal `make`, and `make test` restored and revalidated
the release artifacts after sanitizer testing.

## Test results

| Executable | Passed | Failed |
|---|---:|---:|
| `test_dopri5` | 49 | 0 |
| `test_dual4` | 25 | 0 |
| `test_geodesic_events` | 28 | 0 |
| `test_geodesic_failures` | 10 | 0 |
| `test_geodesics` | 49 | 0 |
| `test_geodesics_kerr` | 11 | 0 |
| `test_geodesics_schwarzschild` | 11 | 0 |
| `test_hamiltonian` | 26 | 0 |
| `test_kerr_bl` | 49 | 0 |
| `test_math` | 15 | 0 |
| `test_metric_cli` | 18 | 0 |
| `test_metric_derivatives` | 3 | 0 |
| `test_metrics` | 70 | 0 |
| `test_types` | 8 | 0 |
| `test_units` | 14 | 0 |
| Relativity subtotal | 386 | 0 |
| Fixture-independent legacy assertions | 66 | 0 |
| Total | 452 | 0 |

All 386 relativity assertions also passed under AddressSanitizer and
UndefinedBehaviorSanitizer with no runtime diagnostic.

The optional `data/de440.asc` fixture was absent, so its eight external-data
assertions were visibly skipped and are not counted. `test_horizons` remains a
successful print-only comparison. The release and sanitizer builds each
reported one pre-existing warning in `src/cr3bp.cpp` for `sub_iter`; new
Phase 1 sources compile without warnings.

## Numerical evidence

DOPRI5 on `y'=y`:

```text
normalized embedded error for h=0.01:  7.2527173583230749e-3
dense midpoint absolute error:          6.5875835675299754e-8
fixed global errors h=.2,.1,.05:         1.6884248621451547e-7
                                         6.3380456438721922e-9
                                         2.1638779656996121e-10
observed orders:                         4.7354961186746989
                                         4.8723470052423368
```

Directed dense roots:

```text
positive-step affine error:  5.5511151231257827e-17
negative-step affine error:  1.1102230246251565e-16
```

Minkowski:

```text
null maximum normalized constraint: 0
forward/backward maximum state error: 4.4408920985006202e-17
```

Schwarzschild, `M=1`:

```text
outgoing radial-null maximum constraint: 3.8555017764255443e-15
weak-field b=100 deflection:              0.041222539797915125
leading 4M/b target:                      0.040000000000000001
relative difference:                      3.06%
weak-field maximum constraint:            3.5239255270923555e-12
```

The same test checks exact photon-sphere radial RHS conditions at `r=3M`.

Kerr BL, `M=1`, `chi=0.7`:

```text
initial H:                    -1.2490009027033011e-16
maximum normalized constraint: 6.1423419174274421e-14
accepted/rejected steps:       14 / 0
monitored E drift:             0
monitored Lz drift:            0
Carter diagnostic:             NaN (unavailable by design)
```

All ordinary-ray constraint maxima are below the v3 `1e-10` gate.

## Failure and mutation evidence

- Negating a nonzero dense-extension coefficient caused both midpoint and
  endpoint assertions to fail for positive and negative steps.
- Removing the Hamilton RHS `1/2` factor caused the time-dependent-metric
  `dp_t` assertion to fail.
- Replacing the RMS divisor by a non-RMS expression caused the
  two-component RMS assertion to fail.
- Reclassifying an invalid metric-domain trial as a horizon caused
  `test_geodesic_failures` to exit nonzero; restoring
  `InvalidMetricPoint` returned it to 10/10.
- Extremely tight but finite tolerances exercise overflow-resistant RMS
  accumulation. Unknown norm enums and overflowing constraint denominators
  are rejected rather than silently accepted.

## Most likely bugs

1. Endpoint-bracket event detection can miss an even number of roots, a
   tangential root, or multiple closely spaced roots inside one accepted
   step. Event-dense trajectories still need convergence sweeps with reduced
   `max_step`.
2. Boyer–Lindquist remains singular at the horizon. `InvalidMetricPoint` is a
   chart-domain outcome, not evidence of capture; physical crossing requires
   Phase 3 Kerr–Schild coordinates and an explicit horizon event.
3. Long bound timelike trajectories use adaptive DOPRI5, not a
   structure-preserving method. Secular invariant drift beyond the short
   ordinary tests is not characterized.
4. Coordinate canonical states are caller-supplied. Until Phase 2 tetrads and
   observer-frequency normalization exist, a caller can provide a null state
   with physically unintended local direction even if the Hamiltonian
   constraint is valid.
5. CPU `double` and ordinary summation remain authoritative. Near-extremal,
   near-margin, very large-momentum, and long-duration cases lack a
   multiprecision oracle.
6. E and Lz monitoring is opt-in and meaningful only when the caller knows
   the metric has the corresponding symmetries. Carter remains unavailable.
7. Only Apple Clang 16 on macOS arm64 was verified locally; GCC/Linux and
   other floating-point environments may expose portability or threshold
   differences.

## Fastest way to falsify

1. `./tests/relativity/test_dopri5`: wrong tableau, RMS, component scale,
   controller sign, dense coefficient, or non-finite handling must fail.
2. `./tests/relativity/test_hamiltonian`: wrong variance/order, target,
   denominator, derivative dimension, sign, or `1/2` must fail.
3. `./tests/relativity/test_geodesic_events` and
   `./tests/relativity/test_geodesic_failures`: direction, bracket, iteration,
   trial-domain, and root-failure semantics must remain explicit.
4. `./tests/relativity/test_geodesics`: analytic Minkowski lines, limits,
   first-event selection, and reversal must pass.
5. `./tests/relativity/test_geodesics_schwarzschild`: radial null,
   photon-sphere, and weak-bending gates must pass. Halve `max_step`; a
   nonconvergent result falsifies the ordinary-ray claim.
6. `./tests/relativity/test_geodesics_kerr`: maximum constraint must remain
   below `1e-10`, E/Lz must remain zero for the tested metric, and Carter must
   remain NaN.
7. Rebuild all relativity tests with ASan/UBSan. Any diagnostic invalidates
   this gate.

## Not completed

- Observer/tetrad construction, local photon or timelike initialization,
  frequency normalization, and future-direction classification.
- Carter constant, separated Kerr/Mino-time integration, and long-time
  structure-preserving timelike integration.
- Kerr–Schild coordinates, physical horizon crossing, interior cutoff
  execution, or singularity treatment.
- Disk/material models, radiative transfer, reference renderer, image/movie
  pipeline, Solar Local Patch, WASM/GPU, UI, or visual regression.
- DE440 external-fixture assertions and GCC/Linux verification.

## Result

Phase 1 passes on the verified source state. The only next allowed
implementation phase is Phase 2: observer and tetrad construction with local
round-trip and frequency-normalization validation.
