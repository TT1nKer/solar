# Validation: Solar relativity Phase 4 Cartesian Kerr–Schild

## Claim

Solar provides an installable, horizon-regular Cartesian ingoing
Kerr–Schild metric, a canonical Boyer–Lindquist/Kerr–Schild overlap
transform, and explicit outer-horizon and interior-cutoff events. For the
tested subextremal Kerr null and timelike families, the Kerr–Schild
Hamiltonian flow agrees with the independent Boyer–Lindquist flow at common
exterior events and continues a timelike plunge from the outer horizon to the
declared interior cutoff.

## Equations and conventions

The metric uses coordinate order `(t,x,y,z)`, signature `(-,+,+,+)`,
geometrized units `G=c=1`, and canonical covariant momentum `p_mu`. The
positive implicit radius is defined by

```text
r^4 - (x^2+y^2+z^2-a^2) r^2 - a^2 z^2 = 0
H = M r^3 / (r^4+a^2 z^2)
l_mu = [1, (r x+a y)/(r^2+a^2),
           (r y-a x)/(r^2+a^2), z/r]
g_mu_nu = eta_mu_nu + 2 H l_mu l_nu
g^mu_nu = eta^mu_nu - 2 H l^mu l^nu
```

The implementation evaluates the positive root with a cancellation-safe
branch, constructs the analytic inverse from the Kerr–Schild null form, and
uses `Dual4` automatic differentiation for all inverse-metric derivatives.
The independent analytic radius gradient is validated separately.

The safe exterior chart transform is

```text
t_KS = t_BL + F_t(r)
phi_tilde = phi_BL + F_phi(r)
x = (r cos(phi_tilde)-a sin(phi_tilde)) sin(theta)
y = (r sin(phi_tilde)+a cos(phi_tilde)) sin(theta)
z = r cos(theta)

dF_t/dr = 2 M r / Delta
dF_phi/dr = a / Delta
```

The master prompt contains a sign conflict: its differential identity requires
`dF_phi/dr=+a/Delta`, while a later displayed antiderivative has a leading
minus and differentiates to the opposite sign. Phase 4 follows the
differential identity and the documented ingoing transform:

```text
F_phi = a/(r_+-r_-) log((r-r_+)/(r-r_-))
```

The full forward Jacobian is differentiated from the coordinate expression.
Canonical momenta use `p_KS=(J^-1)^T p_BL`; no spherical momentum is copied or
treated as a Cartesian vector.

References:

- SpECTRE Kerr–Schild coordinate transformation:
  <https://spectre-code.org/classgr_1_1KerrSchildCoords.html>
- GRay2 covariant ray-tracing formulation:
  <https://arxiv.org/abs/1706.07062>
- Einstein Toolkit exact Kerr–Schild initial data:
  <https://einsteintoolkit.org/thornguide/EinsteinInitialData/Exact/documentation.html>

No external implementation code was copied.

## Numerical gates

All values below are from the clean Release run.

| Gate | Actual | Required |
|---|---:|---:|
| Analytic radius gradient vs independent `Dual4` | `4.44089e-16` | `<1e-12` |
| AD inverse derivative vs five-point difference | `6.92227e-11` | `<3e-8` |
| Position round trip | `8.88178e-16` | `<1e-10` |
| Forward/reverse Jacobian identity | `1.02546e-15` | `<5e-12` |
| Canonical momentum round trip | `9.43690e-16` | `<1e-10` |
| Covector/vector pairing invariance | `8.74301e-16` | `<1e-10` |
| Hamiltonian invariance across charts | `6.66134e-16` | `<1e-10` |
| `dF_t/dr` finite-difference error | `5.85154e-12` | `<1e-9` |
| `dF_phi/dr` finite-difference error | `8.17645e-13` | `<1e-9` |
| Ordinary BL/KS position P95 | `6.64417e-12` | `<1e-8` |
| Ordinary BL/KS momentum P95 | `6.24937e-12` | `<1e-8` |
| Ordinary maximum Hamiltonian constraint | `3.21179e-12` | `<1e-10` |
| Ordinary maximum stationary-energy drift | `0` | `<1e-10` |
| Ordinary maximum axial-angular-momentum drift | `1.31681e-13` | `<1e-10` |
| Near-horizon position error at `r_++0.02M` | `6.13639e-10` | `<1e-7` |
| Near-horizon momentum error | `2.47107e-10` | `<1e-7` |
| Horizon event radius error | `<5e-9 M` | `<5e-9 M` |
| Horizon Hamiltonian constraint | `9.91845e-15` | `<1e-10` |
| Horizon tetrad error | `2.53961e-16` | `<1e-10` |
| Interior cutoff radius error | `1.23903e-12 M` | `<5e-9 M` |
| Interior Hamiltonian constraint | `1.30158e-13` | `<1e-10` |
| Affine interval evolved inside horizon | `0.855515M` | `>0` |

The plunge first terminates as `HorizonCrossing`, restarts from the exact
localized event state, and later terminates as `InteriorCutoff`. The default
cutoff is `max(0.05M, configured)`; it is a declared numerical/model boundary,
not a claim of reaching a physical singularity.

## Commands and results

```bash
make clean
make -j4 test
make test-external-consumer
```

The seven Phase 4 executables were also built and run with:

```text
-fsanitize=address,undefined -fno-omit-frame-pointer
ASAN_OPTIONS=detect_leaks=0:halt_on_error=1
UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1
```

Results:

```text
clean Release relativity: 2267 passed, 0 failed across 34 executables
clean Release legacy:       67 passed, 0 failed across 5 executables
clean Release total:      2334 passed, 0 failed across 39 executables
Phase 4 focused Release:   189 passed, 0 failed across 7 executables
Phase 4 ASan/UBSan:         189 passed, 0 failed across 7 executables
external installed consumer:
  KS accepted steps=7
  KS constraint=4.0906515605215856e-17
  KS inverse error=1.1102230246251565e-16
optional DE440 fixture:
  skipped because data/de440.asc is absent
```

As a mutation check, transposing the forward covector-transform Jacobian index
caused the chart suite to fail 3 assertions and the BL/KS common-event suite
to fail 20 assertions. Restoring the canonical inverse-transpose made both
suites pass again under ASan/UBSan.

Verified implementation commit:
`74087094dc98d2278f17e6350e15cb5c87557d00`.

Local platform: Darwin 23.6.0 arm64, Apple Clang 16.0.0. Phase 4 Linux/GCC CI
is required before merge. No Linux sanitizer, RTX 3080, or DE440-data
validation is claimed.

## Result

LOCAL PASSED. The Cartesian Kerr–Schild metric, canonical overlap transform,
chart-aware invariant monitoring, horizon localization/restart, positive
interior evolution, installed public consumption, and BL/KS common-event
agreement satisfy the local Phase 4 evidence gate.

## Limitations

- Only positive-radius, subextremal Kerr with `abs(chi)<1` is supported.
  Extremal/near-extremal characterization and the negative-radius extension
  are not claimed.
- BL/KS conversion is deliberately limited to a finite exterior overlap away
  from the polar axis. Horizon crossing uses states already in Cartesian
  Kerr–Schild coordinates.
- Interior evolution stops at `max(0.05M, configured)`; there is no ring
  singularity, quantum-gravity, or maximal-extension model.
- Validation covers moderate spins and short null/timelike trajectories.
  Long bound timelike evolution still lacks a structure-preserving integrator.
- Event detection brackets endpoint sign changes; tangencies or multiple
  crossings within one accepted step remain a generic integrator limitation.
- The optional default DE440 data fixture and Linux sanitizers remain
  unverified on this machine.

## Fastest falsification

Run:

```bash
./tests/relativity/test_kerr_schild
./tests/relativity/test_kerr_schild_derivatives
./tests/relativity/test_kerr_chart_transform
./tests/relativity/test_kerr_bl_ks_crosscheck
./tests/relativity/test_geodesics_kerr_schild
make test-external-consumer
```

Any non-null Kerr–Schild one-form, failed inverse identity, derivative error
above its gate, canonical-pairing or Hamiltonian mismatch, common-event error
above its gate, false horizon/interior termination, non-positive inside
affine interval, or installed-symbol failure invalidates Phase 4.
