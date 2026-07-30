# Validation: Solar relativity Phase 3 separated Kerr solver

## Claim

Solar provides a public, installable Kerr separated-geodesic solver in Mino
time. For the tested exterior Boyer–Lindquist null and timelike families, its
worldlines agree with the independent generic Hamiltonian solver at common
events while preserving the Phase 3 constraint and Carter gates.

## Equations and conventions

The solver uses coordinate order `(t,r,theta,phi)`, signature `(-,+,+,+)`,
canonical momentum `p_mu`, `E=-p_t`, `Lz=p_phi`, and
`mu=cos(theta)`. Mino time satisfies

```text
d lambda / d gamma = Sigma
P(r) = E (r^2 + a^2) - a Lz
R(r) = P(r)^2 - Delta [mass_sq r^2 + (Lz-aE)^2 + Q]
U(mu) = Q(1-mu^2)
        - mu^2 [a^2 (mass_sq-E^2) (1-mu^2) + Lz^2]
dr/dgamma = s_r sqrt(R)
dmu/dgamma = s_mu sqrt(U)
```

The azimuth, coordinate-time, and affine equations are the literal v3
separated equations. The ordinary path is a fixed-size DOPRI5 state with no
per-step allocation. A trial entering `R<0` or `U<0` is rejected, the simple
root is bracketed and refined, and the direction is changed only after a
controlled root release. A synchronized smooth phase state carries
`(r,dr/dgamma,mu,dmu/dgamma)` through the non-Lipschitz square-root
neighborhood so a double-precision coordinate does not lose the turning
phase. Subsequent phase crossings are localized from velocity sign changes.
Small derivative roots terminate as `NearCriticalOrbit`; the implementation
never continues with `sqrt(max(V,0))`.

References:

- Fujita and Hikida, analytic Kerr geodesics in Mino time:
  <https://arxiv.org/abs/0906.1420>
- Dexter and Agol, independent Kerr geodesic reference implementation:
  <https://arxiv.org/abs/0903.0620>

No external implementation code was copied.

## Common-event fixtures

All values below are from the clean Release run. Coordinate error combines
`dt/M`, `dr/M`, `dtheta`, and wrapped `dphi`; affine error is normalized by
`M`.

| Fixture | Worldline error | Fine-step error | Affine error | Separated constraint | Separated Carter | Turns `(r,mu)` |
|---|---:|---:|---:|---:|---:|---:|
| null Kerr, `chi=+0.5` | `2.23464e-15` | `5.32936e-15` | `2.22045e-15` | `1.41104e-16` | `6.02254e-16` | `(0,0)` |
| null Kerr, `chi=-0.5` | `2.22047e-15` | `6.31878e-13` | `3.99680e-15` | `1.24952e-16` | `4.51691e-16` | `(0,0)` |
| timelike Kerr, `chi=+0.5` | `2.08750e-14` | `8.02783e-13` | `1.82077e-14` | `1.20526e-16` | `9.12976e-16` | `(0,0)` |
| null Schwarzschild | `1.48972e-15` | `1.24429e-14` | `1.33227e-15` | `1.27375e-16` | `5.53521e-16` | `(0,0)` |
| radial return, Schwarzschild | `1.18898e-12` | `1.22805e-12` | `9.73444e-13` | `1.23681e-16` | `0` | `(1,0)` |
| polar return, Schwarzschild | `2.27726e-11` | `2.06562e-11` | `1.52234e-11` | `1.54540e-16` | `8.26469e-16` | `(0,1)` |

Release worldline P95 and maximum were both `2.27726e-11`. The largest
Hamiltonian reference constraint was `8.50112e-16`; the largest Hamiltonian
Carter error was `1.15706e-14`. Every ordinary error is below `1e-8`, the
maximum is below `1e-7`, and all invariant gates are below `1e-10`.

The radial scattering probe located
`r_min=4.4533631938113549M`, recorded 37 rejected and 971 accepted steps, and
returned through `r=10M`. Halving the maximum Mino step changed coordinate
time by `6.25278e-12M` and azimuth by `5.51115e-13`. The polar return recorded
26 rejected and 469 accepted steps. The exact Schwarzschild photon-sphere
fixture terminates as `NearCriticalOrbit` without incrementing a simple-turn
counter. The phase-crossing unit test also verifies that a root at a step
start is not counted twice.

## Commands and results

```bash
make clean
make -j4 test
make test-external-consumer

cmake -S . -B build-phase3-sanitize \
  -DSOLAR_BUILD_CLI=OFF \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer"
cmake --build build-phase3-sanitize --parallel
```

The five Phase 3 executables were compiled separately against the sanitizer
libraries with the same sanitizer flags and run:

```text
clean Release relativity:  2076 passed, 0 failed across 27 executables
clean Release legacy:        67 passed, 0 failed across 5 executables
Phase 3 focused Release:    163 passed, 0 failed across 5 executables
Phase 3 ASan/UBSan:          163 passed, 0 failed across 5 executables
sanitizer worldline P95/max: 2.32107e-11 / 2.32107e-11
external installed consumer: 4 accepted steps, constraint 1.23943e-16
optional DE440 fixture:      skipped because data/de440.asc is absent
```

The controlled `+1e-4M` radial-coordinate mutation caused all six
common-event coordinate gates and both aggregate gates to fail. The mutation
was reverted before the passing runs.

Verified implementation commit:
`0c67826f222ca9078c49fc64b8c91b7baf9eac5a`.

Platform: Darwin 23.6.0 arm64, Apple Clang 16.0.0. No Linux, GCC, RTX 3080, or
DE440-data validation is claimed.

## Result

PASSED. Phase 3 separated/Mino evolution, simple radial and polar turns,
near-critical refusal, installed public consumption, and Hamiltonian
common-event agreement satisfy the local evidence gate.

## Limitations

- Exterior Boyer–Lindquist chart only; invalid BL points are not physical
  horizon crossings.
- No Kerr–Schild chart, chart transform, interior evolution, or singularity
  treatment.
- No analytic elliptic-function backend or fundamental-frequency API.
- Near-axis trajectories with nonzero `Lz` remain unsupported because the BL
  azimuth term is coordinate-singular.
- The smooth turning phase uses a `1e-10` normalized potential-drift ceiling,
  matching the CPU geodesic constraint gate; near-critical failures remain
  explicit rather than projected away.
- Validation covers six common-event families, moderate spins, and one simple
  turn of each type. Long bound, near-extremal, and repeated-turn scientific
  sweeps remain future work.

## Fastest falsification

Run:

```bash
./tests/relativity/test_kerr_separated
./tests/relativity/test_kerr_separated_crosscheck
./tests/relativity/test_kerr_separated_turning_phase
make test-external-consumer
```

Any missed forbidden-potential rejection, repeated root count, silent
near-critical flip, installed-symbol failure, invariant above `1e-10`, or
common-event error above its gate invalidates Phase 3.
