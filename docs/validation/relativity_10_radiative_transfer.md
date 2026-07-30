# Validation: Solar relativity Phase 5 invariant radiative transfer

## Claim

Solar provides an installable, unpolarized covariant radiative-transfer
kernel for observer-to-past ray tracing, explicit frequency/redshift
validation, analytic Kerr circular-disk and optically thin torus matter
models, and bounded thin-disk surface composition. The tested constant
coefficient solutions agree with their analytic values, and common
Boyer–Lindquist/Kerr–Schild disk, torus, velocity, redshift, and surface
normal values agree within the declared gates.

## Model boundary

This phase is fixed-background general relativity in geometrized spacetime
units with caller-selected material, frequency, and intensity units. It is
unpolarized and does not implement scattering, returning radiation, radiation
feedback, self-gravity, magnetic evolution, GRMHD, an absolute Planck
spectrum, or a detector bandpass.

`AnalyticCircularDiskFluid` is a surface-bounded circular flow with a
controlled zero-torque temperature shape, not a full Page–Thorne or GRMHD
disk. `AnalyticOpticallyThinTorus` is a compact kinematic Gaussian test model,
not a Fishbone–Moncrief equilibrium. Thin-disk event localization remains the
Phase 6 renderer's responsibility. No renderer, image, movie, RGB transform,
CLI transfer command, WASM, or GPU path is claimed here.

## Equations and conventions

Solar keeps signature `(-,+,+,+)`, canonical photon momentum `p_mu`, and
future-directed photon momentum while the ray position is integrated toward
the past. The observer supplies positive segment length
`ds=-d lambda`. Local frequency and redshift are

```text
nu_emit_normalized = -p_mu u_emit^mu
nu_emit = nu_observer nu_emit_normalized
g = nu_observer / nu_emit = 1 / nu_emit_normalized
```

The invariant transfer quantities are

```text
I = I_nu / nu^3
J = j_nu / nu^2
A = nu alpha_nu
dI/dlambda = J - A I
```

For one constant observer-to-past segment:

```text
tau = A ds
e = exp(-tau)
phi(tau) = (1-exp(-tau))/tau
delta_I = T J ds phi(tau)
I <- I + delta_I
T <- T e
```

The implementation uses `expm1` and a series branch near zero. The negative
exponent is required by foreground-to-background formal-solution
composition; it is not obtained by directly changing the sign of the
future-affine differential equation.

For grey volume emission:

```text
J = emissivity_scale density / nu_emit^2
A = absorption_scale density nu_emit
```

At a thin-disk surface:

```text
T_obs = g T_emit
I_nu,obs = g^3 I_nu,emit
I_bol,obs = g^4 I_bol,emit
```

Semi-transparent sheets compose in linear intensity with
`1-exp(-tau_surface)` and shared foreground transmission. No sRGB alpha
blend enters the physical kernel.

## Analytic fixtures

The flat constant-coefficient fixtures are:

| Fixture | Analytic result |
|---|---:|
| `J=2.5, A=0, ds=4` | `I=10`, `T=1` |
| `J=3, A=0.4, ds=5` | `I=6.4849853757254046` |
| same source transmission | `T=0.1353352832366127` |
| `J=4, A=2, ds=500` | `I=2`, `T=0` by physical saturation |
| foreground `I=1.25, T=0.4`, then `J=2, A=0.5, ds=3` | `I=2.4929917437625124`, `T=0.08925206405937193` |

Subdividing the `J=3, A=0.4, ds=5` interval into 1, 2, 7, and
100 segments produces the same state within the measured gate.

The Doppler fixture uses Minkowski momentum `p_mu=(-1,1,0,0)` and a
`v=0.3` emitter. It checks `nu_emit=gamma(1-v) nu_observer`, its reciprocal
redshift, and grey coefficients at `nu_observer=230 GHz`.

The thin-disk fixture has `T_emit=8`, emitted specific intensity `6`,
emitted bolometric intensity `10`, and normalized `g=0.5`. It therefore
requires `T_obs=4`, specific intensity `0.75`, and bolometric intensity
`0.625`. Two semi-transparent sheets with `tau=log(2)` require specific
intensity `0.3984375`, bolometric intensity `0.322265625`, and transmission
`0.25`.

## Numerical gates

All actual values below are from the clean local Release run.

| Gate | Actual | Required |
|---|---:|---:|
| Constant/formal-solution maximum normalized error | `1.7813505037117663e-16` | `<5e-14` |
| 1/2/7/100 subdivision maximum error | `1.3695920163902828e-15` | `<5e-14` |
| Boosted-emitter redshift maximum error | `1.6293619469317347e-16` | `<2e-14` |
| BL/KS disk/torus velocity component error | `1.7826011424196474e-16` | `<1e-10` |
| Fluid four-velocity norm error | `2.2204460492503131e-16` | `<1e-10` |
| BL/KS surface-normal component/norm error | `1.5407439555097887e-33` | `<1e-10` |
| Surface redshift/composition maximum error | `5.5511151231257827e-16` | `<1e-10` |
| Installed-consumer constant-transfer error | `2.2204460492503131e-16` | `<5e-14` |

Disk tests cover the ISCO default, profile values, inner/outer/surface
vacuum support, both orbit senses, spin signs, and parameter mismatch.
Torus tests cover center values, Gaussian profile, north/south symmetry,
compact cutoff, and BL/KS agreement. Every accepted fluid velocity is checked
future-directed and unit timelike.

Surface tests cover all recorded fields, front/back classification, opaque
closure, two ordered semi-transparent sheets, exactly eight accepted
crossings, explicit ninth-crossing refusal, vacuum samples that do not
consume image order, failure atomicity, and BL/KS normal/redshift agreement.

## Commands and results

Local platform:

```text
Darwin 23.6.0 arm64
Apple Clang 16.0.0
```

Release commands:

```bash
make clean
make -j4 test
make test-external-consumer
git diff --check
```

Sanitizer build:

```bash
make clean
make -j4 \
  CXXFLAGS='-std=c++17 -O1 -g -Wall -Wextra -Iinclude \
  -fsanitize=address,undefined -fno-omit-frame-pointer' \
  tests/relativity/test_radiative_transfer \
  tests/relativity/test_emission_models \
  tests/relativity/test_transfer_failures \
  tests/relativity/test_fluid_models \
  tests/relativity/test_thin_disk
```

Each focused executable ran with:

```text
ASAN_OPTIONS=detect_leaks=0:halt_on_error=1
UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1
```

Results:

```text
clean Release relativity: 2462 passed, 0 failed across 39 executables
clean Release legacy:       67 passed, 0 failed across 5 reporting executables
clean Release total:      2529 passed, 0 failed
Phase 5 focused Release:   195 passed, 0 failed across 5 executables
Phase 5 ASan/UBSan:         195 passed, 0 failed across 5 executables
optional DE440 fixture:
  skipped because data/de440.asc is absent
```

The seven legacy executables include the two non-assertion-reporting
`test_de` and `test_horizons` programs. The clean build retains one pre-existing
`src/cr3bp.cpp:388` unused-variable warning outside the Phase 5 call path; no
Phase 5 source emits a compiler warning.

The clean installed consumer reported:

```json
{"solar_version":"0.2.0-alpha.1","physics_contract":"relativity-v3-phase2","samples":128,"left_edge":-4.0962666587138692,"right_edge":6.1381557247154532,"separated_steps":4,"separated_constraint":1.2394276503472693e-16,"ks_steps":7,"ks_constraint":4.0906515605215856e-17,"ks_inverse_error":1.1102230246251565e-16,"transfer_intensity":3.1074793594062804,"transfer_transmission":0.22313016014842982,"disk_temperature":8,"torus_density":3,"surface_specific_intensity":6.0000000000000036,"surface_crossings":1}
```

Verified implementation commit:
`9f9c144f86bb296e1887522c0b1bb9317f4f41ad`.

The PR-head Linux/GCC release candidate passed the full Release build, test
suite, installed consumer, CLI, and sample-command gates on Ubuntu 24.04
x86_64 with GCC 13.3.0 at commit
`a65befbc1ae1462b52297b2a2e5c9e82e040eb3c`:
<https://github.com/TT1nKer/solar/actions/runs/30536006406>.

Linux retained the constant/subdivision maxima above and measured
`max_surface_normal_error=7.7037197775489434e-34`. Its installed consumer
reported transfer intensity `3.1074793594062804`, transmission
`0.22313016014842982`, disk temperature `8`, torus density `3`, surface
specific intensity `6.0000000000000036`, and one surface crossing. No Linux
sanitizer, RTX 3080, or DE440-data validation is claimed.

## Mutation sensitivity

Two temporary, uncommitted mutations were applied and restored:

1. Replacing `exp(-tau)` and `expm1(-tau)` with the wrong positive exponent
   caused `test_radiative_transfer` to fail 9 assertions and
   `test_thin_disk` to fail 3 assertions.
2. Replacing the full `u_KS=J u_BL` transform with a component copy caused
   `test_fluid_models` to terminate on its explicit non-timelike velocity
   guard and `test_thin_disk` to fail its BL/KS event check.

Both restored implementations passed their focused Release and sanitizer
sets. Neither mutation is present in the repository.

## References

- Bronzwaer et al., RAPTOR I, arbitrary-spacetime general-relativistic
  radiative transfer: <https://arxiv.org/abs/1801.10452>
- Mościbrodzka and Gammie, `ipole`, analytic constant-coefficient covariant
  transport: <https://arxiv.org/abs/1712.03057>
- Event Horizon Telescope Collaboration, radiative-transfer code
  verification: <https://doi.org/10.3847/1538-4357/ab96c6>

These are equation and verification references only. No external
implementation code was copied.

## Result

`PASSED`. All local Phase 5 physics, failure, installation, mutation,
Release, and sanitizer gates pass, and the PR-head Linux/GCC Release,
installed-consumer, CLI, and sample-command gates pass.

## Limitations

- Coefficients are constant per caller-supplied segment; adaptive sampling and
  accepted-step coupling belong to Phase 6.
- The disk profile is a controlled analytic approximation and the torus is a
  kinematic Gaussian. Neither represents a calibrated astrophysical plasma.
- Material and intensity units are caller-selected. There is no absolute
  physical spectral calibration or camera bandpass.
- Thin-disk surface event localization and repeated ray/disk intersection
  detection are not part of this phase.
- Polarization, scattering, returning radiation, radiation feedback,
  self-gravity, magnetic fields, and GRMHD are not implemented.
- Near-extremal, long-path, and dense multi-crossing transfer need broader
  independent characterization before scientific production use.

## Fastest falsification

Run:

```bash
./tests/relativity/test_radiative_transfer
./tests/relativity/test_emission_models
./tests/relativity/test_transfer_failures
./tests/relativity/test_fluid_models
./tests/relativity/test_thin_disk
make test-external-consumer
```

Any constant-solution or subdivision error above its gate, non-physical
frequency accepted as valid, invalid sample converted to vacuum, non-timelike
fluid velocity, BL/KS vector mismatch, wrong surface power of `g`, incorrect
foreground transmission, unbounded crossing growth, or missing installed
symbol invalidates Phase 5.
