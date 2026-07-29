# Solar API Audit for Relativity

Audit base: `main` at `1268d89f5f84779d2fc66a1487280af586d1a245`.
Audit platform: macOS 14.8.7 arm64, Apple Clang 16.0.0.

## Repository tree (three levels)

Generated from the real checkout, excluding `.git` and generated
objects/binaries:

```text
.
├── .github/workflows/ci.yml
├── CHANGELOG.md
├── FEATURE_STATUS.md
├── LIMITATIONS.md
├── Makefile
├── README.md
├── STRATEGY.md
├── cli/main.cpp
├── docs
│   ├── readme-cover.svg
│   └── validation
│       ├── 01_kepler_solver.md
│       ├── 02_hohmann_earth_mars.md
│       ├── 03_energy_conservation.md
│       ├── 04_de440_earth_j2000.md
│       ├── 05_lagrange_points.md
│       ├── 06_halo_orbit_jacobi.md
│       ├── 07_frame_transforms.md
│       └── 08_lambert_solver.md
├── include/solar
│   ├── atmosphere.h
│   ├── body.h
│   ├── constants.h
│   ├── cr3bp.h
│   ├── drag.h
│   ├── ephemeris.h
│   ├── force_model.h
│   ├── frame.h
│   ├── gr_correction.h
│   ├── gravity.h
│   ├── iau_rotation.h
│   ├── integrator.h
│   ├── j2.h
│   ├── jpl_de.h
│   ├── kepler.h
│   ├── mission.h
│   ├── montecarlo.h
│   ├── moons.h
│   ├── nbody.h
│   ├── network.h
│   ├── spacecraft.h
│   ├── spherical_harmonics.h
│   ├── srp.h
│   ├── time_scale.h
│   ├── trajectory.h
│   ├── transfer.h
│   └── vehicle.h
├── src
│   ├── atmosphere.cpp
│   ├── cr3bp.cpp
│   ├── drag.cpp
│   ├── ephemeris.cpp
│   ├── force_model.cpp
│   ├── frame.cpp
│   ├── gr_correction.cpp
│   ├── gravity.cpp
│   ├── iau_rotation.cpp
│   ├── integrator.cpp
│   ├── j2.cpp
│   ├── jpl_de.cpp
│   ├── kepler.cpp
│   ├── mission.cpp
│   ├── montecarlo.cpp
│   ├── moons.cpp
│   ├── nbody.cpp
│   ├── network.cpp
│   ├── spacecraft.cpp
│   ├── spherical_harmonics.cpp
│   ├── srp.cpp
│   ├── time_scale.cpp
│   ├── trajectory.cpp
│   └── transfer.cpp
└── tests
    ├── test_de.cpp
    ├── test_horizons.cpp
    ├── test_kepler.cpp
    ├── test_montecarlo.cpp
    ├── test_network.cpp
    └── test_validation.cpp
```

## Current math and state types

`include/solar/body.h` owns a field-based Euclidean `solar::Vec3` with
arithmetic, norm, dot, cross, and unchecked normalization. `solar::State`
contains Cartesian `pos` in km and `vel` in km/s. `Body` combines physical,
hierarchy, orbital-element, and mutable state data.

`include/solar/frame.h` owns row-major `Mat3` rotations and transforms a
`State`. There is no four-vector, four-covector, 4×4 matrix, tensor variance,
geometrized-unit wrapper, or automatic differentiation type.

The legacy `Vec3` is reusable only in the Newtonian boundary. Its implicit
km/km·s⁻¹ semantics and Euclidean dot make it unsafe as a spacetime API.

## Integrator APIs

`include/solar/integrator.h` exposes:

- `rk4_step(states, masses, accel, dt)` for second-order N-body states;
- `verlet_step(states, masses, accel, prev_accels, dt)` for in-place
  velocity Verlet;
- `dopri5_step(...) -> AdaptiveResult` for N-body `State`;
- `dopri5_generic_step(vector<double>, t, dt, f, atol, rtol)
  -> GenericAdaptiveResult` for flat first-order ODEs.

The generic DOPRI5 stages are potentially reusable as a numerical reference,
but its public shape is not sufficient for the planned solver: it uses one
absolute and relative tolerance across all components, allocates vectors at
every stage, provides no dense output/event localization, exposes no failure
reason, and does not validate derivative size or finiteness. It must be adapted
behind a typed relativity integrator rather than used directly as the public GR
API.

`NBodySim::step` advances `time` even when its DOPRI5 step is rejected. This is
a correctness risk in the legacy path but is outside this additive Phase 0A
change. A later geodesic solver must not inherit that state transition.

## Force model API

`ForceModel::compute(bodies, time, acc)` adds three-dimensional acceleration to
a pre-sized vector; `potential_energy` defaults to zero. `NBodySim` sums these
models. `GRCorrection` is explicitly a 1PN Schwarzschild acceleration
correction in Newtonian coordinates, not a metric or geodesic model.

No `ForceModel` implementation can represent covariant Hamilton equations,
null motion, horizon events, or coordinate charts. It is an adapter-side
legacy capability only and must not become the fixed-background GR abstraction.

## Frame and time APIs

`Frame` distinguishes Ecliptic J2000, ICRF, and body-fixed 3D frames.
`transform_state` rotates Cartesian position/velocity with optional IAU body
rotation. It has no four-vector/covector Jacobian contract.

`Epoch` stores Julian date plus `TimeScale` (`TDB`, `TT`, `UTC`, `TAI`) and
provides conversions. These types can be reused at the eventual Solar adapter
boundary, but not as affine parameter or spacetime coordinate time without an
explicit conversion policy.

## Trajectory and mission APIs

`trajectory.h` contains porkchop, Lambert-based launch window, patched flyby,
and multi-flyby request/result structs. `mission.h` contains event-list mission
plans and delta-v/fuel reports. All use Solar's km, km/s, days, and Julian-date
conventions.

These are L2/L3 Solar workflows. They may supply boundary conditions to a
future local strong-field adapter, but their state types cannot be reused
inside the GR kernel.

## CLI dispatch

`cli/main.cpp` is a 1,244-line manual dispatcher. Each command is a static
`cmd_*` function; `main` rewrites argv after optional `--de`, then selects
commands with an `if` chain inside one `try/catch`. Commands print directly to
stdout/stderr and return `0` or `1`; top-level `std::exception` becomes
`Error: <message>`.

No CLI registry or structured result boundary exists. A later relativity CLI
should add one narrow command entry without moving unrelated commands in the
same change. The file already exceeds the soft review size and merits a future
dispatch split, but Phase 0A does not touch it.

## Tests and assertions

Tests are independent C++ executables with local counters/check functions and
integer exit status. There is no external framework. `test_horizons` is a
print-only comparison and always returns success. `test_de` requires an
external DE440 fixture for its eight assertions.

The original Make recipe ran every executable in a shell loop but returned only
the last command's status, allowing an early failing test to be hidden. Phase
−1 changes the aggregate recipe to retain failure.

## Error handling

The codebase mixes:

- exceptions for parsing/CLI boundaries;
- `bool valid` or `converged` result flags;
- `nullptr` for lookup failure;
- silent early return in some force models;
- unchecked division/normalization in low-level vector math;
- stderr plus exit codes in tests and CLI.

Phase 0A uses exceptions only for programmer/input-domain violations in small
value operations. Later numerical solvers need explicit termination and
diagnostic enums rather than exceptions for ordinary physical outcomes.

## Build discovery

The original Makefile used `$(wildcard src/*.cpp)` and
`$(wildcard tests/*.cpp)`, so nested relativity code was invisible. Phase −1
uses sorted recursive `find` expressions. Existing targets remain
`libsolar.a`, `solar`, `make test`, and `make clean`; no dependency or build
system is added.

## Reuse, adapt, and add

Directly reusable:

- physical constants after explicit unit conversion;
- `Epoch`/`TimeScale` at an external adapter boundary;
- build/test conventions and static-library packaging;
- documented validation-report style.

Adapt only:

- generic DOPRI5 coefficients/algorithm, after typed state, component scales,
  finite validation, dense output, and event roots exist;
- `Frame`/`Mat3` at the Solar-to-local-chart boundary;
- `Body`/ephemeris and trajectory/mission results as external initial
  conditions;
- CLI dispatch for a later thin entry point.

Must be new:

- geometrized-unit boundary;
- four-dimensional vector/matrix operations;
- covariant and contravariant wrappers;
- phase-space and diagnostic sample types;
- forward-mode four-variable automatic differentiation;
- later: metric/chart, Hamiltonian, event, observer/tetrad, and radiative
  transfer contracts.

## Compatibility risks

1. `solar::Vec3` and future `solar::relativity::Vec3` share a short name but
   different storage/units. Namespace qualification and no implicit conversion
   are required.
2. Recursive test discovery changes CI truthfulness: previously invisible or
   masked failing tests will now fail `make test`.
3. Recursive source discovery can accidentally archive an experimental `.cpp`
   placed below `src`; source placement now implies build participation.
4. `PhaseSpaceState::x` is typed as a contravariant vector by the v3 contract,
   although coordinates are mathematically points. Changing this later would
   be an API migration.
5. Legacy compilation has seven warnings (one CR3BP unused variable and six
   unused DOPRI constants); Phase 0A records but does not mix their cleanup into
   the relativity work.
6. The docs claim five test files/64 assertions while the repository contains
   six executables including a print-only Horizons program. Documentation
   counts need separate cleanup.

## Minimum directory change

```text
include/solar/relativity/
  dual4.h
  math.h
  types.h
  units.h
src/relativity/
  units.cpp
tests/relativity/
  test_dual4.cpp
  test_math.cpp
  test_types.cpp
  test_units.cpp
docs/relativity/
  CONVENTIONS.md
  SOLAR_API_AUDIT.md
docs/validation/
  relativity_00_baseline.md
```

No empty future module or application directory is created.
