# Solar Relativity Phase −1/0A Design

## Scope

This is the first bounded implementation slice of
`SOLAR_RELATIVITY_KERR_完整实施主提示词_v3.md`. It establishes the audited
Solar baseline and the L0 conventions required by later fixed-background
relativity work. It does not implement a metric, geodesic equation, black-hole
classification, renderer, adapter, or UI.

The authoritative implementation remains C++17 `double`. No runtime dependency
is added.

## Baseline and build closure

The untouched `main` commit is `1268d89f5f84779d2fc66a1487280af586d1a245`.
Its library and CLI build on macOS 14.8.7 with Apple Clang 16, but
`tests/test_validation.cpp` omits its direct `<sstream>` dependency. After that
is fixed, the suite also reveals that the Makefile ignores a failed child test
when a later test succeeds. The optional DE440 test must skip only when its
default fixture is absent; an explicitly supplied invalid path remains a
failure.

The Makefile will discover `src/**/*.cpp` and `tests/**/test_*.cpp`
recursively, preserve `libsolar.a` and `solar`, and return nonzero if any child
test fails.

## L0 module boundaries

- `units.h/.cpp` owns SI-to-geometrized conversions and input-domain
  validation. It reuses Solar's existing `G` and `C_LIGHT` constants after
  explicit kilometre-to-metre conversion.
- `math.h` owns fixed-size numeric vectors, 4×4 matrices, finite/max norms,
  explicit metric contraction, multiplication, and checked inversion.
- `types.h` owns tensor-variance wrappers and geodesic state/sample value
  types. `Contravariant4` and `Covariant4` are intentionally non-interchangeable.
- `dual4.h` owns a scalar value and its four coordinate derivatives, including
  the arithmetic and elementary functions needed to differentiate a metric in
  Phase 0B.

All four modules depend only on the C++ standard library, except `units.cpp`,
which reads existing Solar constants. No legacy `solar::Vec3`, `State`,
`NBodySim`, or force-model type is changed.

## Mathematical contract

- Four-dimensional spacetime with signature `(-,+,+,+)`.
- Coordinate component zero is always time.
- Supported chart order declarations are Minkowski Cartesian `(t,x,y,z)`,
  Boyer–Lindquist `(t,r,theta,phi)`, and Cartesian Kerr–Schild `(t,x,y,z)`.
- The physical spin axis is world `+z`; display-axis conversion is outside the
  physics layer.
- Integration units use `G=c=1`.
- `M_length = G M_kg / c^2` and `M_time = G M_kg / c^3`.
- A spacetime dot product is unavailable without an explicit metric. The only
  convenience contraction in Phase 0A is explicitly named
  `minkowski_dot_minus_plus_plus_plus`.
- `PhaseSpaceState::x` follows the v3 public contract and uses
  `Contravariant4`; a separate affine-point type is deferred to an intentional
  API revision rather than silently diverging from the contract.
- Null samples represent unavailable proper time with quiet NaN. Timelike
  unit-mass states use Hamiltonian target `-1/2`; null states use `0`.

## Failure behavior

- Unit construction rejects non-finite or non-positive mass with
  `std::invalid_argument`.
- Unit conversions reject non-finite inputs with `std::invalid_argument`.
- Vector division rejects zero or non-finite divisors.
- Matrix inversion rejects non-finite and singular/ill-conditioned inputs with
  `std::domain_error`; it uses partial pivoting.
- `Dual4::variable` rejects derivative indices outside `[0,3]`.
- `sqrt` rejects negative values, `log` rejects non-positive values, and
  division rejects a zero denominator. Singular derivatives are not hidden.

These checks protect mathematical domains; they do not clamp or alter a
physical solution.

## Verification

Tests are standalone C++ executables using the repository's existing
pass/fail-counter style:

- `test_units`: solar-mass reference values, round trips, velocity scaling, and
  invalid input.
- `test_math`: vector arithmetic, finite/max norm, explicit Minkowski
  contraction, matrix multiplication/inversion, and singular rejection.
- `test_types`: compile-time variance separation and null/timelike sample
  conventions.
- `test_dual4`: product, quotient, chain rule, elementary functions,
  multivariable derivatives, and domain failures.

Each production behavior is introduced only after its focused test is observed
failing for the missing feature. Final verification is `make clean`, `make`,
`make test`, direct Phase 0A test executables, README smoke commands, and
`git diff --check`.

## Compatibility and deferred risks

This slice is additive. Existing public headers and binaries retain their names
and behavior. The only test behavior change is honest optional-fixture handling
and reliable failure propagation.

The legacy generic DOPRI5 remains unsuitable for later Hamiltonian evolution:
it has scalar tolerances, no dense output/event roots, and no typed phase-space
contract. It may be adapted later but is not reused as the Phase 0A foundation.
Legacy `solar::Vec3` remains the Newtonian kilometre-based type; automatic
conversion to relativity `Vec3` is deliberately omitted to prevent unit mixing.
