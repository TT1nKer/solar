# Solar Relativity Phase 2 Observer and Kerr Benchmark Design

**Status:** Approved for implementation on 2026-07-30  
**Source contract:** `SOLAR_RELATIVITY_KERR_完整实施主提示词_v3.md`  
**Previous gate:** Phase 1 `PASSED` at
`1bf7442c141e5440642af18677dc8680f96bfc9d`

## Goal

Complete the v3 Phase 2 foundation needed to initialize physically meaningful
Kerr rays and validate them against independent invariants and analytic
benchmarks:

- orthonormal observer tetrads in the `(-,+,+,+)` convention;
- Kerr Boyer-Lindquist ZAMO, static, arbitrary/look-at, and equatorial
  circular observers;
- local photon and unit-mass timelike initialization into canonical
  `(x^mu,p_mu)` state;
- observer-frequency normalization and explicit future-direction checks;
- Kerr `E`, `Lz`, and Carter `Q`, including geodesic drift diagnostics;
- ISCO, equatorial photon-orbit, marginally-bound, and circular-orbit
  quantities;
- the asymptotic Bardeen critical curve and a CPU backward-ray cross-check.

Phase 2 remains a library and validation phase. It does not add a renderer,
matter, radiative transfer, a separated/Mino-time solver, Kerr-Schild
coordinates, or physical horizon crossing.

## Governing conventions

- Signature is `(-,+,+,+)`.
- Boyer-Lindquist order is `(t,r,theta,phi)`.
- Coordinates and four-velocities are contravariant; canonical momenta are
  covariant.
- Units are geometrized with `G=c=1`; `mass_M` is the length scale.
- The observer four-velocity is tetrad leg `e_(0)^mu`.
- A physical photon momentum is always future-directed. Backward tracing uses
  a negative affine step while retaining future-directed momentum.
- Observer frequency is exactly `nu=-p_mu u_obs^mu`.
- Locally initialized photon momentum is scaled so `nu_obs=1`.
- A Boyer-Lindquist inner-boundary event used by the shadow test is a capture
  proxy only. It is never reported as physical horizon crossing.
- No default Hamiltonian-constraint projection is introduced.

## Approaches considered

### One observer class owning all Phase 2 behavior

This has a small initial API but couples tetrad algebra, observer construction,
local initialization, Kerr orbit formulas, invariant monitoring, and the
shadow benchmark. It would force Phase 3 and Phase 4 changes through one large
module and make focused tests difficult.

### Polymorphic observer and attitude framework

Virtual observer, orbit, and attitude interfaces would accommodate future
Fermi-Walker and parallel-transported cameras, but Phase 2 has only one
implemented attitude law and a few deterministic constructors. Adding a
factory hierarchy now would create unused indirection.

### Focused value types and free functions

This is the selected approach. Small L0 value/error types describe frames and
results. Focused L1 modules own tensor operations, observer construction,
initialization, Kerr constants, orbit formulas, and the analytic shadow.
The existing L2 geodesic flow receives an optional Carter evaluator and never
depends on the concrete Kerr metric.

## Layer and dependency design

### L1: spacetime algebra

`include/solar/relativity/spacetime_algebra.h` provides:

```cpp
double metric_inner_product(
    const Mat4& covariant,
    const Contravariant4& left,
    const Contravariant4& right) noexcept;

double covector_vector_pairing(
    const Covariant4& covector,
    const Contravariant4& vector) noexcept;

Covariant4 lower_index(
    const Mat4& covariant,
    const Contravariant4& vector);

Contravariant4 raise_index(
    const Mat4& contravariant,
    const Covariant4& covector);
```

The module has no knowledge of observers or Kerr. It prevents independent
modules from reproducing tensor contractions with inconsistent variance.

### L0/L1: tetrads and observers

`include/solar/relativity/observer.h` defines:

```cpp
struct Tetrad {
    // basis[a] stores e_(a)^mu
    std::array<Contravariant4, 4> basis;
};

struct ObserverFrame {
    Contravariant4 x;
    Tetrad tetrad;
};

enum class ObserverError {
    None,
    NonFiniteInput,
    InvalidMetricPoint,
    FourVelocityNotUnitTimelike,
    DegenerateSpatialSeed,
    StaticWorldlineNotTimelike,
    CircularWorldlineNotTimelike,
    TetradValidationFailure,
};

struct ObserverResult {
    ObserverError error = ObserverError::None;
    std::optional<ObserverFrame> frame;
    std::string message;

    explicit operator bool() const noexcept;
};

struct LookAtAttitude {
    Contravariant4 look_direction;
    Contravariant4 up_reference;
};
```

The public operations are:

```cpp
double tetrad_orthonormality_error(
    const Metric& metric,
    const ObserverFrame& observer);

Contravariant4 tetrad_to_coordinate(
    const Tetrad& tetrad,
    const Vec4& local_components);

Vec4 coordinate_to_tetrad(
    const Mat4& covariant,
    const Tetrad& tetrad,
    const Contravariant4& coordinate_vector);

ObserverResult make_zamo_observer(
    const KerrBoyerLindquistMetric& metric,
    const Contravariant4& x);

ObserverResult make_static_observer(
    const Metric& metric,
    const Contravariant4& x);

ObserverResult make_arbitrary_observer(
    const Metric& metric,
    const Contravariant4& x,
    const Contravariant4& four_velocity,
    const std::array<Contravariant4, 3>& spatial_seeds);

ObserverResult make_look_at_observer(
    const Metric& metric,
    const Contravariant4& x,
    const Contravariant4& four_velocity,
    const LookAtAttitude& attitude);
```

ZAMO uses the explicit v3 `alpha`, `omega`, and spatial legs. Static and
arbitrary frames use Lorentzian Gram-Schmidt:

```text
v_perp = v + g(v,e0)e0 - sum_i g(v,ei)ei
```

Every normalized spatial vector must have finite positive norm. The completed
spatial basis is made right-handed relative to the `(r,theta,phi)` coordinate
orientation by flipping only the final leg when necessary. Static construction
uses coordinate `r`, `theta`, and `phi` seeds. Look-at construction preserves
the projected look direction as leg 1 and projected up direction as leg 2;
leg 3 completes the right-handed basis. It is a deterministic first attitude
law, not transport between events.

Construction succeeds only when the final maximum orthonormality error is
strictly below the v3 ordinary-exterior gate `1e-12`. Static construction returns
`ObserverError::StaticWorldlineNotTimelike` whenever `g_tt>=0`; it does not
clamp the ergosphere. Arbitrary construction requires a unit timelike
four-velocity to the same `1e-12` gate and treats the caller's selected
four-velocity as the time orientation. A general metric has no independent
coordinate-sign rule that can prove which timelike branch the caller intended.

### L1: local physical initialization

`include/solar/relativity/local_initialization.h` defines:

```cpp
enum class InitialStateError {
    None,
    NonFiniteInput,
    InvalidObserverFrame,
    InvalidLocalDirection,
    SuperluminalLocalVelocity,
    NonFutureDirected,
    ConstraintViolation,
};

struct InitialStateResult {
    InitialStateError error = InitialStateError::None;
    std::optional<PhaseSpaceState> state;
    double measured_frequency =
        std::numeric_limits<double>::quiet_NaN();
    std::string message;

    explicit operator bool() const noexcept;
};

double observer_measured_frequency(
    const Covariant4& momentum,
    const Contravariant4& observer_velocity) noexcept;

InitialStateResult initialize_local_photon(
    const Metric& metric,
    const ObserverFrame& observer,
    const Vec3& local_direction,
    double affine = 0.0);

InitialStateResult initialize_local_timelike(
    const Metric& metric,
    const ObserverFrame& observer,
    const Vec3& local_velocity,
    double affine = 0.0);
```

A nonzero finite photon direction is normalized before forming
`k^(a)=(1,n)`. This makes the API accept a direction rather than requiring
callers to pre-normalize it. The coordinate vector is formed from the tetrad,
lowered with the metric, and scaled by the positive measured frequency.
Successful photon output must satisfy:

```text
abs(nu_obs-1) <= 1e-12
Hamiltonian normalized error <= 1e-10
```

Timelike velocity requires finite `|v|^2<1` and forms
`u^(a)=gamma(1,v)`. It is lowered without photon rescaling and must satisfy the
unit-mass Hamiltonian gate. Both paths reject a non-finite or non-orthonormal
observer frame and require positive measured frequency.

### L1: Kerr constants and L2 monitoring

`include/solar/relativity/kerr_constants.h` defines:

```cpp
struct KerrConstants {
    double E;
    double Lz;
    double Q;
    double mass_sq;
};

KerrConstants evaluate_kerr_constants(
    const KerrBoyerLindquistMetric& metric,
    const PhaseSpaceState& state,
    GeodesicKind kind);
```

It implements the v3 covariant-momentum convention:

```text
E  = -p_t
Lz =  p_phi
Q  = p_theta^2
     + cos(theta)^2 [
         a^2 (mass_sq-E^2) + Lz^2/sin(theta)^2
       ]
```

Unknown kind, non-finite state, invalid metric point, or non-finite result is
rejected. No general metric is falsely assigned a Carter invariant.

`GeodesicIntegrationConfig` gains an empty-by-default:

```cpp
using InvariantEvaluator =
    std::function<double(const PhaseSpaceState&)>;

InvariantEvaluator carter_evaluator;
```

`IntegrationDiagnostics` gains `max_carter_abs_error`. When the evaluator is
empty, both Carter fields remain NaN. When present, the initial value is
evaluated before any trajectory work, every accepted state is checked, and
diagnostics use:

```text
relative = abs(C-C0) / max(1,abs(C0))
absolute = abs(C-C0)
```

The same v3 denominator replaces the Phase 1 E/Lz drift denominator. An
evaluator exception or non-finite result terminates explicitly as a non-finite
invariant evaluation; it is never ignored.

### L1: Kerr analytic orbit quantities

`include/solar/relativity/kerr_orbits.h` defines:

```cpp
enum class OrbitSense {
    Prograde,
    Retrograde,
};

enum class CircularOrbitStability {
    Stable,
    Unstable,
};

struct CircularTimelikeOrbit {
    double radius;
    double angular_velocity;
    double specific_energy;
    double specific_lz;
    CircularOrbitStability stability;
};

struct CircularOrbitResult {
    std::optional<CircularTimelikeOrbit> orbit;
    std::string message;

    explicit operator bool() const noexcept;
};
```

Functions provide the v3 ISCO, equatorial photon-orbit, marginally-bound, and
circular timelike quantities, plus:

```cpp
ObserverResult make_equatorial_circular_observer(
    const KerrBoyerLindquistMetric& metric,
    const Contravariant4& x,
    OrbitSense sense);
```

`OrbitSense` is relative to the black-hole spin. For negative spin, prograde
coordinate `Omega` and `Lz` are negative. At zero spin, prograde uses positive
`phi` and retrograde negative `phi`; the radii remain equal. A timelike
circular orbit may be unstable. It is rejected only when the timelike
normalization is non-positive or the radius is at/below the corresponding
photon orbit. Stability is determined separately from the corresponding ISCO.

Analytic `E` and `Lz` are validated against lowering the constructed observer
four-velocity with the metric.

### L1 validation reference: asymptotic Kerr shadow

`include/solar/relativity/kerr_shadow.h` defines:

```cpp
struct ShadowCriticalPoint {
    double alpha;
    double beta;
    double photon_radius;
};

std::vector<ShadowCriticalPoint> bardeen_shadow_curve(
    const KerrBoyerLindquistMetric& metric,
    double inclination,
    std::size_t samples_per_branch);
```

For nonzero spin it solves the two inclination-dependent
`beta(r_p)^2=0` visible tips inside the equatorial spherical-photon interval,
then samples uniformly between their screen `alpha` values and solves the
corresponding spherical-photon radius. The solve uses the cancellation-safe
equation
`xi_numerator - alpha*a*sin(i)*(r_p-M)=0`; `beta^2` is evaluated as
`eta+a^2*cos(i)^2-alpha^2*cos(i)^2`, avoiding division by a vanishing
`sin(i)`. This keeps the curve closed and its interior accurate even for
near-axis views. The v3 quantities use dimensionless intermediates; only a
small negative interior radicand within the documented roundoff allowance is
clamped to zero. Inclination must be finite and strictly inside `(0,pi)`. A
numerically unresolvable sub-ULP visible interval fails explicitly.

For `|a/M| <= 64*sqrt(epsilon)`, the numerically singular rotating formula is
replaced by the Schwarzschild circle `alpha^2+beta^2=27M^2`. This branch is
documented as a double-precision limit, not new physics.

The production API is analytic only. Numerical backward-ray comparison stays
in `tests/relativity/test_kerr_shadow_raytrace.cpp`:

- construct equatorial ZAMOs at `r=1000M` and `2000M`;
- lower the observer time and azimuthal tetrad legs and solve the local
  azimuthal direction so every input screen coordinate satisfies the conserved
  asymptotic relation `alpha=-Lz/E`;
- integrate with a negative affine step;
- classify a safe exterior inner-boundary event as a capture proxy and a
  larger-radius event as escape;
- binary-search the two horizontal shadow edges;
- require both Kerr radii to meet the v3 sampled-screen gates and remain
  mutually stable;
- require the Schwarzschild CPU critical root to match `3*sqrt(3)M` with
  relative error below `1e-6`.

The test never maps `InvalidMetricPoint` to `HorizonCrossing`.

The public construction, normalized sampling, and root/formula machinery are
separated into `kerr_shadow.cpp`, private `kerr_shadow_sampling`, and private
`kerr_shadow_geometry`, respectively. This keeps numerical stabilization out
of the public API and preserves the dependency direction
`construction -> sampling -> geometry`.

## Correctness and failure gates

- All enum inputs are validated; unknown values never choose fallback physics.
- All input vectors, coordinates, radii, inclinations, and derived square-root
  arguments are checked for finite valid domains.
- Gram-Schmidt rejects degenerate seeds instead of inventing axes.
- Tetrad error is measured from all 16 `g(e_a,e_b)-eta_ab` components.
- Static observers fail inside the ergosphere.
- Circular timelike existence and stability are separate outcomes.
- Photon future direction is measured by `-p.u`, never by a coordinate-time
  sign.
- Backward rays retain future-directed momentum.
- Carter diagnostics remain opt-in and Kerr-specific through the supplied
  evaluator.
- Bardeen finite-distance disagreement is treated as convergence behavior,
  not a hard equality.

## Validation design

Focused executable tests will cover:

- tensor lowering/raising and metric contractions in Minkowski and Kerr;
- arbitrary boosted Minkowski frames and tetrad round trips;
- Kerr ZAMO orthonormality, zero angular momentum, and far-field limit;
- static observer success outside and explicit ergosphere failure;
- look/up preservation and degenerate-seed rejection;
- photon nullness, `nu_obs=1`, future direction, and local-coordinate-local
  round trip;
- timelike unit norm, Lorentz gamma behavior, and `|v|>=1` rejection;
- Kerr constants for equatorial `Q=0` and generic non-equatorial states;
- Carter relative and absolute drift on an integrated generic Kerr null ray;
- Schwarzschild and signed-spin ISCO/photon/marginally-bound limits;
- stable and unstable circular timelike orbits;
- circular analytic `E/Lz/Omega` against the metric-lowered observer;
- Schwarzschild shadow radius and Kerr critical-curve symmetry/scaling;
- CPU backward-ray horizontal edges converging toward the analytic Kerr
  boundary;
- malformed inputs and unknown enums for every new public entry point.

After focused RED/GREEN cycles, the gate is:

```bash
make clean
make
make test
git diff --check
```

All relativity tests and the shared legacy DOPRI adapter will also be rebuilt
and run under AddressSanitizer and UndefinedBehaviorSanitizer.

## Independent references

- The v3 contract remains authoritative for API scope and formulas.
- Cunha et al., *Shadows of Kerr black holes with and without scalar hair*,
  provides an independent ZAMO, local-frequency, backward-ray, and numerical
  Kerr-shadow reference: <https://arxiv.org/abs/1605.08293>.
- The Black Hole Perturbation Toolkit documents the same `(-,+,+,+)`,
  geometrized-unit convention and maintained Kerr constants/special-orbit
  capabilities: <https://bhptoolkit.org/KerrGeodesics/>.
- No external package implementation is copied into Solar.

## Compatibility

- Existing `Metric`, canonical state, Hamiltonian, event, and geodesic
  signatures remain source-compatible.
- The optional Carter evaluator and absolute diagnostic default to
  unavailable. The absolute field is appended after all pre-existing
  diagnostic members so Phase 1 aggregate initializers remain source-compatible.
- No new dependency is added.
- The existing Phase 1 branch and draft PR remain unchanged.

## Explicitly deferred

- Fermi-Walker and parallel-transported attitude laws.
- Separated Kerr/Mino-time integration and turning-point handling.
- Cartesian Kerr-Schild and physical horizon crossing.
- Ray bundles, Jacobi fields, polarization, and redshift along materials.
- Disk/material models, radiative transfer, renderer, images, movies, UI,
  WASM, GPU, and Solar Local Patch integration.
