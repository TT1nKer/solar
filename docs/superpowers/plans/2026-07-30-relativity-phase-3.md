# Solar Relativity Phase 3 Kerr Separated/Mino Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a validated, allocation-bounded Kerr separated/Mino-time CPU
geodesic path with explicit turning points and Hamiltonian worldline
cross-checks.

**Architecture:** A public L3 contract composes focused L1 potential, state,
turning, and event blocks inside one L2 integration flow. The flow evolves the
fixed state `[t,r,mu,phi,lambda]` with signed first-order equations, treats a
simple potential root as an explicit state transition, and leaves the generic
Hamiltonian BL solver authoritative for common-event validation.

**Tech Stack:** C++17, Solar fixed-size DOPRI5, existing Kerr BL/Hamiltonian
types, Make and CMake/CTest, AddressSanitizer, UndefinedBehaviorSanitizer.

## Global Constraints

- Signature is `(-,+,+,+)` and BL order is `(t,r,theta,phi)`.
- Canonical momentum is covariant `p_mu`; `E=-p_t`, `Lz=p_phi`.
- Units are geometrized `G=c=1`; `mass_M` is a length and Mino time has
  inverse-length units.
- Backward tracing uses a negative Mino step without reversing physical
  momentum.
- Never continue a production trajectory with `sqrt(max(V,0))`.
- A negative-potential trial is rejected; a simple root is located and
  crossed explicitly; a double root terminates as `NearCriticalOrbit`.
- Generic Hamiltonian BL remains the authority; comparisons use common
  physical events rather than equal integration parameters.
- Keep the hot integration loop free of per-step heap allocation.
- Existing public geodesic behavior is unchanged; append new enum values.
- Do not add Kerr-Schild, matter, transfer, renderer, CUDA, C ABI, WASM, or UI
  work in this plan.

---

## File map

### Public L0/L3 contract

- `include/solar/relativity/geodesic_types.h`: append
  `TerminationReason::NearCriticalOrbit`.
- `include/solar/relativity/kerr_separated.h`: configuration, diagnostics,
  result, and `KerrSeparatedIntegrator`.

### L1 building blocks

- `src/relativity/kerr_separated_potentials.h`
- `src/relativity/kerr_separated_potentials.cpp`
  own finite potential evaluation, derivatives, and normalized scales.
- `src/relativity/kerr_separated_state.h`
- `src/relativity/kerr_separated_state.cpp`
  own Hamiltonian-to-Mino initialization and public-state reconstruction.
- `src/relativity/kerr_separated_turning.h`
- `src/relativity/kerr_separated_turning.cpp`
  own bracketed scalar roots and simple/critical classification.
- `src/relativity/kerr_separated_events.h`
- `src/relativity/kerr_separated_events.cpp`
  own fixed-state dense event localization and public `EventHit` conversion.

### L2 flow

- `src/relativity/kerr_separated.cpp`: validates configuration and events,
  attempts steps, performs turning transitions, applies limits, and records
  diagnostics.

### Behavior tests

- `tests/relativity/test_kerr_separated_potentials.cpp`
- `tests/relativity/test_kerr_separated_state.cpp`
- `tests/relativity/test_kerr_separated.cpp`
- `tests/relativity/test_kerr_separated_crosscheck.cpp`
- `tests/relativity/test_types.cpp`
- `tests/external_consumer/main.cpp`

### Evidence

- `docs/validation/relativity_08_kerr_separated.md`
- `RELATIVITY_STATUS.md`

The existing CMake and Make source globs already discover new `.cpp` and
`test_*.cpp` files. CMake installs `include/solar`, so no build-file change is
needed for the new public header.

---

### Task 1: Public contract and literal separated potentials

**Files:**
- Modify: `include/solar/relativity/geodesic_types.h`
- Create: `include/solar/relativity/kerr_separated.h`
- Create: `src/relativity/kerr_separated_potentials.h`
- Create: `src/relativity/kerr_separated_potentials.cpp`
- Create: `tests/relativity/test_kerr_separated_potentials.cpp`
- Modify: `tests/relativity/test_types.cpp`

**Interfaces:**
- Consumes: `KerrConstants`, `GeodesicKind`,
  `numerics::Dopri5Config<5>`, `GeodesicEvent`, and `EventHit`.
- Produces:

```cpp
namespace solar::relativity::detail {

struct KerrSeparatedPotentialValues {
    double delta;
    double sigma;
    double radial;
    double polar;
    double radial_derivative;
    double polar_derivative;
    double radial_scale;
    double polar_scale;
};

class KerrSeparatedPotentials {
public:
    KerrSeparatedPotentials(
        double mass_M,
        double spin_length_M,
        KerrConstants constants);

    KerrSeparatedPotentialValues evaluate(
        double radius_M,
        double mu) const;
};

} // namespace solar::relativity::detail
```

- [ ] **Step 1: Write the failing literal-potential tests**

Use hand-derived fixtures rather than production helpers:

```cpp
const KerrConstants null_constants{1.0, 2.0, 3.0, 0.0};
const detail::KerrSeparatedPotentials null_potential(
    1.0, 0.5, null_constants);
const auto values = null_potential.evaluate(4.0, 0.25);

check_near("literal Delta", values.delta, 8.25, 1.0e-15);
check_near("literal Sigma", values.sigma, 16.015625, 1.0e-15);
check_near("literal radial R", values.radial, 189.25, 1.0e-13);
check_near("literal polar U", values.polar, 2.5771484375, 1.0e-14);
check_near(
    "literal radial derivative", values.radial_derivative,
    212.5, 1.0e-13);
check_near(
    "literal polar derivative", values.polar_derivative,
    -3.390625, 1.0e-14);
```

Add a unit-mass fixture, `a -> -a` fixture with unchanged `E,Lz,Q`, a
five-point finite-difference derivative comparison, and rejection checks for
non-positive mass, non-finite constants, `abs(mu)>1`, and overflowing
coordinates. In `test_types.cpp`, assert that the appended enum is
representable without changing aggregate defaults.

- [ ] **Step 2: Run the new test and verify RED**

Run:

```bash
make tests/relativity/test_kerr_separated_potentials
```

Expected: compilation fails because `kerr_separated.h` and the internal
potential header do not exist.

- [ ] **Step 3: Add the minimal public types and potential implementation**

Define the public configuration and results exactly as follows:

```cpp
struct KerrSeparatedConfig {
    GeodesicKind kind;
    numerics::Dopri5Config<5> dopri5;
    double initial_mino_step;
    double min_mino_step;
    double max_mino_step;
    std::size_t max_rejections_per_step;
    std::size_t max_total_steps;
    double max_affine;
    double max_coordinate_time;
    double potential_tolerance;
    double root_tolerance;
    double critical_derivative_tolerance;
    double polar_axis_tolerance;

    static KerrSeparatedConfig cpu_reference(
        GeodesicKind kind,
        double mass_scale,
        double initial_mino_step,
        double max_mino_step,
        double max_affine);
};

struct KerrSeparatedDiagnostics {
    std::size_t accepted_steps = 0;
    std::size_t rejected_steps = 0;
    std::size_t radial_turns = 0;
    std::size_t polar_turns = 0;
    double min_mino_step = std::numeric_limits<double>::quiet_NaN();
    double max_mino_step = std::numeric_limits<double>::quiet_NaN();
    double min_radius_M = std::numeric_limits<double>::quiet_NaN();
    double azimuthal_advance = 0.0;
    double winding = 0.0;
    double max_radial_residual = 0.0;
    double max_polar_residual = 0.0;
    double max_constraint_error = 0.0;
    double max_carter_rel_error = 0.0;
    TerminationReason reason = TerminationReason::NonFiniteState;
    std::string message;
};

struct KerrSeparatedIntegrationResult {
    PhaseSpaceState final_state;
    KerrConstants constants;
    KerrSeparatedDiagnostics diagnostics;
    std::optional<EventHit> event;
};

class KerrSeparatedIntegrator {
public:
    explicit KerrSeparatedIntegrator(
        const KerrBoyerLindquistMetric& metric) noexcept;

    KerrSeparatedIntegrationResult integrate(
        const PhaseSpaceState& initial,
        const KerrSeparatedConfig& config,
        const std::vector<GeodesicEvent>& events = {}) const;

private:
    const KerrBoyerLindquistMetric* metric_;
};
```

The public class declaration matches the design spec. Append
`NearCriticalOrbit` after `UserEvent`.

Implement the equations from the design literally. Compute scales as:

```cpp
radial_scale = std::max(
    {std::pow(mass_M, 4), std::fabs(radial)});
polar_scale = std::max(
    {mass_M * mass_M, std::fabs(polar)});
```

Throw `std::invalid_argument` for invalid constructor inputs and
`std::domain_error` for an invalid evaluation point or non-finite result.

- [ ] **Step 4: Run the focused tests and verify GREEN**

Run:

```bash
make tests/relativity/test_kerr_separated_potentials \
     tests/relativity/test_types
./tests/relativity/test_kerr_separated_potentials
./tests/relativity/test_types
```

Expected: all assertions pass; finite-difference normalized derivative error
is below `1e-8`.

- [ ] **Step 5: Commit the potential contract**

```bash
git add include/solar/relativity/geodesic_types.h \
        include/solar/relativity/kerr_separated.h \
        src/relativity/kerr_separated_potentials.h \
        src/relativity/kerr_separated_potentials.cpp \
        tests/relativity/test_kerr_separated_potentials.cpp \
        tests/relativity/test_types.cpp
git commit -m "feat(relativity): add separated Kerr potentials"
```

---

### Task 2: Hamiltonian/Mino state conversion

**Files:**
- Create: `src/relativity/kerr_separated_state.h`
- Create: `src/relativity/kerr_separated_state.cpp`
- Create: `tests/relativity/test_kerr_separated_state.cpp`

**Interfaces:**
- Consumes: `HamiltonGeodesicRhs`, `evaluate_kerr_constants`, and
  `KerrSeparatedPotentials`.
- Produces:

```cpp
namespace solar::relativity::detail {

using KerrMinoState = numerics::StateN<5>;

enum class SeparatedDirection : int {
    Negative = -1,
    Locked = 0,
    Positive = 1,
};

struct KerrSeparatedState {
    KerrMinoState values;
    SeparatedDirection radial_direction;
    SeparatedDirection polar_direction;
};

struct KerrSeparatedInitialState {
    KerrSeparatedState state;
    KerrConstants constants;
};

class KerrSeparatedCriticalInitialState : public std::domain_error {
public:
    using std::domain_error::domain_error;
};

KerrSeparatedInitialState initialize_kerr_separated_state(
    const KerrBoyerLindquistMetric& metric,
    const PhaseSpaceState& initial,
    GeodesicKind kind,
    double normalized_potential_tolerance,
    double normalized_critical_tolerance,
    double integration_direction);

PhaseSpaceState reconstruct_kerr_phase_space(
    const KerrBoyerLindquistMetric& metric,
    const KerrConstants& constants,
    const KerrSeparatedState& state);

} // namespace solar::relativity::detail
```

State indices are constants named `kTime=0`, `kRadius=1`, `kMu=2`,
`kAzimuth=3`, and `kAffine=4` in the internal header.

- [ ] **Step 1: Write round-trip and direction tests**

Construct literal canonical states whose `p_r` and `p_theta` are nonzero.
Assert:

```cpp
check_near("round-trip affine", rebuilt.affine, initial.affine, 1.0e-14);
for (std::size_t i = 0; i < 4; ++i) {
    check_near("round-trip coordinate", rebuilt.x.v[i],
               initial.x.v[i], 1.0e-12);
    check_near("round-trip momentum", rebuilt.p.v[i],
               initial.p.v[i], 1.0e-11);
}
check("inward radial direction",
      separated.state.radial_direction ==
          detail::SeparatedDirection::Negative);
check("mu direction follows -sin(theta) theta-dot",
      separated.state.polar_direction ==
          detail::SeparatedDirection::Negative);
```

Add:

- exact equatorial `theta=pi/2,p_theta=0` becomes `Locked`;
- a simple radial-root initial state selects the departure sign from
  `sign(integration_direction * R')`;
- a radial double root throws a typed internal
  `KerrSeparatedCriticalInitialState`;
- inconsistent constants/potential beyond tolerance are rejected;
- near-axis nonzero-`Lz` reconstruction is rejected.

- [ ] **Step 2: Run the state test and verify RED**

Run:

```bash
make tests/relativity/test_kerr_separated_state
```

Expected: compilation fails because `kerr_separated_state.h` does not exist.

- [ ] **Step 3: Implement conversion without fallback signs**

Compute the Hamiltonian tangent once:

```cpp
const auto rhs = HamiltonGeodesicRhs(metric)(initial);
const double sigma = potential_values.sigma;
const double radial_velocity = sigma * rhs.dx.v[1];
const double polar_velocity =
    -std::sin(initial.x.v[2]) * sigma * rhs.dx.v[2];
```

Use a nonzero finite velocity sign directly. At a potential root with zero
velocity, use `sign(integration_direction * V')` only when the normalized
derivative exceeds the critical tolerance. Preserve exact equatorial locked
motion. Throw the critical marker for any other `V≈0,V'≈0` initial state.

Reconstruct `theta=acos(mu)` and:

```cpp
p_t = -E;
p_phi = Lz;
p_r = radial_direction * sqrt(R) / Delta;
p_theta = -polar_direction * sqrt(U) / sqrt(1-mu*mu);
```

For a locked direction, set the corresponding square-root momentum to zero.
Reject a negative potential instead of clamping it.

- [ ] **Step 4: Run focused tests and verify GREEN**

Run:

```bash
make tests/relativity/test_kerr_separated_state
./tests/relativity/test_kerr_separated_state
./tests/relativity/test_kerr_separated_potentials
```

Expected: all assertions pass and round-trip normalized Hamiltonian error is
below `1e-12`.

- [ ] **Step 5: Commit state conversion**

```bash
git add src/relativity/kerr_separated_state.h \
        src/relativity/kerr_separated_state.cpp \
        tests/relativity/test_kerr_separated_state.cpp
git commit -m "feat(relativity): convert Kerr states to Mino form"
```

---

### Task 3: Basic separated flow, limits, and terminal events

**Files:**
- Create: `src/relativity/kerr_separated_events.h`
- Create: `src/relativity/kerr_separated_events.cpp`
- Create: `src/relativity/kerr_separated.cpp`
- Create: `tests/relativity/test_kerr_separated.cpp`

**Interfaces:**
- Consumes: public Phase 3 API, state conversion, potentials, DOPRI5, and
  existing event contracts.
- Produces:

```cpp
namespace solar::relativity::detail {

struct KerrSeparatedEventSelection {
    enum class Status { None, Found, Failed } status;
    std::optional<EventHit> hit;
    TerminationReason reason;
    std::string message;
};

KerrSeparatedEventSelection select_initial_kerr_event(
    const PhaseSpaceState& initial,
    const std::vector<GeodesicEvent>& events);

KerrSeparatedEventSelection select_first_kerr_step_event(
    const KerrBoyerLindquistMetric& metric,
    const KerrConstants& constants,
    const KerrSeparatedState& directions,
    const numerics::Dopri5DenseOutput<5>& dense,
    const std::vector<GeodesicEvent>& events);

} // namespace solar::relativity::detail
```

- [ ] **Step 1: Write basic integration and event tests**

Use an equatorial outgoing Schwarzschild null state and assert:

```cpp
const auto result = integrator.integrate(initial, config, {escape_event});
check("escape event reason",
      result.diagnostics.reason == TerminationReason::Escaped);
check("escape event exists", result.event.has_value());
check_near("escape radius", result.final_state.x.v[1], 20.0, 1.0e-9);
check("no radial turn", result.diagnostics.radial_turns == 0);
check("equatorial plane remains locked",
      std::fabs(result.final_state.x.v[2] - pi/2.0) < 1.0e-13);
```

Add tests for:

- initial any-direction event accepts zero steps;
- two roots in one accepted step select the first in integration direction;
- directed event does not fire at the initial point;
- `MaxAffine`, `MaxCoordinateTime`, and `MaxSteps`;
- negative Mino integration decreases affine while momentum remains
  future-directed;
- malformed event and unknown enum fail before a step;
- zero/non-finite/wrong-sign steps and non-positive tolerances fail
  deterministically.

- [ ] **Step 2: Run the flow test and verify RED**

Run:

```bash
make tests/relativity/test_kerr_separated
```

Expected: link failure because `KerrSeparatedIntegrator::integrate` and the
factory are not defined.

- [ ] **Step 3: Implement event localization and the no-turn flow**

The event locator evaluates start/end public states, checks direction,
brackets the scalar event value, and performs at most 100 safeguarded
secant/bisection iterations over the dense Mino interval. It terminates when:

```cpp
std::fabs(candidate_state.affine - bracket_state.affine)
    <= event.root_tolerance
```

The RHS computes the five equations from the design. If `R<0`, `U<0`, the BL
point is invalid, or a result is non-finite, it returns NaNs and records the
specific stage failure. In this task, a forbidden potential trial is rejected
and shrunk; exhausting the rejection budget returns `StepUnderflow`.

Use these factory defaults:

```cpp
relative_tolerance = 1.0e-11;
absolute_tolerance = {
    1.0e-12 * mass_scale,
    1.0e-12 * mass_scale,
    1.0e-13,
    1.0e-13,
    1.0e-12 * mass_scale,
};
min_mino_step = 1.0e-14 / mass_scale;
max_rejections_per_step = 32;
max_total_steps = 1'000'000;
max_coordinate_time = infinity;
potential_tolerance = 1.0e-12;
root_tolerance = 1.0e-12;
critical_derivative_tolerance = 1.0e-10;
polar_axis_tolerance = 1.0e-12;
```

Validate every public field before initialization. Apply limits in integration
direction, preserve the sign of proposed steps, and reject a step that cannot
advance the Mino parameter.

- [ ] **Step 4: Run basic flow tests and verify GREEN**

Run:

```bash
make tests/relativity/test_kerr_separated
./tests/relativity/test_kerr_separated
./tests/relativity/test_geodesic_events
./tests/relativity/test_geodesics
```

Expected: all assertions pass and existing generic event behavior is
unchanged.

- [ ] **Step 5: Commit the basic flow**

```bash
git add src/relativity/kerr_separated_events.h \
        src/relativity/kerr_separated_events.cpp \
        src/relativity/kerr_separated.cpp \
        tests/relativity/test_kerr_separated.cpp
git commit -m "feat(relativity): integrate separated Kerr rays"
```

---

### Task 4: Simple turning transitions and critical roots

**Files:**
- Create: `src/relativity/kerr_separated_turning.h`
- Create: `src/relativity/kerr_separated_turning.cpp`
- Modify: `src/relativity/kerr_separated.cpp`
- Modify: `tests/relativity/test_kerr_separated.cpp`

**Interfaces:**
- Consumes: scalar `R/R'` or `U/U'` evaluations and the accepted pre-turn
  separated state.
- Produces:

```cpp
namespace solar::relativity::detail {

enum class TurningCoordinate { Radial, Polar };
enum class TurningStatus { Simple, NearCritical, Failed };

struct TurningRoot {
    TurningStatus status;
    double coordinate;
    double normalized_derivative;
    std::size_t iterations;
    std::string message;
};

TurningRoot locate_kerr_turning_root(
    TurningCoordinate coordinate,
    double allowed_coordinate,
    double forbidden_coordinate,
    const KerrSeparatedPotentials& potentials,
    double fixed_other_coordinate,
    double normalized_root_tolerance,
    double normalized_potential_tolerance,
    double normalized_critical_derivative_tolerance);

} // namespace solar::relativity::detail
```

- [ ] **Step 1: Add failing radial, polar, and critical tests**

Use literal canonical fixtures and assert:

```cpp
check("scattering ray returns through outer sphere",
      result.diagnostics.reason == TerminationReason::Escaped);
check("one radial turn", result.diagnostics.radial_turns == 1);
check("minimum radius is finite",
      std::isfinite(result.diagnostics.min_radius_M));

check("polar return event reached",
      polar.diagnostics.reason == TerminationReason::UserEvent);
check("one polar turn", polar.diagnostics.polar_turns == 1);

check("spherical photon is explicitly critical",
      critical.diagnostics.reason ==
          TerminationReason::NearCriticalOrbit);
check("critical ray does not count a simple radial turn",
      critical.diagnostics.radial_turns == 0);
```

Also assert that halving `max_mino_step` changes the return-crossing state by
less than the coarser run's error, and that no accepted state reports negative
normalized potential beyond `potential_tolerance`.

- [ ] **Step 2: Run the turn tests and verify RED**

Run:

```bash
make tests/relativity/test_kerr_separated
./tests/relativity/test_kerr_separated
```

Expected: the scattering and polar cases exhaust rejection/step limits near
their first turning point; the spherical case is not yet classified.

- [ ] **Step 3: Implement bracketed roots and bounded release**

`locate_kerr_turning_root` requires one allowed and one forbidden potential
value, then uses safeguarded secant/bisection for at most 100 iterations.
Radial bracket width is normalized by `M`; polar width is already
dimensionless. Classify `NearCritical` only when both normalized potential and
derivative are within their configured tolerances.

When a trial first records a forbidden stage:

1. locate the scalar root between the accepted coordinate and the first
   forbidden coordinate;
2. shrink and retry until the accepted coordinate is within
   `root_tolerance` of the root;
3. advance the remaining complete state with a bounded local root step using
   `q''=V'/2` and midpoint values for `t`, `phi`, and `lambda`;
4. set the root coordinate exactly, flip only the affected direction, and
   increment the matching turn count;
5. release with:

```cpp
const double coordinate_tolerance =
    coordinate == TurningCoordinate::Radial
        ? config.root_tolerance * metric.mass()
        : config.root_tolerance;
const double release_h = std::copysign(
    std::min(
        config.max_mino_step,
        std::max(
            config.min_mino_step,
            2.0 * std::sqrt(
                coordinate_tolerance /
                std::fabs(physical_derivative)))),
    attempted_step);
const double released_coordinate =
    root + 0.25 * physical_derivative * release_h * release_h;
```

Cap `release_h` by `max_mino_step` and current affine/time limits. Accept the
release only if its normalized potential is nonnegative, set the outgoing
direction to `sign(physical_derivative * release_h)`, and require the
reconstructed Hamiltonian state to be finite. A critical root returns
`NearCriticalOrbit`; an unresolved bracket returns `EventRootFailure`; an
unrepresentable release returns `StepUnderflow`.

- [ ] **Step 4: Run turn and regression tests and verify GREEN**

Run:

```bash
make tests/relativity/test_kerr_separated
./tests/relativity/test_kerr_separated
./tests/relativity/test_geodesics_kerr
./tests/relativity/test_kerr_shadow_raytrace
```

Expected: all tests pass; simple-turn counters are exact and the existing
Hamiltonian shadow gate remains unchanged.

- [ ] **Step 5: Commit turning support**

```bash
git add src/relativity/kerr_separated_turning.h \
        src/relativity/kerr_separated_turning.cpp \
        src/relativity/kerr_separated.cpp \
        tests/relativity/test_kerr_separated.cpp
git commit -m "feat(relativity): handle Kerr turning points"
```

---

### Task 5: Render-relevant invariant diagnostics and failure semantics

**Files:**
- Modify: `src/relativity/kerr_separated.cpp`
- Modify: `tests/relativity/test_kerr_separated.cpp`

**Interfaces:**
- Consumes: reconstructed public states, `hamiltonian_constraint_error`,
  `evaluate_kerr_constants`, and separated potential values.
- Produces: every `KerrSeparatedDiagnostics` field with explicit availability
  and no silent repair.

- [ ] **Step 1: Add failing diagnostic and failure tests**

For a successful winding ray, assert:

```cpp
check_near(
    "winding derives from unwrapped azimuth",
    result.diagnostics.winding,
    result.diagnostics.azimuthal_advance / (2.0 * pi),
    1.0e-15);
check("minimum radius includes initial and accepted states",
      result.diagnostics.min_radius_M <= initial.x.v[1]);
check("Hamiltonian diagnostic meets reference gate",
      result.diagnostics.max_constraint_error < 1.0e-10);
check("Carter diagnostic meets reference gate",
      result.diagnostics.max_carter_rel_error < 1.0e-10);
```

Add explicit cases for:

- non-finite initial state -> `NonFiniteState`;
- BL-invalid initial point -> `InvalidMetricPoint`;
- near-axis nonzero `Lz` -> `InvalidMetricPoint` with zero accepted steps;
- event callback throw/non-finite return -> `EventRootFailure`;
- stage non-finite unrelated to a potential root -> `NonFiniteState`;
- rejection budget and total-step budget report distinct reasons;
- `min_mino_step` and `max_mino_step` remain NaN when no step is accepted.

- [ ] **Step 2: Run the focused test and verify RED**

Run:

```bash
make tests/relativity/test_kerr_separated
./tests/relativity/test_kerr_separated
```

Expected: at least the constraint, Carter, winding, and no-step availability
assertions fail because the diagnostics are not yet updated on every accepted
transition.

- [ ] **Step 3: Centralize diagnostic updates**

After initial validation and each accepted ordinary/root/release state:

```cpp
diagnostics.min_radius_M =
    std::min(diagnostics.min_radius_M, state.x.v[1]);
diagnostics.azimuthal_advance =
    state.x.v[3] - initial.x.v[3];
diagnostics.winding =
    diagnostics.azimuthal_advance / (2.0 * pi);
diagnostics.max_constraint_error = std::max(
    diagnostics.max_constraint_error,
    hamiltonian_constraint_error(metric, state, config.kind));
```

Compute Carter drift with denominator
`max(abs(Q0), mass_M*mass_M*1e-14)` and retain the absolute numerator for a
zero-`Q` internal check. Compute radial/polar residuals from reconstructed
Mino velocities and the independently evaluated potentials. If any diagnostic
evaluation is non-finite, terminate as `NonFiniteState`; do not preserve a
successful event classification.

- [ ] **Step 4: Run focused and failure regressions and verify GREEN**

Run:

```bash
make tests/relativity/test_kerr_separated
./tests/relativity/test_kerr_separated
./tests/relativity/test_geodesic_failures
./tests/relativity/test_kerr_constants
```

Expected: all assertions pass and failed paths retain their original state and
zero accepted-step count when failure precedes integration.

- [ ] **Step 5: Commit diagnostics**

```bash
git add src/relativity/kerr_separated.cpp \
        tests/relativity/test_kerr_separated.cpp
git commit -m "feat(relativity): report separated ray diagnostics"
```

---

### Task 6: Common-event Hamiltonian cross-validation

**Files:**
- Create: `tests/relativity/test_kerr_separated_crosscheck.cpp`
- Modify only if a failing fixture exposes a real defect:
  `src/relativity/kerr_separated.cpp`,
  `src/relativity/kerr_separated_events.cpp`,
  `src/relativity/kerr_separated_turning.cpp`

**Interfaces:**
- Consumes: `KerrSeparatedIntegrator`, `GeodesicIntegrator`, identical
  canonical initial states, and identical public events.
- Produces: independent worldline and convergence evidence for null/timelike,
  spin sign, radial turn, and polar turn families.

- [ ] **Step 1: Write the failing cross-check executable**

Build literal initial states by fixing `(E,Lz,p_theta)` and solving only the
Hamiltonian quadratic for `p_r`; do not derive expected crossings from the
separated equations. For each solver, terminate at the same event.

Measure coordinate disagreement with wrapped azimuth:

```cpp
const double dphi = std::remainder(
    separated.x.v[3] - hamiltonian.x.v[3], 2.0 * pi);
const double position_error_M = std::sqrt(
    std::pow((separated.x.v[1] - hamiltonian.x.v[1]) / mass_M, 2) +
    std::pow(separated.x.v[2] - hamiltonian.x.v[2], 2) +
    dphi * dphi);
```

Include:

- ordinary null Kerr `chi=+0.5`;
- ordinary null Kerr `chi=-0.5`;
- ordinary timelike Kerr `chi=+0.5`;
- Schwarzschild null;
- one-return radial scattering null;
- one-return off-equatorial polar null.

Assert every ordinary error `<1e-8`, maximum error `<1e-7`, each solver's
constraint `<1e-10`, and each separated Carter drift `<1e-10`. Run each
ordinary fixture with `max_mino_step` and half that value; the fine result
must not be farther from Hamiltonian by more than `10% + 1e-11`.

- [ ] **Step 2: Verify the cross-check catches a controlled mutation**

Temporarily add `1.0e-4` to the separated radial coordinate inside the local
`position_error_M` calculation, build, and run:

```bash
make tests/relativity/test_kerr_separated_crosscheck
./tests/relativity/test_kerr_separated_crosscheck
```

Expected: the `<1e-8` common-event assertion fails. Revert only that test
mutation and rerun. This proves the test detects a materially wrong worldline
even if the unmutated production implementation already satisfies the gate.

- [ ] **Step 3: Run the unmutated cross-check and correct demonstrated defects**

For each failing fixture, record whether the error comes from event root
conversion, ordinary DOPRI stepping, root transition, or state
reconstruction. Make the smallest corresponding production change. Do not
loosen the `1e-8/1e-7/1e-10` gates and do not special-case fixture values.
If every unmutated fixture passes, make no production change in this task.

- [ ] **Step 4: Run all Phase 3 and authority tests and verify GREEN**

Run:

```bash
make tests/relativity/test_kerr_separated_potentials \
     tests/relativity/test_kerr_separated_state \
     tests/relativity/test_kerr_separated \
     tests/relativity/test_kerr_separated_crosscheck
./tests/relativity/test_kerr_separated_potentials
./tests/relativity/test_kerr_separated_state
./tests/relativity/test_kerr_separated
./tests/relativity/test_kerr_separated_crosscheck
./tests/relativity/test_geodesics_kerr
./tests/relativity/test_kerr_shadow_raytrace
```

Expected: all assertions pass. Print compact maximum/P95 position,
constraint, Carter, and convergence errors for the validation record.

- [ ] **Step 5: Commit cross-validation**

```bash
git add tests/relativity/test_kerr_separated_crosscheck.cpp \
        src/relativity/kerr_separated.cpp \
        src/relativity/kerr_separated_events.cpp \
        src/relativity/kerr_separated_turning.cpp
git commit -m "test(relativity): cross-check separated Kerr worldlines"
```

Omit unchanged production paths from `git add`.

---

### Task 7: Installed consumer and Phase 3 evidence gate

**Files:**
- Modify: `tests/external_consumer/main.cpp`
- Create: `docs/validation/relativity_08_kerr_separated.md`
- Modify: `RELATIVITY_STATUS.md`

**Interfaces:**
- Consumes: installed `Solar::Relativity` and the final focused test output.
- Produces: proof that the new API is public/installable and a truthful
  Phase 3 status handoff.

- [ ] **Step 1: Add an installed-consumer call before changing docs**

In the external consumer, construct the ray through already installed public
observer/initialization APIs:

```cpp
KerrBoyerLindquistMetric metric(1.0, 0.5);
Contravariant4 x{{0.0, 20.0, pi / 2.0, 0.0}};
const auto observer = make_zamo_observer(metric, x);
const auto photon = initialize_local_photon(
    metric, *observer.frame, Vec3{{-1.0, 0.0, 0.0}});
const auto config = KerrSeparatedConfig::cpu_reference(
    GeodesicKind::Null, 1.0, 1.0e-5, 1.0e-4, 0.1);
const auto result = KerrSeparatedIntegrator(metric).integrate(
    *photon.state, config);
```

Return failure unless:

```cpp
result.diagnostics.reason == TerminationReason::MaxAffine &&
result.diagnostics.accepted_steps > 0 &&
result.diagnostics.max_constraint_error < 1.0e-10
```

- [ ] **Step 2: Run the installed consumer and verify RED if export is broken**

Run:

```bash
make test-external-consumer
```

Expected before an export defect fix: compile/link failure naming the missing
installed header or symbol. If it passes immediately, the existing recursive
header install and source glob are proven sufficient; do not edit CMake.

- [ ] **Step 3: Run the complete local release and sanitizer gates**

Run:

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

Compile and run the four Phase 3 tests against the sanitizer library, then
run:

```bash
git diff --check
git status --short
```

Expected: all fixture-independent tests pass, the optional DE440 fixture is
the only allowed skip, sanitizers report no issue, and the diff check is
clean.

- [ ] **Step 4: Write evidence and advance the status gate**

`docs/validation/relativity_08_kerr_separated.md` must contain:

- equations and conventions;
- exact tested fixture table;
- common-event error, P95/max, constraint, Carter, and convergence numbers;
- simple-turn and critical-root evidence;
- commands, platform, compiler, and verified commit;
- explicit limitations: BL exterior only, no physical horizon crossing,
  no analytic elliptic backend, and no near-axis nonzero-`Lz`.

Update `RELATIVITY_STATUS.md` to:

```text
CURRENT_PHASE: 3
PHASE_STATE: PASSED
NEXT_ALLOWED_ACTION:
- Phase 4 only: Kerr-Schild coordinates, transforms, and physical horizon
  crossing validation.
```

Update completed work, fastest falsification, likely bugs, assertion counts,
and verified commands from actual output. Do not claim Linux, RTX 3080, or
DE440 validation.

- [ ] **Step 5: Commit the verified Phase 3 gate**

```bash
git add tests/external_consumer/main.cpp \
        docs/validation/relativity_08_kerr_separated.md \
        RELATIVITY_STATUS.md
git commit -m "docs(relativity): record Phase 3 validation"
```

---

## Final branch review

- [ ] Compare `git diff --stat origin/main...HEAD` against the file map.
- [ ] Confirm no Phase 4+, Gargantua, renderer, CUDA, or unrelated classic
  Solar files changed.
- [ ] Re-run all four focused tests and `make test-external-consumer` after
  the final commit.
- [ ] Use `superpowers:requesting-code-review`, then
  `superpowers:verification-before-completion`.
- [ ] Publish `codex/relativity-phase-3` and open one Solar PR targeting
  `main`; do not merge the superseded PR #2.
- [ ] Only after Solar Phase 3 is merged, update Gargantua Studio's pinned
  Solar commit and add a consumer regression in its own repository/PR.
