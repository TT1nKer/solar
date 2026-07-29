# Solar Relativity Phase 1 Hamiltonian Geodesic Design

**Status:** Approved for implementation on 2026-07-29  
**Source contract:** `SOLAR_RELATIVITY_KERR_完整实施主提示词_v3.md`  
**Previous gate:** Phase 0B `PASSED`

## Goal

Build the authoritative CPU foundation for fixed-background Hamiltonian
geodesics:

- evaluate the canonical Hamiltonian and its normalized constraint error;
- evaluate the general eight-dimensional Hamilton equations for any `Metric`;
- adapt Solar's DOPRI5 implementation into a shared fixed-dimension-capable
  kernel with component tolerances and high-order dense output;
- locate directed events inside accepted steps with a bracketed safeguarded
  root solver;
- integrate null and unit-mass timelike canonical states with explicit,
  truthful termination diagnostics;
- validate Minkowski and Schwarzschild geodesics plus stationary/axisymmetric
  behavior without entering the observer/tetrad phase.

Phase 1 is a library and validation phase. It does not add a ray CLI because
the v3 ray CLI requires a physical observer and screen-to-tetrad
initialization, both owned by Phase 2.

## Governing conventions

- Signature is `(-,+,+,+)`.
- Coordinates are contravariant `x^mu`; canonical momenta are covariant
  `p_mu`.
- Coordinate zero is time.
- Units are geometrized, `G=c=1`.
- The authoritative photon momentum is future-directed physical momentum.
  Backward tracing later uses a negative affine step or the negated RHS; it
  does not replace the state with a past-directed photon.
- The state integrated by DOPRI5 is
  `Y=(x^0,x^1,x^2,x^3,p_0,p_1,p_2,p_3)`.
- The affine parameter is the independent variable and remains in
  `PhaseSpaceState::affine`, outside the eight-component array.
- Null target Hamiltonian is `0`; unit-mass timelike target is `-1/2`.
- No default per-step constraint projection is permitted.

## Scope resolution for v3 section 7.5

Section 7.5 describes local tetrad initialization, while the `Tetrad` type,
observer implementations, measured-frequency normalization, and local
round-trip tests are defined in Phase 2. Phase 1 therefore:

- accepts canonical coordinate states supplied by tests or later callers;
- validates finite state, valid metric point, and Hamiltonian constraint;
- does not define a duplicate or hidden tetrad type;
- does not claim future direction from a coordinate-time sign, which is not a
  general observer-independent test;
- records tetrad-to-coordinate conversion, observer frequency normalization,
  and local photon/timelike initialization as explicitly not completed until
  Phase 2.

This preserves the phase gate instead of either fabricating observer physics
or duplicating the Phase 2 API.

## Approaches considered

### Keep the existing dynamic-vector DOPRI5 unchanged

This would route the eight-dimensional geodesic through
`std::vector<double>`. It has the smallest initial diff but allocates at every
stage, exposes only one absolute tolerance, uses a maximum component norm
instead of the required RMS norm, has no dense output, and cannot carry the
required controller configuration. It is not suitable as the authoritative
GR path.

### Add a separate relativity-only DOPRI5

This would isolate Phase 1 from legacy Solar code, but it would duplicate the
Butcher tableau, error estimator, and controller. The two implementations
would drift and contradict the v3 instruction to generalize the audited
integrator rather than copy an unrelated solver.

### Shared container-generic stage engine

This is the selected approach. One template implementation performs all seven
Dormand–Prince stages for array-like containers. The new fixed-size public API
uses `std::array<double,N>`. The existing `dopri5_generic_step` becomes a
compatibility adapter over the same stage engine with its existing scalar
tolerance and maximum-norm policy. The specialized N-body `dopri5_step`
remains unchanged because its `State`/`Vec3` acceleration interface and
position/velocity norm semantics are a separate established contract.

This removes duplicated generic stage logic while avoiding a rewrite of the
unrelated N-body integrator.

## Layer and dependency design

### L0: geodesic contracts

`include/solar/relativity/geodesic_types.h` owns:

```cpp
enum class EventDirection {
    Any,
    Increasing,
    Decreasing,
};

enum class TerminationReason {
    HorizonCrossing,
    Escaped,
    DiskSurfaceHit,
    MaterialSurfaceHit,
    RadialTurningPoint,
    PolarTurningPoint,
    MaxAffine,
    MaxProperTime,
    MaxCoordinateTime,
    MaxSteps,
    StepUnderflow,
    InvalidMetricPoint,
    NonFiniteState,
    ConstraintViolation,
    EventRootFailure,
    UserEvent,
};

struct IntegrationDiagnostics {
    std::size_t accepted_steps = 0;
    std::size_t rejected_steps = 0;
    double min_step;
    double max_step;
    double max_constraint_error = 0.0;
    double max_energy_rel_error;
    double max_lz_rel_error;
    double max_carter_rel_error;
    TerminationReason reason;
    std::string message;
};
```

`min_step` and `max_step` are quiet NaN until an accepted step exists.
Unavailable invariant diagnostics are quiet NaN, never fabricated as zero.
Carter monitoring remains unavailable in Phase 1.

The same header defines:

```cpp
using EventFunction =
    std::function<double(const PhaseSpaceState&)>;

struct GeodesicEvent {
    std::string name;
    EventFunction function;
    EventDirection direction = EventDirection::Any;
    TerminationReason reason = TerminationReason::UserEvent;
    double root_tolerance = 1.0e-10;
};

struct EventHit {
    std::size_t event_index;
    double affine;
    PhaseSpaceState state;
    double value;
    std::size_t root_iterations;
};
```

Event directions are measured from the accepted step's start toward its end.
This definition is independent of whether the affine step is positive or
negative.

### L1: shared DOPRI5 numerical building block

`include/solar/numerics/dopri5.h` owns:

```cpp
template <std::size_t N>
using StateN = std::array<double, N>;

enum class ErrorNorm {
    RootMeanSquare,
    Maximum,
};

template <std::size_t N>
struct Dopri5Config {
    StateN<N> absolute_tolerance;
    double relative_tolerance;
    double safety;
    double min_factor;
    double max_factor;
    ErrorNorm error_norm;
};

template <std::size_t N>
class Dopri5DenseOutput {
public:
    double start() const noexcept;
    double end() const noexcept;
    StateN<N> evaluate(double independent_variable) const;
};

template <std::size_t N>
struct Dopri5StepResult {
    enum class Status {
        Completed,
        NonFiniteState,
        NonFiniteDerivative,
    };

    Status status;
    StateN<N> state;
    double step_used;
    double next_step;
    double error;
    bool accepted;
    std::optional<Dopri5DenseOutput<N>> dense_output;
};

template <std::size_t N, typename Rhs>
Dopri5StepResult<N> dopri5_step(
    const StateN<N>& state,
    double independent_variable,
    double step,
    const Rhs& rhs,
    const Dopri5Config<N>& config);
```

The implementation requirements are:

- the standard seven-stage Dormand–Prince 5(4) tableau;
- fifth-order trial state and embedded fourth-order error;
- v3 scaling
  `s_i=atol_i+rtol*max(abs(y_i),abs(y_trial_i))`;
- RMS error for the relativity path;
- acceptance only for finite `error<=1`;
- sign-preserving next step;
- explicit `safety`, `min_factor`, and `max_factor`;
- rejected-step shrink never increases the step magnitude;
- all tolerances, state values, stage derivatives, step values, and controller
  values checked for finite/valid ranges;
- invalid configuration throws `std::invalid_argument`, while a non-finite
  state or stage derivative returns the corresponding non-success status;
- a fourth-order Dormand–Prince continuous extension over the complete
  accepted step, including both endpoints;
- no extrapolation outside the step interval.

`dense_output` is present only when all seven stages and the trial state are
finite. The geodesic flow uses it only for an accepted step.

The internal stage evaluator is container-generic. The legacy
`dopri5_generic_step` adapter uses a vector state, uniform absolute tolerance,
and `ErrorNorm::Maximum` to preserve its public behavior. Its signatures in
`include/solar/integrator.h` remain unchanged.

### L1: Hamiltonian physics

`include/solar/relativity/hamiltonian.h` and
`src/relativity/hamiltonian.cpp` own:

```cpp
struct PhaseSpaceDerivative {
    Contravariant4 dx;
    Covariant4 dp;
};

double hamiltonian(
    const Metric& metric,
    const PhaseSpaceState& state);

double hamiltonian_constraint_error(
    const Metric& metric,
    const PhaseSpaceState& state,
    GeodesicKind kind);

class HamiltonGeodesicRhs {
public:
    explicit HamiltonGeodesicRhs(const Metric& metric);

    PhaseSpaceDerivative operator()(
        const PhaseSpaceState& state) const;
};

numerics::StateN<8> pack_phase_space(
    const PhaseSpaceState& state);

PhaseSpaceState unpack_phase_space(
    double affine,
    const numerics::StateN<8>& packed);
```

The implementation follows:

```text
H       = 1/2 g^mu,nu p_mu p_nu
dx^mu   = g^mu,nu p_nu
dp_mu   = -1/2 (partial_mu g^alpha,beta) p_alpha p_beta
```

The normalized constraint error is exactly:

```text
abs(H-H0) /
(1 + abs(H0) + 1/2 sum_mu,nu abs(g^mu,nu p_mu p_nu))
```

The RHS always uses the metric's general derivative interface. It does not
force `dp_t` or `dp_phi` to zero. Existing stationary/axisymmetric metrics
produce exact zero through their exact derivative dimensions. A test-only
time-dependent metric proves that the general path can produce nonzero
`dp_t`.

### L1: bracketed event root

`include/solar/relativity/event_root.h` and
`src/relativity/event_root.cpp` own a one-step root locator:

```cpp
enum class EventRootStatus {
    NoRoot,
    Found,
    Failed,
};

struct EventRootResult {
    EventRootStatus status;
    std::optional<EventHit> hit;
    std::string message;
};

EventRootResult locate_event(
    std::size_t event_index,
    const GeodesicEvent& event,
    const numerics::Dopri5DenseOutput<8>& dense_output);
```

The root locator:

- evaluates the event on `PhaseSpaceState` reconstructed at dense points;
- accepts exact roots at either endpoint;
- requires a sign bracket for an interior root;
- applies `Any`, `Increasing`, or `Decreasing` along step progression;
- combines a secant proposal with a bisection safeguard;
- terminates when affine bracket width is at most `root_tolerance`;
- rejects a non-finite event value;
- has a finite iteration cap and reports failure rather than using pure
  Newton iteration.

`NoRoot` and `Failed` are distinct. A direction mismatch or unbracketed event
is `NoRoot`; a non-finite event value or failure to reach the requested
tolerance is `Failed`.

When several events occur in one accepted step, the flow selects the smallest
dense fraction in `[0,1]`, meaning the first event encountered in the actual
integration direction.

### L2: geodesic integration flow

`include/solar/relativity/geodesic_integrator.h` and
`src/relativity/geodesic_integrator.cpp` own:

```cpp
struct GeodesicIntegrationConfig {
    GeodesicKind kind;
    numerics::Dopri5Config<8> dopri5;
    double initial_step;
    double min_step;
    double max_step;
    std::size_t max_rejections_per_step;
    std::size_t max_total_steps;
    double max_affine;
    double max_proper_time;
    double max_coordinate_time;
    double constraint_tolerance;
    bool monitor_energy;
    bool monitor_lz;

    static GeodesicIntegrationConfig cpu_reference(
        GeodesicKind kind,
        double mass_scale,
        double initial_step,
        double max_step,
        double max_affine);
};

struct GeodesicIntegrationResult {
    PhaseSpaceState final_state;
    IntegrationDiagnostics diagnostics;
    std::optional<EventHit> event;
};

class GeodesicIntegrator {
public:
    explicit GeodesicIntegrator(const Metric& metric);

    GeodesicIntegrationResult integrate(
        const PhaseSpaceState& initial,
        const GeodesicIntegrationConfig& config,
        const std::vector<GeodesicEvent>& events = {}) const;
};
```

`cpu_reference` uses the v3 defaults:

- `rtol=1e-11`;
- `atol_t=1e-11*M`;
- `atol_x=1e-12*M`;
- `atol_p=1e-12`;
- `min_step=1e-12*M`;
- `max_rejections_per_step=12`;
- `max_total_steps=2,000,000`;
- controller `safety=0.9`, `min_factor=0.2`,
  `max_factor=5.0`;
- ordinary constraint tolerance `1e-10`.

`mass_scale`, `initial_step`, `max_step`, and `max_affine` must be finite and
positive in magnitude where applicable. `initial_step` carries the integration
direction. `max_proper_time` and `max_coordinate_time` default to infinity and
become active only when the caller supplies finite positive limits.

All three limits are displacements from the initial state:

- affine limit: `abs(lambda-lambda_initial)`;
- proper-time limit for a normalized unit-mass timelike state:
  `abs(lambda-lambda_initial)`;
- coordinate-time limit: `abs(x^0-x^0_initial)`.

The affine step is capped to land exactly on its configured limit.
Proper-time and coordinate-time crossings inside an accepted step are located
through the same dense-output bracketed root machinery, so their returned
states do not silently overshoot the requested limit.

The flow:

1. validates config, initial finite state, metric point, and initial
   constraint;
2. records initial E/Lz only when their monitoring flags are explicitly
   enabled;
3. caps each attempted step by `max_step` and the remaining affine limit;
4. evaluates one DOPRI5 trial;
5. treats an invalid trial metric point as a rejected step and shrinks, so a
   large trial cannot falsely terminate a still-valid trajectory;
6. terminates transparently after repeated invalid-domain trials, non-finite
   stages, rejection exhaustion, or step underflow;
7. after acceptance, checks all events through dense output and terminates at
   the first hit;
8. updates finite-state, constraint, optional invariant, proper-time,
   coordinate-time, and step diagnostics;
9. returns an explicit `TerminationReason` and message on every path.

For a normalized unit-mass timelike state, affine parameter may represent
proper time. A finite `max_proper_time` is therefore implemented as absolute
affine displacement for `TimelikeUnitMass` only and rejected for null
configuration.

The integrator does not retain every accepted state. This prevents a default
two-million-step reference run from allocating an unbounded trajectory.
Future sampling, radiative transfer, and rendering consume accepted steps
through later phase-specific flows.

## Failure semantics

- Invalid initial metric point: `InvalidMetricPoint`, zero accepted steps.
- Invalid trial point: reject and shrink; repeated failure becomes
  `InvalidMetricPoint`, never `HorizonCrossing` in BL.
- Non-finite initial state, stage derivative, dense state, or accepted state:
  `NonFiniteState`.
- Constraint above configured gate initially or after acceptance:
  `ConstraintViolation`.
- Step magnitude below `min_step`: `StepUnderflow`.
- More than `max_rejections_per_step` consecutive rejections:
  `StepUnderflow` for numerical error or `InvalidMetricPoint` for repeated
  metric-domain rejection, with a diagnostic message distinguishing them.
- Event evaluation/root iteration failure: `EventRootFailure`.
- Configured affine/proper/coordinate/step limits use their matching
  termination reasons.
- A user event returns its declared reason and event payload.
- No exception is swallowed. Invalid API configuration throws
  `std::invalid_argument`; trajectory-domain outcomes return diagnostics.

## Invariant monitoring

Energy `E=-p_t` and axial angular momentum `Lz=p_phi` are not universal
invariants for arbitrary metrics. The integrator therefore updates their drift
only when explicitly requested. For a nonzero initial invariant, drift is
`abs(current-initial)/abs(initial)`; for an exactly zero initial invariant, the
field records absolute drift. Unmonitored fields are NaN.

`max_carter_rel_error` remains NaN in Phase 1 because the general metric
interface does not define a Carter constant. Kerr Carter validation belongs to
Phase 2/3.

## Validation design

### Shared DOPRI5

- Test `y'=y` against `exp(t)` and observe fifth-order convergence.
- Force one rejection with tight tolerance and verify sign-preserving shrink.
- Verify component-specific absolute tolerances and the v3 RMS norm with
  hand-computed stage-error scaling.
- Integrate with a negative step and preserve direction.
- Verify dense output at both endpoints and internal points against an
  analytic polynomial/exponential solution.
- Reject zero/non-finite step, invalid tolerances, invalid controller factors,
  and non-finite RHS output.
- Record a golden legacy `dopri5_generic_step` result before refactoring and
  require the compatibility wrapper to retain acceptance, error, state, and
  next-step behavior.

### Hamiltonian

- Minkowski null state: exact `H=0`, `dp=0`, and constant tangent.
- Minkowski timelike state: exact `H=-1/2`.
- Hand-compute the normalized constraint denominator, including momentum
  scaling.
- Compare general RHS against a direct independent contraction at
  Schwarzschild and Kerr points.
- Require exact `dp_t=dp_phi=0` for current stationary/axisymmetric metrics.
- Use a test-only time-dependent metric and require nonzero `dp_t`, proving no
  universal conservation assumption was inserted.
- Reject invalid metric points and non-finite canonical momenta.

### Events

- Increasing, decreasing, and any-direction roots.
- Exact start and end roots.
- No root for a direction mismatch.
- Negative affine step with direction defined along step progression.
- First of multiple roots inside one accepted step.
- Non-finite event function and exhausted root iterations.
- Root affine error below `1e-10*M` for external geometric events.

### Geodesics

- Minkowski null line agrees with its analytic line to floating-point scale.
- Minkowski timelike inertial line agrees with its analytic line and retains
  `H=-1/2`.
- Forward integration followed by equal backward integration recovers the
  initial canonical state.
- Schwarzschild outgoing radial null motion has linear `r(lambda)` and the
  analytic tortoise-coordinate time relation.
- Schwarzschild weak-field null deflection at large impact parameter converges
  toward `4M/b`.
- At `r=3M`, `b=3*sqrt(3)*M`, the equatorial photon-sphere canonical RHS has
  zero radial velocity and radial momentum derivative.
- Current Kerr BL metrics conserve monitored E and Lz to the v3 ordinary-ray
  gate where the general RHS provides exact zero derivatives.
- Ordinary-ray maximum normalized Hamiltonian error is below `1e-10`.
- Initial constraint violation, max steps, max affine, max proper time, max
  coordinate time, step underflow, invalid BL domain, non-finite state, and
  user event return their exact reasons.
- A deliberately perturbed DOPRI coefficient or Hamiltonian derivative must
  make an analytic validation fail before the mutation is restored.

### Regression and runtime checks

- All existing Phase 0A/0B and legacy tests remain green.
- AddressSanitizer and UndefinedBehaviorSanitizer run all relativity tests.
- Release build uses the existing C++17 warning flags.
- `git diff --check` is clean.

## Compatibility and migration

- No existing public Solar integrator signature is removed or renamed.
- Existing CR3BP callers continue using `dopri5_generic_step`.
- The specialized N-body adaptive solver is not behaviorally rewritten.
- `Metric` remains unchanged.
- `PhaseSpaceState` and `GeodesicSample` retain their Phase 0A layouts.
- New types live under `solar::relativity` or `solar::numerics`, preventing
  collision with the legacy `solar::State` and adaptive result types.
- No dependency beyond the C++17 standard library is added.

## Explicitly not completed

- Tetrad/local observer initialization and observer-frequency normalization.
- General future-direction classification relative to an observer.
- Automatic constraint projection.
- Carter constant and Kerr separated/Mino-time integration.
- Kerr–Schild coordinates or physical horizon crossing.
- Disk/material models, radiative transfer, ray CLI, rendering, WASM/GPU, and
  UI.
- Claiming BL invalid-domain termination as physical horizon capture.

## Gate

Phase 1 may be marked `PASSED` only when:

- Hamilton RHS, adapted DOPRI5, dense output, directed events, diagnostics,
  Minkowski/Schwarzschild geodesics, and null/timelike constraints are present;
- every required test and full regression exits successfully;
- ordinary normalized Hamiltonian error is below `1e-10`;
- failures and unavailable invariants are visible rather than silently
  classified;
- the validation report records platform, commands, thresholds, model
  boundary, likely bugs, and fastest falsification paths;
- no Phase 2 observer/tetrad implementation has entered the branch.
