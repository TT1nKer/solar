# Solar Relativity Phase 2 Observer and Kerr Benchmark Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build and independently validate the observer, tetrad, local
initialization, Kerr invariant, analytic orbit, and shadow-benchmark
foundation required by v3 Phase 2.

**Architecture:** Focused L1 value/function modules own spacetime algebra,
observer frames, canonical initialization, Kerr constants, orbit formulas,
and the Bardeen curve. The existing L2 geodesic flow receives an optional
Carter evaluator and remains independent of the concrete Kerr metric.
Numerical shadow comparison is validation-only and never treats a
Boyer-Lindquist chart failure as physical horizon capture.

**Tech Stack:** C++17, GNU Make, `std::array`, `std::optional`,
`std::function`, existing `solar::relativity::Metric`, Kerr BL metric,
Hamiltonian geodesic integrator, and standalone executable tests.

## Global Constraints

- Signature is `(-,+,+,+)`; coordinate zero is time.
- Boyer-Lindquist order is `(t,r,theta,phi)`.
- Coordinates/four-velocities are contravariant and canonical momenta are
  covariant.
- Units are geometrized with `G=c=1`.
- Authoritative photon momentum is future-directed; backward tracing uses a
  negative affine step.
- Observer frequency is `nu=-p_mu u_obs^mu`.
- Locally initialized photon momentum is normalized to `nu_obs=1`.
- Tetrad and Hamiltonian normalized errors must remain below `1e-10`.
- Carter relative error uses `abs(Q-Q0)/max(1,abs(Q0))`; absolute error is also
  reported.
- Static observers do not exist where `g_tt>=0`.
- Circular timelike existence and stability are distinct outcomes.
- `OrbitSense` is relative to black-hole spin; zero spin uses positive `phi`
  as the prograde convention.
- The asymptotic Bardeen formula is not a finite-distance hard oracle.
- Boyer-Lindquist inner-boundary events are capture proxies only and are never
  `HorizonCrossing`.
- No projection, external dependency, renderer, matter, transfer, separated
  solver, Kerr-Schild chart, GPU/WASM, or UI is added.
- Existing Phase 1 public call signatures remain source-compatible.

---

### Task 1: Spacetime algebra primitives

**Files:**
- Create: `include/solar/relativity/spacetime_algebra.h`
- Create: `tests/relativity/test_spacetime_algebra.cpp`

**Interfaces:**
- Consumes: `Mat4`, `Contravariant4`, and `Covariant4`.
- Produces: `metric_inner_product`, `covector_vector_pairing`, `lower_index`,
  and `raise_index`.

- [x] **Step 1: Write the failing real-behavior test**

Create `tests/relativity/test_spacetime_algebra.cpp` with a local pass/fail
counter. Use the literal Minkowski matrix and vectors:

```cpp
const Mat4 eta{{
    {{-1.0, 0.0, 0.0, 0.0}},
    {{ 0.0, 1.0, 0.0, 0.0}},
    {{ 0.0, 0.0, 1.0, 0.0}},
    {{ 0.0, 0.0, 0.0, 1.0}},
}};
const Contravariant4 u{Vec4{{2.0, 0.5, -1.0, 3.0}}};
const Contravariant4 v{Vec4{{1.5, -2.0, 4.0, 0.25}}};

check_near("literal Minkowski inner product",
           metric_inner_product(eta, u, v), -7.25, 0.0);

const Covariant4 lowered = lower_index(eta, u);
check_near("lowered time component", lowered.v[0], -2.0, 0.0);
check_near("lowered spatial component", lowered.v[2], -1.0, 0.0);

const Contravariant4 raised = raise_index(eta, lowered);
check_near("raise/lower time round trip", raised.v[0], 2.0, 0.0);
check_near("covector-vector pairing",
           covector_vector_pairing(lowered, v), -7.25, 0.0);
```

Add a non-diagonal finite matrix case with hand-derived literal output and a
non-finite input case that propagates NaN rather than silently fabricating a
finite contraction.

Mutation target: swapping covariant indices, omitting the time sign, or pairing
two covectors must fail a literal assertion.

- [x] **Step 2: Run and verify RED**

```bash
make tests/relativity/test_spacetime_algebra
```

Expected: compilation fails because
`solar/relativity/spacetime_algebra.h` does not exist.

- [x] **Step 3: Implement the minimal algebra**

Implement direct four-component sums. `lower_index` uses
`p_mu=g_mu,nu v^nu`; `raise_index` uses
`v^mu=g^mu,nu p_nu`. Do not add a tensor class or cache metric evaluations.

- [x] **Step 4: Verify GREEN**

```bash
make tests/relativity/test_spacetime_algebra
./tests/relativity/test_spacetime_algebra
```

Expected: all assertions pass with no compiler warning.

- [x] **Step 5: Commit**

```bash
git add include/solar/relativity/spacetime_algebra.h \
  tests/relativity/test_spacetime_algebra.cpp
git commit -m "feat: add spacetime algebra primitives"
```

### Task 2: Generic, static, and look-at observer tetrads

**Files:**
- Create: `include/solar/relativity/observer.h`
- Create: `src/relativity/observer.cpp`
- Create: `tests/relativity/test_observers.cpp`

**Interfaces:**
- Consumes: Task 1 algebra and any `Metric`.
- Produces: `Tetrad`, `ObserverFrame`, `ObserverError`, `ObserverResult`,
  `LookAtAttitude`, tetrad transforms, arbitrary/static/look-at constructors.

- [x] **Step 1: Add failing observer-contract tests**

Create `tests/relativity/test_observers.cpp`. For Minkowski:

```cpp
MinkowskiMetric metric;
const Contravariant4 x{Vec4{{0.0, 0.0, 0.0, 0.0}}};

const ObserverResult stationary = make_static_observer(metric, x);
check("Minkowski static observer exists", bool(stationary));
check_near("static tetrad orthonormal",
           tetrad_orthonormality_error(metric, *stationary.frame),
           0.0, 1.0e-15);

const double gamma = 1.25; // v=0.6
const Contravariant4 boosted_u{Vec4{{gamma, 0.75, 0.0, 0.0}}};
const LookAtAttitude attitude{
    Contravariant4{Vec4{{0.0, 0.0, 0.0, -2.0}}},
    Contravariant4{Vec4{{0.0, 0.0, 3.0, 0.0}}},
};
const ObserverResult look_at =
    make_look_at_observer(metric, x, boosted_u, attitude);
check("boosted look-at observer exists", bool(look_at));
check_near("boosted tetrad orthonormal",
           tetrad_orthonormality_error(metric, *look_at.frame),
           0.0, 2.0e-15);
```

Project coordinate vectors back to local components and assert a complete
round trip using literals. Require the projected look leg to have positive
inner product with the supplied look direction and the up leg with the
supplied up reference. Assert the spatial-component determinant is positive.

Add failures for:

```cpp
const Contravariant4 non_unit_u{Vec4{{1.0, 0.5, 0.0, 0.0}}};
check("non-unit four-velocity rejected",
      make_look_at_observer(
          metric, x, non_unit_u, attitude).error ==
          ObserverError::FourVelocityNotUnitTimelike);

const LookAtAttitude parallel{
    attitude.look_direction,
    attitude.look_direction,
};
check("parallel look/up rejected",
      make_look_at_observer(
          metric, x, boosted_u, parallel).error ==
          ObserverError::DegenerateSpatialSeed);
```

For Kerr `M=1`, `chi=0.8`, `r=1.8`, `theta=pi/2`, require the metric point to
be valid and `make_static_observer` to return exactly
`ObserverError::StaticWorldlineNotTimelike`.

Mutation target: the wrong Lorentzian projection sign, Euclidean
Gram-Schmidt, missing ergosphere test, or omitted handedness flip must fail.

- [x] **Step 2: Run and verify RED**

```bash
make tests/relativity/test_observers
```

Expected: compilation fails because `observer.h` is absent.

- [x] **Step 3: Implement observer contracts and Lorentzian Gram-Schmidt**

Implement the exact types from the design. Use:

```cpp
candidate = seed;
candidate = candidate +
    metric_inner_product(g, candidate, e0) * e0.v;
for (const spatial_leg : completed_spatial_legs) {
    candidate = candidate -
        metric_inner_product(g, candidate, spatial_leg) *
            spatial_leg.v;
}
```

Normalize only for a finite positive norm. Build the right-handed spatial
determinant explicitly. Re-evaluate all 16 tetrad products and reject error
above `1e-10`.

Static observer uses `u^t=1/sqrt(-g_tt)` and coordinate `r`, `theta`, `phi`
seeds. Check `g_tt<0` before the square root.

- [x] **Step 4: Verify GREEN and focused regressions**

```bash
make tests/relativity/test_observers tests/relativity/test_metrics
./tests/relativity/test_observers
./tests/relativity/test_metrics
```

Expected: observer and existing metric tests pass.

- [x] **Step 5: Commit**

```bash
git add include/solar/relativity/observer.h \
  src/relativity/observer.cpp \
  tests/relativity/test_observers.cpp
git commit -m "feat: add observer tetrad construction"
```

### Task 3: Kerr ZAMO observer

**Files:**
- Modify: `include/solar/relativity/observer.h`
- Modify: `src/relativity/observer.cpp`
- Modify: `tests/relativity/test_observers.cpp`

**Interfaces:**
- Consumes: `KerrBoyerLindquistMetric` and Task 2 observer contracts.
- Produces: `make_zamo_observer`.

- [x] **Step 1: Add failing ZAMO tests**

For `M=1`, `chi=0.7`, `x=(0,8,1.1,0.3)`:

```cpp
const ObserverResult zamo = make_zamo_observer(kerr, x);
check("ordinary Kerr ZAMO exists", bool(zamo));
check("ZAMO tetrad error below gate",
      tetrad_orthonormality_error(kerr, *zamo.frame) < 1.0e-12);

const Mat4 g = kerr.covariant(x);
const Covariant4 zamo_momentum =
    lower_index(g, zamo.frame->tetrad.basis[0]);
check_near("ZAMO axial angular momentum",
           zamo_momentum.v[3], 0.0, 2.0e-15);
check("ZAMO future coordinate time",
      zamo.frame->tetrad.basis[0].v[0] > 0.0);
```

At `r=1.0e6`, compare the ZAMO time leg to `(1,0,0,0)` and the spatial legs to
the asymptotic normalized spherical-coordinate legs. Require convergence
within `3e-6`. Require invalid BL points and the polar axis to return
`ObserverError::InvalidMetricPoint`.

Mutation target: omitting `omega/alpha`, using `omega` with the wrong sign, or
using `sqrt(Sigma/A)*sin(theta)` instead of division must fail.

- [x] **Step 2: Run and verify RED**

```bash
make tests/relativity/test_observers
```

Expected: compilation fails because `make_zamo_observer` is undeclared.

- [x] **Step 3: Implement the explicit v3 ZAMO basis**

Compute `Sigma`, `Delta`, `A`, `alpha`, `omega`, and the four legs from v3.
Reject every non-finite/non-positive square-root input. Reuse final tetrad
validation from Task 2; do not route the explicit ZAMO through generic
Gram-Schmidt.

- [x] **Step 4: Verify GREEN**

```bash
make tests/relativity/test_observers
./tests/relativity/test_observers
```

Expected: all observer tests pass.

- [x] **Step 5: Commit**

```bash
git add include/solar/relativity/observer.h \
  src/relativity/observer.cpp \
  tests/relativity/test_observers.cpp
git commit -m "feat: add Kerr ZAMO observer"
```

### Task 4: Local photon and timelike initialization

**Files:**
- Create: `include/solar/relativity/local_initialization.h`
- Create: `src/relativity/local_initialization.cpp`
- Create: `tests/relativity/test_local_initialization.cpp`

**Interfaces:**
- Consumes: observer frames, spacetime algebra, Hamiltonian constraint.
- Produces: `InitialStateError`, `InitialStateResult`,
  `observer_measured_frequency`, `initialize_local_photon`, and
  `initialize_local_timelike`.

- [x] **Step 1: Add failing local-initialization tests**

Use a Minkowski static observer. For photon direction `(2,0,0)`:

```cpp
const InitialStateResult photon =
    initialize_local_photon(
        metric, *observer.frame, Vec3{{2.0, 0.0, 0.0}});
check("local photon initializes", bool(photon));
check_near("photon p_t", photon.state->p.v[0], -1.0, 1.0e-15);
check_near("photon p_x", photon.state->p.v[1], 1.0, 1.0e-15);
check_near("photon observer frequency",
           photon.measured_frequency, 1.0, 1.0e-15);
check("photon null constraint",
      hamiltonian_constraint_error(
          metric, *photon.state, GeodesicKind::Null) < 1.0e-14);
```

Raise and project the photon momentum and require local components
`(1,1,0,0)` within `2e-15`.

For local velocity `(0.6,0,0)`, require:

```text
p_t = -1.25
p_x =  0.75
H   = -0.5
measured local energy = 1.25
```

Require zero photon direction, non-finite direction, `|v|=1`, `|v|>1`, a
non-orthonormal observer, and a non-finite affine parameter to return the
specific error enum without a state.

Add a real backward-integration test: initialize a future outward Minkowski
photon at `x=10`, integrate with negative `initial_step`, require final
`x<10`, and require its measured frequency against the original observer
velocity to stay positive.

Mutation target: forgetting to normalize direction/frequency, lowering with
the inverse metric, using `+p.u`, reversing photon momentum for backward
tracing, or omitting gamma must fail.

- [x] **Step 2: Run and verify RED**

```bash
make tests/relativity/test_local_initialization
```

Expected: compilation fails because `local_initialization.h` is absent.

- [x] **Step 3: Implement minimal initialization**

Form local photon `(1,n/|n|)` and timelike
`gamma(1,v)`. Expand through the tetrad, lower with `g_mu,nu`, and evaluate
frequency with the covector-vector pairing. Scale only photon momentum by
`1/nu`. Validate the observer tetrad and Hamiltonian constraint before
returning success.

- [x] **Step 4: Verify GREEN**

```bash
make tests/relativity/test_local_initialization \
  tests/relativity/test_geodesics
./tests/relativity/test_local_initialization
./tests/relativity/test_geodesics
```

Expected: local and existing geodesic assertions pass.

- [x] **Step 5: Commit**

```bash
git add include/solar/relativity/local_initialization.h \
  src/relativity/local_initialization.cpp \
  tests/relativity/test_local_initialization.cpp
git commit -m "feat: initialize local physical states"
```

### Task 5: Kerr constants and opt-in Carter diagnostics

**Files:**
- Create: `include/solar/relativity/kerr_constants.h`
- Create: `src/relativity/kerr_constants.cpp`
- Create: `tests/relativity/test_kerr_constants.cpp`
- Modify: `include/solar/relativity/geodesic_types.h`
- Modify: `include/solar/relativity/geodesic_integrator.h`
- Modify: `src/relativity/geodesic_config.cpp`
- Modify: `src/relativity/geodesic_integrator.cpp`
- Modify: `tests/relativity/test_geodesic_events.cpp`
- Modify: `tests/relativity/test_geodesics_kerr.cpp`

**Interfaces:**
- Consumes: canonical state, Kerr metric, geodesic kind.
- Produces: `KerrConstants`, `evaluate_kerr_constants`, optional
  `carter_evaluator`, relative and absolute Carter diagnostics.

- [x] **Step 1: Add failing literal Kerr-constant tests**

For `M=1`, `chi=0.5`, `theta=pi/3`,
`p_t=-1`, `p_theta=3`, `p_phi=2`, assert:

```cpp
const KerrConstants null_constants =
    evaluate_kerr_constants(kerr, state, GeodesicKind::Null);
check_near("Kerr E from covariant p_t",
           null_constants.E, 1.0, 0.0);
check_near("Kerr Lz from covariant p_phi",
           null_constants.Lz, 2.0, 0.0);
check_near("literal null Carter Q",
           null_constants.Q,
           10.270833333333333, 2.0e-14);
check_near("null mass squared",
           null_constants.mass_sq, 0.0, 0.0);

const KerrConstants timelike_constants =
    evaluate_kerr_constants(
        kerr, state, GeodesicKind::TimelikeUnitMass);
check_near("literal timelike Carter Q",
           timelike_constants.Q,
           10.333333333333333, 2.0e-14);
```

Require equatorial `p_theta=0` to give exact `Q=0`. Reject unknown kind,
non-finite state, invalid metric point, and the polar axis.

Mutation target: using contravariant momentum, omitting `cos^2`, wrong mass
term sign, or using `Lz^2*sin^2` must fail.

- [x] **Step 2: Run and verify RED for the constants API**

```bash
make tests/relativity/test_kerr_constants
```

Expected: compilation fails because `kerr_constants.h` is absent.

- [x] **Step 3: Implement constants and verify focused GREEN**

Implement the exact v3 equation with finite/domain validation:

```bash
make tests/relativity/test_kerr_constants
./tests/relativity/test_kerr_constants
```

Expected: literal constant tests pass.

- [x] **Step 4: Add failing geodesic diagnostic tests**

Append `max_carter_abs_error` to `IntegrationDiagnostics`; existing default
diagnostic tests must require both Carter fields to be NaN when no evaluator
is supplied.

In `test_geodesics_kerr.cpp`, initialize a generic non-equatorial null ray
from a ZAMO with local direction `(0.3,0.4,0.5)`. Configure:

```cpp
config.monitor_energy = true;
config.monitor_lz = true;
config.carter_evaluator =
    [&metric](const PhaseSpaceState& state) {
        return evaluate_kerr_constants(
            metric, state, GeodesicKind::Null).Q;
    };
```

Require an ordinary termination, maximum normalized Hamiltonian error below
`1e-10`, Carter relative and absolute error below `1e-10`, and finite Carter
diagnostics.

Add a Minkowski monitor fixture:

```cpp
config.carter_evaluator =
    [](const PhaseSpaceState& state) {
        return 1.0e-6 + 1.0e-9 * state.affine;
    };
```

For affine displacement `2`, require both absolute and v3 normalized error to
be `2e-9` within integration roundoff. This catches accidental division by
the small initial value. Require a NaN-returning or throwing evaluator to
terminate explicitly without accepting a fabricated diagnostic.

- [x] **Step 5: Run and verify RED for monitoring**

```bash
make tests/relativity/test_geodesic_events \
  tests/relativity/test_geodesics_kerr
```

Expected: compilation fails because the absolute field/evaluator is absent.

- [x] **Step 6: Implement the generic evaluator path**

Add:

```cpp
using InvariantEvaluator =
    std::function<double(const PhaseSpaceState&)>;
```

and an empty `carter_evaluator` to config. Evaluate its initial value after
state/metric/constraint validation but before event or step work. Update on
accepted states. Use:

```cpp
const double difference = std::fabs(current - initial);
const double normalized =
    difference / std::max(1.0, std::fabs(initial));
```

Apply the same denominator to E and Lz. Catch evaluator exceptions and reject
non-finite output with an explicit diagnostic message.

- [x] **Step 7: Verify GREEN and Phase 1 compatibility**

```bash
make tests/relativity/test_kerr_constants \
  tests/relativity/test_geodesic_events \
  tests/relativity/test_geodesics_kerr \
  tests/relativity/test_geodesics
./tests/relativity/test_kerr_constants
./tests/relativity/test_geodesic_events
./tests/relativity/test_geodesics_kerr
./tests/relativity/test_geodesics
```

Expected: new Carter behavior and all existing geodesic behavior pass.

- [x] **Step 8: Commit**

```bash
git add include/solar/relativity/kerr_constants.h \
  src/relativity/kerr_constants.cpp \
  tests/relativity/test_kerr_constants.cpp \
  include/solar/relativity/geodesic_types.h \
  include/solar/relativity/geodesic_integrator.h \
  src/relativity/geodesic_config.cpp \
  src/relativity/geodesic_integrator.cpp \
  tests/relativity/test_geodesic_events.cpp \
  tests/relativity/test_geodesics_kerr.cpp
git commit -m "feat: monitor Kerr Carter invariant"
```

### Task 6: Kerr special and circular timelike orbits

**Files:**
- Create: `include/solar/relativity/kerr_orbits.h`
- Create: `src/relativity/kerr_orbits.cpp`
- Create: `tests/relativity/test_kerr_orbits.cpp`

**Interfaces:**
- Consumes: Kerr metric, observer construction, spacetime algebra.
- Produces: `OrbitSense`, `CircularOrbitStability`,
  `CircularTimelikeOrbit`, `CircularOrbitResult`, analytic radii/properties,
  and `make_equatorial_circular_observer`.

- [x] **Step 1: Add failing special-radius tests**

For Schwarzschild limit `M=2`, require:

```text
ISCO                = 12
photon orbit        = 6
marginally bound    = 8
```

for both senses. For `M=1`, `chi=0.5`, require:

```text
prograde ISCO       = 4.233002529530826
retrograde ISCO     = 7.554584714512358
prograde photon     = 2.347296355333861
retrograde photon   = 3.532088886237956
prograde mb         = 2.914213562373095
retrograde mb       = 4.949489742783178
```

Repeat with `chi=-0.5` and require radii to retain prograde/retrograde
meaning relative to spin. Reject unknown sense.

Mutation target: signed-spin radii, swapped prograde/retrograde square-root
branch, or a missing mass scale must fail literal fixtures.

- [x] **Step 2: Run and verify RED**

```bash
make tests/relativity/test_kerr_orbits
```

Expected: compilation fails because `kerr_orbits.h` is absent.

- [x] **Step 3: Implement special radii**

Use `abs(chi)` in relative-spin radius formulas and multiply all
dimensionless radii by `M`. Use `std::cbrt` for the ISCO cube roots. Validate
finite outputs and sense.

- [x] **Step 4: Add failing circular-orbit and observer tests**

At Schwarzschild `M=1`, `r=6`:

```text
specific E          = 0.9428090415820634
prograde Lz         = 3.4641016151377544
retrograde Lz       = -3.4641016151377544
abs(Omega)          = 0.06804138174397717
stability           = Stable
```

At `r=5`, require a successful `Unstable` timelike orbit. At `r=3`, require
failure rather than incorrectly calling it unstable timelike.

For `chi=-0.5`, require prograde `Omega<0` and `Lz<0`; retrograde signs must
be positive. Build circular observers at `theta=pi/2`, lower their time legs,
and require:

```cpp
check_near("circular E matches lowered observer",
           -lowered.v[0], orbit.specific_energy, 2.0e-13);
check_near("circular Lz matches lowered observer",
           lowered.v[3], orbit.specific_lz, 2.0e-13);
check("circular tetrad gate",
      tetrad_orthonormality_error(metric, *observer.frame) < 1.0e-10);
```

Reject non-equatorial input, invalid metric point, and a radius without a
timelike circular orbit using
`ObserverError::CircularWorldlineNotTimelike`.

Mutation target: conflating `r<ISCO` with nonexistence, the wrong coordinate
rotation sign, or mismatched analytic/lowered E/Lz must fail.

- [x] **Step 5: Implement circular quantities and observer**

Use dimensionless `x=r/M` and coordinate rotation sign `s`:

```text
M Omega = s / (x^(3/2) + s chi)
E = (x^(3/2)-2sqrt(x)+s chi) /
    (x^(3/4) sqrt(x^(3/2)-3sqrt(x)+2s chi))
Lz/M = s (x^2-2s chi sqrt(x)+chi^2) /
       (x^(3/4) sqrt(x^(3/2)-3sqrt(x)+2s chi))
```

Require positive finite normalization. Construct
`u=(u^t,0,0,Omega*u^t)` from the metric and route its spatial basis through
the generic observer constructor.

- [x] **Step 6: Verify GREEN**

```bash
make tests/relativity/test_kerr_orbits \
  tests/relativity/test_observers
./tests/relativity/test_kerr_orbits
./tests/relativity/test_observers
```

Expected: all orbit and observer tests pass.

- [x] **Step 7: Commit**

```bash
git add include/solar/relativity/kerr_orbits.h \
  src/relativity/kerr_orbits.cpp \
  tests/relativity/test_kerr_orbits.cpp
git commit -m "feat: add Kerr circular orbit observers"
```

### Task 7: Analytic Bardeen shadow curve

**Files:**
- Create: `include/solar/relativity/kerr_shadow.h`
- Create: `src/relativity/kerr_shadow.cpp`
- Create: `tests/relativity/test_kerr_shadow.cpp`

**Interfaces:**
- Consumes: Kerr special photon radii.
- Produces: `ShadowCriticalPoint` and `bardeen_shadow_curve`.

- [x] **Step 1: Add failing analytic-curve tests**

For `M=2`, `chi=0`, require every point to satisfy:

```cpp
const double radius_squared =
    point.alpha * point.alpha + point.beta * point.beta;
check_near("Schwarzschild shadow circle",
           radius_squared, 108.0, 5.0e-13);
check_near("Schwarzschild photon radius",
           point.photon_radius, 6.0, 0.0);
```

For `M=1`, `chi=0.5`, equatorial inclination, require points at the two
zero-beta edges:

```text
left alpha          = -4.096266658713869
right alpha         =  6.138155724715452
```

Require every non-axis point to have a reflected partner `(alpha,-beta)`,
all components finite, and all sampled photon radii inside the analytic
physical interval. With `chi=-0.5`, require the curve to mirror horizontally.
With `M=3`, require alpha, beta, and photon radius to scale by 3.

Require exact zero and `|chi|<=64*sqrt(epsilon)` to use the Schwarzschild
limit. Reject inclination `0`, `pi`, NaN, and `samples_per_branch<2`.

Mutation target: wrong xi denominator sign, omitted cotangent term, accepting
a truly negative beta radicand, or direct division by tiny spin must fail.

- [x] **Step 2: Run and verify RED**

```bash
make tests/relativity/test_kerr_shadow
```

Expected: compilation fails because `kerr_shadow.h` is absent.

- [x] **Step 3: Implement the analytic curve**

For rotating Kerr, sample the closed interval between the relative-spin
prograde and retrograde photon radii. Compute `xi`, `eta`, and the beta
radicand in dimensionless long-double intermediates, then restore mass scale.
Clamp a negative radicand only when:

```cpp
radicand >= -128.0L *
    std::numeric_limits<long double>::epsilon() *
    std::max(1.0L, radicand_term_scale)
```

Return upper branch in increasing photon radius and lower branch in reverse,
without duplicating zero-beta endpoints.

- [x] **Step 4: Verify GREEN**

```bash
make tests/relativity/test_kerr_shadow
./tests/relativity/test_kerr_shadow
```

Expected: all analytic curve assertions pass.

- [x] **Step 5: Commit**

```bash
git add include/solar/relativity/kerr_shadow.h \
  src/relativity/kerr_shadow.cpp \
  tests/relativity/test_kerr_shadow.cpp
git commit -m "feat: add analytic Kerr shadow curve"
```

### Task 8: CPU backward-ray shadow cross-check and Phase 2 gate

**Files:**
- Modify: `tests/relativity/test_kerr_shadow.cpp`
- Create: `docs/validation/relativity_03_tetrad.md`
- Create: `docs/validation/relativity_05_carter.md`
- Create: `docs/validation/relativity_06_shadow.md`
- Create: `docs/validation/relativity_07_timelike.md`
- Modify: `RELATIVITY_STATUS.md`

**Interfaces:**
- Consumes: all Phase 2 modules and the Phase 1 Hamiltonian integrator.
- Produces: independent numerical/analytic shadow evidence and Phase 2
  validation reports.

- [x] **Step 1: Add the numerical shadow acceptance benchmark**

In `test_kerr_shadow.cpp`, use `M=1`, `chi=0.5`, equatorial observer radius
`r_obs=1000`. Build a ZAMO. For horizontal Bardeen screen coordinate
`alpha`, initialize:

```cpp
const double local_phi = -alpha / r_obs;
const Vec3 local_direction{{
    std::sqrt(1.0 - local_phi * local_phi),
    0.0,
    local_phi,
}};
```

This is the future photon arriving radially outward at the camera. Integrate
with a negative affine step and two explicit events:

```cpp
const double inner_radius =
    metric.outer_horizon_radius() + 1.0e-3;
const double escape_radius = 1.1 * r_obs;

GeodesicEvent inner{
    "BL capture proxy",
    [inner_radius](const PhaseSpaceState& state) {
        return state.x.v[1] - inner_radius;
    },
    EventDirection::Decreasing,
    TerminationReason::UserEvent,
    1.0e-10,
};
GeodesicEvent escape{
    "outer escape",
    [escape_radius](const PhaseSpaceState& state) {
        return state.x.v[1] - escape_radius;
    },
    EventDirection::Increasing,
    TerminationReason::Escaped,
    1.0e-10,
};
```

Use `max_step=2`, `max_affine=4000`, Hamiltonian tolerance `1e-10`, and Carter
monitoring. Classify only inner event index 0 as captured and escape index 1
as escaped. Any metric failure, root failure, constraint failure, or affine
limit fails the benchmark.

Binary-search left bracket `[-8,0]` and right bracket `[0,8]` to screen
tolerance `1e-3`. Require:

```text
abs(numerical_left  - analytic_left)  < 3e-2
abs(numerical_right - analytic_right) < 3e-2
max Hamiltonian error                < 1e-10
max Carter relative error            < 1e-10
```

Also require every initialized ray to have positive observer frequency while
being integrated backward with negative affine steps.

Mutation target: past-directed camera momentum, wrong screen horizontal sign,
mapping `InvalidMetricPoint` to capture, or a wrong analytic boundary must
fail.

- [x] **Step 2: Run the independent acceptance benchmark**

```bash
make tests/relativity/test_kerr_shadow
./tests/relativity/test_kerr_shadow
```

This task adds no production API: Tasks 1-7 were each developed through RED
before implementation. The cross-module acceptance benchmark may therefore
pass on its first run. If it fails, record the actual mismatch and add a
narrow failing regression before changing production code. A timeout or
non-event outcome is a failure, not a skip.

- [x] **Step 3: Make only benchmark-supported corrections**

If the first run exposes a defect, reproduce it with the narrowest failing
assertion before changing production code. Do not loosen the `1e-10`
Hamiltonian/Carter gates. Finite-distance edge tolerance may be tightened
after a two-radius convergence sweep, but may not be enlarged beyond `3e-2`.

- [x] **Step 4: Run focused Phase 2 tests**

```bash
make \
  tests/relativity/test_spacetime_algebra \
  tests/relativity/test_observers \
  tests/relativity/test_local_initialization \
  tests/relativity/test_kerr_constants \
  tests/relativity/test_kerr_orbits \
  tests/relativity/test_kerr_shadow \
  tests/relativity/test_geodesics_kerr
./tests/relativity/test_spacetime_algebra
./tests/relativity/test_observers
./tests/relativity/test_local_initialization
./tests/relativity/test_kerr_constants
./tests/relativity/test_kerr_orbits
./tests/relativity/test_kerr_shadow
./tests/relativity/test_geodesics_kerr
```

Expected: every focused assertion passes.

- [x] **Step 5: Run complete release verification**

```bash
make clean
make
make test
git diff --check
```

Expected: all relativity and fixture-independent legacy assertions pass.
The optional DE440 fixture may remain a visible skip and is not counted.

- [x] **Step 6: Run sanitizers**

```bash
make clean
make CXXFLAGS='-std=c++17 -O1 -g -Wall -Wextra -Iinclude \
  -fsanitize=address,undefined -fno-omit-frame-pointer' \
  $(find tests/relativity -name 'test_*.cpp' -type f | sort | \
    sed 's/\.cpp$//') \
  tests/test_integrator
```

Run every built relativity executable and `tests/test_integrator` with:

```bash
ASAN_OPTIONS=detect_leaks=0:halt_on_error=1
UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1
```

Expected: no sanitizer diagnostic.

- [x] **Step 7: Write validation reports from actual output**

Each report must contain the v3 headings:

```text
Claim
Model boundary
Reference
Command
Inputs
Expected
Actual
Error
Result
Limitations
Fastest falsification
```

Record exact tetrad errors, photon/timelike constraints and frequencies,
Carter relative/absolute drift, special-orbit errors, analytic shadow
fixtures, numerical edge errors, compiler/platform, assertion totals, and
sanitizer results. Do not paste complete terminal logs.

- [x] **Step 8: Advance the status gate**

Only after Steps 4-7 pass, update:

```text
CURRENT_PHASE: 2
PHASE_STATE: PASSED
LAST_VERIFIED_COMMIT: exact output of `git rev-parse HEAD` immediately before
the validation-only documentation changes
NEXT_ALLOWED_ACTION: Phase 3 only
```

`RELATIVITY_STATUS.md` must preserve explicit missing work and list at least
three likely bugs plus fastest falsification commands.

- [x] **Step 9: Restore release artifacts and reverify documentation tree**

```bash
make clean
make
make test
git diff --check
git status --short
```

Expected: release build/test success and only intended documentation/status
changes before the evidence commit.

- [x] **Step 10: Commit the gate evidence**

```bash
git add tests/relativity/test_kerr_shadow.cpp \
  docs/validation/relativity_03_tetrad.md \
  docs/validation/relativity_05_carter.md \
  docs/validation/relativity_06_shadow.md \
  docs/validation/relativity_07_timelike.md \
  RELATIVITY_STATUS.md
git commit -m "test: pass relativity phase 2 gate"
```

## Plan self-review

- Every v3 Phase 2 deliverable maps to a task: tetrad/observer Tasks 2-3,
  local initialization Task 4, Kerr constants Task 5, analytic orbits Task 6,
  analytic and CPU shadow Tasks 7-8.
- Every production module is introduced only after a failing test.
- Public type and function names are consistent with the approved design.
- No task enters separated/Mino, Kerr-Schild, matter, transfer, renderer,
  frontend, or GPU scope.
- Physical nonexistence is explicit; no task fabricates fallback observer or
  capture behavior.
- No placeholder, hidden second solver, or external dependency is required.
