# Solar Relativity Phase 1 Hamiltonian Geodesic Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build and independently validate the general eight-dimensional
Hamiltonian geodesic foundation required by v3 Phase 1.

**Architecture:** A shared container-generic Dormand–Prince stage engine serves
the new fixed-size numerical API and the existing vector compatibility wrapper.
Focused L1 modules own Hamiltonian physics and bracketed event roots; one L2
flow owns adaptive integration and truthful termination diagnostics. Tetrad
initialization remains in Phase 2.

**Tech Stack:** C++17, GNU Make, `std::array`, existing
`solar::relativity::Metric` and tensor wrappers, standalone executable tests.

## Global Constraints

- Signature is `(-,+,+,+)`; coordinate zero is time.
- Coordinates are contravariant and canonical momenta are covariant.
- Units are geometrized with `G=c=1`.
- State order is `(x0,x1,x2,x3,p0,p1,p2,p3)`.
- Null target Hamiltonian is `0`; unit-mass timelike target is `-1/2`.
- CPU defaults are `rtol=1e-11`, `atol_t=1e-11 M`,
  `atol_x=1e-12 M`, `atol_p=1e-12`, and `min_step=1e-12 M`.
- Error uses the v3 component scale and RMS norm; acceptance is `error<=1`.
- Default maximum rejection count is 12 and maximum total steps is 2,000,000.
- Ordinary normalized Hamiltonian error must remain below `1e-10`.
- Dense event roots use a bracketed safeguarded method, never pure Newton.
- External geometric-event root tolerance is `1e-10 M`.
- No default constraint projection is allowed.
- No observer, tetrad, Carter constant, separated solver, Kerr–Schild,
  radiative transfer, renderer, or ray CLI is added.
- Existing public N-body and dynamic-vector integrator signatures remain
  compatible.

---

### Task 1: Fixed-size DOPRI5 step and controller

**Files:**
- Create: `include/solar/numerics/dopri5.h`
- Create: `tests/relativity/test_dopri5.cpp`

**Interfaces:**
- Consumes: a fixed `std::array<double,N>` state and RHS callable.
- Produces: `StateN`, `ErrorNorm`, `Dopri5Config`,
  `Dopri5StepResult`, and `dopri5_step`.

- [x] **Step 1: Read the test-quality rules before changing tests**

Read completely:

```bash
sed -n '1,400p' \
  /Users/hostsjim/.codex/plugins/cache/openai-curated-remote/superpowers/6.2.0/skills/test-driven-development/writing-good-tests.md
```

Apply its mutation question to every new assertion: identify the production
expression whose change would make that assertion fail.

- [x] **Step 2: Write the missing fixed-size API test**

Create `tests/relativity/test_dopri5.cpp` with a local `check` counter and:

```cpp
using solar::numerics::Dopri5Config;
using solar::numerics::ErrorNorm;
using solar::numerics::StateN;

Dopri5Config<1> config;
config.absolute_tolerance = StateN<1>{{1.0e-12}};
config.relative_tolerance = 1.0e-11;
config.safety = 0.9;
config.min_factor = 0.2;
config.max_factor = 5.0;
config.error_norm = ErrorNorm::RootMeanSquare;

const StateN<1> initial{{1.0}};
const auto exponential_rhs =
    [](double, const StateN<1>& state) { return state; };
const auto step = solar::numerics::dopri5_step(
    initial, 0.0, 0.01, exponential_rhs, config);

check("fixed DOPRI step accepted", step.accepted);
check_near("fifth-order exponential step",
           step.state[0], std::exp(0.01), 3.0e-14);
check("accepted next step preserves direction", step.next_step > 0.0);
```

Add separate assertions that:

- step `1.0` with `atol=rtol=1e-14` is rejected;
- its rejected `abs(next_step)` is smaller than `abs(step_used)`;
- a negative step produces negative `next_step`;
- zero/non-finite step, non-positive tolerance, non-finite state, invalid
  safety, and factors outside `0<min_factor<=1<=max_factor` do not report a
  successful completed step.

- [x] **Step 3: Run the test and verify RED**

```bash
make tests/relativity/test_dopri5
```

Expected: compilation fails because `solar/numerics/dopri5.h` does not exist.

- [x] **Step 4: Implement the shared stage engine and fixed wrapper**

In `include/solar/numerics/dopri5.h`, define the public API from the design and
one `detail::dopri5_step_impl<State,Rhs>` used by the fixed wrapper.

Use the standard coefficients:

```text
a21=1/5
a31=3/40, a32=9/40
a41=44/45, a42=-56/15, a43=32/9
a51=19372/6561, a52=-25360/2187,
    a53=64448/6561, a54=-212/729
a61=9017/3168, a62=-355/33, a63=46732/5247,
    a64=49/176, a65=-5103/18656
b1=35/384, b3=500/1113, b4=125/192,
    b5=-2187/6784, b6=11/84
e1=71/57600, e3=-71/16695, e4=71/1920,
    e5=-17253/339200, e6=22/525, e7=-1/40
c2=1/5, c3=3/10, c4=4/5, c5=8/9, c6=1
```

For every component compute:

```cpp
const double scale =
    absolute_tolerance[i] +
    relative_tolerance *
        std::max(std::fabs(state[i]), std::fabs(trial[i]));
const double normalized = error_estimate[i] / scale;
```

RMS mode is:

```cpp
error = std::sqrt(sum_squared / static_cast<double>(state.size()));
```

Maximum mode is the maximum absolute normalized component.

Return explicit status for non-finite input/stage output. Invalid configuration
throws `std::invalid_argument`. Compute:

```cpp
factor = safety * std::pow(1.0 / error, 0.2);
factor = std::clamp(factor, min_factor, max_factor);
if (error > 1.0) {
    factor = std::min(factor, 1.0);
}
next_step = std::copysign(std::fabs(step) * factor, step);
```

Do not implement dense output in this task.

- [x] **Step 5: Verify GREEN**

```bash
make tests/relativity/test_dopri5
./tests/relativity/test_dopri5
```

Expected: all fixed-step/controller assertions pass with no new compiler
warning.

- [x] **Step 6: Commit**

```bash
git add include/solar/numerics/dopri5.h \
  tests/relativity/test_dopri5.cpp
git commit -m "feat: add fixed-size DOPRI5 kernel"
```

### Task 2: Fourth-order dense output

**Files:**
- Modify: `include/solar/numerics/dopri5.h`
- Modify: `tests/relativity/test_dopri5.cpp`

**Interfaces:**
- Consumes: the seven finite DOPRI5 stages from Task 1.
- Produces: `Dopri5DenseOutput<N>` and optional dense output in each completed
  step.

- [x] **Step 1: Add failing dense-output tests**

For `y'=y`, accepted step `h=0.2`, assert:

```cpp
const auto dense = step.dense_output.value();
check_near("dense start", dense.evaluate(0.0)[0], 1.0, 0.0);
check_near("dense midpoint", dense.evaluate(0.1)[0],
           std::exp(0.1), 2.0e-7);
check_near("dense end", dense.evaluate(0.2)[0],
           step.state[0], 2.0e-15);
```

Require `evaluate(-1e-6)` and `evaluate(0.200001)` to throw
`std::out_of_range`. Repeat with a negative step and require its endpoints and
midpoint to use the decreasing independent-variable interval correctly.

- [x] **Step 2: Run and verify RED**

```bash
make tests/relativity/test_dopri5
```

Expected: compilation fails because `Dopri5DenseOutput`/`dense_output` is
absent.

- [x] **Step 3: Implement the Dormand–Prince continuous extension**

Store `t0`, `h`, `y0`, and four coefficient vectors `q1..q4`. Form them from
the stage derivatives and the Shampine matrix:

```text
P1 = [1,
      -8048581381/2820520608,
       8663915743/2820520608,
      -12715105075/11282082432]
P2 = [0, 0, 0, 0]
P3 = [0, 131558114200/32700410799,
         -68118460800/10900136933,
          87487479700/32700410799]
P4 = [0, -1754552775/470086768,
          14199869525/1410260304,
         -10690763975/1880347072]
P5 = [0, 127303824393/49829197408,
         -318862633887/49829197408,
          701980252875/199316789632]
P6 = [0, -282668133/205662961,
          2019193451/616988883,
         -1453857185/822651844]
P7 = [0, 40617522/29380423,
         -110615467/29380423,
          69997945/29380423]
```

For normalized `theta=(t-t0)/h`:

```cpp
y(t) = y0 + h * (
    q1*theta + q2*theta*theta +
    q3*theta*theta*theta +
    q4*theta*theta*theta*theta);
```

Allow only the closed interval between `t0` and `t0+h`, independent of step
sign. Construct dense output only when every stage and trial value is finite.

- [x] **Step 4: Verify GREEN and run a mutation**

```bash
make tests/relativity/test_dopri5
./tests/relativity/test_dopri5
```

Temporarily negate one nonzero continuous-extension coefficient and rerun.
Expected: the midpoint assertion fails. Restore the coefficient and rerun to
green before committing.

- [x] **Step 5: Commit**

```bash
git add include/solar/numerics/dopri5.h \
  tests/relativity/test_dopri5.cpp
git commit -m "feat: add DOPRI5 dense output"
```

### Task 3: Preserve the legacy generic integrator API

**Files:**
- Modify: `src/integrator.cpp`
- Create: `tests/test_integrator.cpp`
- Modify: `include/solar/numerics/dopri5.h`

**Interfaces:**
- Consumes: Task 1's container-generic internal stage engine.
- Produces: unchanged `dopri5_generic_step` public behavior without duplicated
  generic stage equations.

- [x] **Step 1: Add the pre-refactor golden compatibility test**

Create `tests/test_integrator.cpp`. For:

```cpp
std::vector<double> initial{1.0, -2.0};
auto rhs = [](double t, const std::vector<double>& y) {
    return std::vector<double>{
        y[0] + t,
        -0.5 * y[1] + 2.0 * t,
    };
};
const auto result = solar::dopri5_generic_step(
    initial, 0.3, 0.1, rhs, 1.0e-9, 1.0e-9);
```

Require the existing literals:

```text
state[0] = 1.1418931121666667
state[1] = -1.834098762375
error    = 8.3354999831097523
next     = 0.058891983053898568
accepted = false
```

Use tolerances of `3e-15` for state, `2e-13` for error, and `2e-15` for next
step. Run the test before refactoring and record that it is green; this is a
characterization gate, not a new behavior claim.

- [x] **Step 2: Expose a vector-compatible detail adapter**

Make the internal stage engine operate on any contiguous array-like state with
`size()` and indexed access. It must create work states by copying the input,
so both fixed arrays and sized vectors retain their shape.

Keep the public fixed-size API unchanged.

- [x] **Step 3: Replace only `dopri5_generic_step` internals**

In `src/integrator.cpp`:

- include `solar/numerics/dopri5.h`;
- remove the duplicate dynamic-vector generic stages;
- construct a vector of uniform `atol`;
- invoke the shared detail engine with `ErrorNorm::Maximum`;
- adapt state/error/accepted/next-step into `GenericAdaptiveResult`.

Do not change `dopri5_step` for N-body `State` objects or any signature in
`include/solar/integrator.h`.

- [x] **Step 4: Verify compatibility**

```bash
make tests/test_integrator tests/test_validation
./tests/test_integrator
./tests/test_validation
```

Expected: the golden literals and all legacy validation assertions pass.

- [x] **Step 5: Commit**

```bash
git add include/solar/numerics/dopri5.h \
  src/integrator.cpp tests/test_integrator.cpp
git commit -m "refactor: share generic DOPRI5 stages"
```

### Task 4: Hamiltonian and general geodesic RHS

**Files:**
- Create: `include/solar/relativity/hamiltonian.h`
- Create: `src/relativity/hamiltonian.cpp`
- Create: `tests/relativity/test_hamiltonian.cpp`

**Interfaces:**
- Consumes: `Metric`, `PhaseSpaceState`, `GeodesicKind`, and `StateN<8>`.
- Produces: `PhaseSpaceDerivative`, Hamiltonian/constraint functions,
  pack/unpack, and `HamiltonGeodesicRhs`.

- [x] **Step 1: Write exact Minkowski Hamiltonian tests**

Add:

```cpp
const MinkowskiMetric metric;
const PhaseSpaceState photon{
    0.0,
    Contravariant4{Vec4{{0.0, 0.0, 0.0, 0.0}}},
    Covariant4{Vec4{{-1.0, 1.0, 0.0, 0.0}}},
};
check_near("null Hamiltonian",
           hamiltonian(metric, photon), 0.0, 0.0);
const auto photon_rhs = HamiltonGeodesicRhs(metric)(photon);
check_near("future time tangent", photon_rhs.dx.v[0], 1.0, 0.0);
check_near("null spatial tangent", photon_rhs.dx.v[1], 1.0, 0.0);
check("Minkowski momentum derivative zero",
      max_norm(photon_rhs.dp.v) == 0.0);
```

For velocity `v=0.6`, use `gamma=1.25` and
`p=(-1.25,0.75,0,0)`; require `H=-0.5` exactly.

Add pack/unpack round-trip and a hand-computed normalized constraint error for
`p=(-2,2.001,0,0)`.

- [x] **Step 2: Add the no-false-symmetry test**

In the test, define a minimal `TimeDependentMetric`:

```cpp
g^tt = -(1 + 0.1*t)
partial_t g^tt = -0.1
```

with covariant inverse diagonal and finite validity. At `p_t=-1`, require:

```text
dp_t = +0.05
```

Also require exact `dp_t=dp_phi=0` for Schwarzschild and Kerr BL at valid
points.

- [x] **Step 3: Run and verify RED**

```bash
make tests/relativity/test_hamiltonian
```

Expected: compilation fails because `hamiltonian.h` is absent.

- [x] **Step 4: Implement the exact contractions**

Use two explicit four-dimensional loops:

```cpp
H += 0.5 * inverse[row][column] *
     state.p.v[row] * state.p.v[column];

dx.v[row] += inverse[row][column] * state.p.v[column];

dp.v[mu] -= 0.5 * derivatives[mu][alpha][beta] *
            state.p.v[alpha] * state.p.v[beta];
```

Compute the constraint denominator independently from `H` using the absolute
value of every matrix-momentum term as required by v3. Reject non-finite
affine/state/momentum and invalid metric points; do not clamp or return a
fallback derivative.

- [x] **Step 5: Verify GREEN and derivative mutation**

```bash
make tests/relativity/test_hamiltonian
./tests/relativity/test_hamiltonian
```

Temporarily remove the `0.5` factor from `dp` and require the time-dependent
metric assertion to fail. Restore and rerun green.

- [x] **Step 6: Commit**

```bash
git add include/solar/relativity/hamiltonian.h \
  src/relativity/hamiltonian.cpp \
  tests/relativity/test_hamiltonian.cpp
git commit -m "feat: add Hamiltonian geodesic RHS"
```

### Task 5: Event contracts and bracketed dense root

**Files:**
- Create: `include/solar/relativity/geodesic_types.h`
- Create: `include/solar/relativity/event_root.h`
- Create: `src/relativity/event_root.cpp`
- Create: `tests/relativity/test_geodesic_events.cpp`

**Interfaces:**
- Consumes: `Dopri5DenseOutput<8>` and event functions over
  `PhaseSpaceState`.
- Produces: exact v3 event/termination/diagnostic types and
  `locate_event`.

- [x] **Step 1: Write event-direction and endpoint tests**

Build exact dense output from an accepted constant RHS step:

```cpp
StateN<8> initial{};
initial[1] = 0.0;
auto rhs = [](double, const StateN<8>&) {
    StateN<8> derivative{};
    derivative[1] = 1.0;
    return derivative;
};
```

For event `x^1-0.3`, require:

- `Increasing` finds affine `0.3` within `1e-12`;
- `Decreasing` returns `NoRoot`;
- `Any` finds the same root;
- exact roots at dense start and end are returned;
- a negative step from `x^1=1` to `0` finds `x^1=0.3` with
  `Decreasing`, defined along step progression.

Use a function that returns NaN only at internal trial points to require
`EventRootStatus::Failed`, distinct from `NoRoot`.

- [x] **Step 2: Run and verify RED**

```bash
make tests/relativity/test_geodesic_events
```

Expected: missing geodesic/event headers.

- [x] **Step 3: Implement L0 types**

Copy every `TerminationReason` field and every diagnostics field from v3.
Initialize unavailable invariant/min-step fields with quiet NaN. Define
`GeodesicEvent`, `EventHit`, `EventRootStatus`, and `EventRootResult` exactly as
the design.

- [x] **Step 4: Implement safeguarded secant/bisection**

Represent the bracket by dense fraction `theta in [0,1]`. At each iteration:

```cpp
theta_secant =
    right - f_right * (right-left) / (f_right-f_left);
```

Use it only when finite and strictly inside the middle 80% of the bracket;
otherwise use the midpoint. Preserve the sign bracket after each finite
evaluation. Stop when:

```cpp
std::fabs(dense.end() - dense.start()) *
    (right - left) <= event.root_tolerance
```

Use at most 100 iterations. Return `Failed` on non-finite values or exhausted
iterations. Do not use a derivative/Newton update.

- [x] **Step 5: Verify GREEN**

```bash
make tests/relativity/test_geodesic_events
./tests/relativity/test_geodesic_events
```

Expected: all direction, endpoint, negative-step, and failure distinctions
pass.

- [x] **Step 6: Commit**

```bash
git add include/solar/relativity/geodesic_types.h \
  include/solar/relativity/event_root.h \
  src/relativity/event_root.cpp \
  tests/relativity/test_geodesic_events.cpp
git commit -m "feat: add directed geodesic events"
```

### Task 6: Adaptive geodesic flow and Minkowski validation

**Files:**
- Create: `include/solar/relativity/geodesic_integrator.h`
- Create: `src/relativity/geodesic_config.cpp`
- Create: `src/relativity/geodesic_config_internal.h`
- Create: `src/relativity/geodesic_step_attempt.cpp`
- Create: `src/relativity/geodesic_step_attempt.h`
- Create: `src/relativity/geodesic_integrator.cpp`
- Create: `tests/relativity/test_geodesics.cpp`
- Create: `tests/relativity/test_geodesic_failures.cpp`

**Interfaces:**
- Consumes: Tasks 1–5.
- Produces: `GeodesicIntegrationConfig`,
  `GeodesicIntegrationResult`, and `GeodesicIntegrator`.

- [x] **Step 1: Write Minkowski null and timelike line tests**

For the null state `p=(-1,1,0,0)`, use:

```cpp
auto config = GeodesicIntegrationConfig::cpu_reference(
    GeodesicKind::Null, 1.0, 0.1, 0.5, 10.0);
const auto result = GeodesicIntegrator(metric).integrate(
    photon, config);
```

Require:

```text
reason = MaxAffine
affine = 10
x0 = 10
x1 = 10
p unchanged
max_constraint_error < 1e-14
```

For timelike `p=(-1.25,0.75,0,0)`, integrate affine `4` and require
`x0=5`, `x1=3`, `H=-1/2`, and `MaxProperTime` when a proper-time limit of `4`
is configured.

- [x] **Step 2: Write limit, event, and reversibility tests**

Add:

- forward null integration for affine `3`, then a new config with negative
  `initial_step` and affine displacement `3`; require every canonical
  component returns within `2e-12`;
- user event `x1-2.25` with `Increasing`; require returned state/affine
  `2.25` and its declared reason;
- two events at `x1=1.5` and `x1=1.0` in one large accepted step; require the
  earlier step-progression root;
- exact `MaxSteps`, `MaxAffine`, `MaxProperTime`, and `MaxCoordinateTime`
  reasons;
- invalid initial constraint, NaN state, zero/underflowing step, and invalid
  config return/throw according to the design.

- [x] **Step 3: Run and verify RED**

```bash
make tests/relativity/test_geodesics
```

Expected: `geodesic_integrator.h` is absent.

- [x] **Step 4: Implement config validation and CPU defaults**

Set per-component tolerances:

```cpp
absolute_tolerance[0] = 1.0e-11 * mass_scale;
for (std::size_t i = 1; i < 4; ++i) {
    absolute_tolerance[i] = 1.0e-12 * mass_scale;
}
for (std::size_t i = 4; i < 8; ++i) {
    absolute_tolerance[i] = 1.0e-12;
}
```

Set the exact controller, step, rejection, total-step, and constraint defaults
from Global Constraints. Reject null `max_proper_time`.

- [x] **Step 5: Implement the accepted/rejected loop**

Use a small private adapter that unpacks `(affine,state8)`, invokes
`HamiltonGeodesicRhs`, and repacks the derivative.

On each attempt:

1. preserve the sign of `initial_step`;
2. cap magnitude by `max_step` and remaining affine displacement;
3. invoke fixed DOPRI5;
4. increment rejected count on error rejection or invalid trial metric point;
5. shrink without accepting the trial;
6. stop at configured rejection/underflow conditions;
7. on acceptance, locate all events and select the smallest dense fraction;
8. evaluate exact dense roots for proper/coordinate-time limits;
9. update constraint and diagnostic maxima;
10. return one explicit reason/message.

Do not allocate a trajectory vector and do not project momentum.

- [x] **Step 6: Verify GREEN**

```bash
make tests/relativity/test_geodesics
./tests/relativity/test_geodesics
```

Expected: all analytic lines, limits, events, error paths, and reversibility
pass.

- [x] **Step 7: Commit**

```bash
git add include/solar/relativity/geodesic_integrator.h \
  src/relativity/geodesic_integrator.cpp \
  tests/relativity/test_geodesics.cpp
git commit -m "feat: integrate Hamiltonian geodesics"
```

### Task 7: Schwarzschild analytic validation

**Files:**
- Create: `tests/relativity/test_geodesics_schwarzschild.cpp`
- Modify only for a demonstrated test failure:
  `src/relativity/geodesic_integrator.cpp`
- Modify only for a demonstrated test failure:
  `include/solar/numerics/dopri5.h`

**Interfaces:**
- Consumes: the completed general flow.
- Produces: independent radial-null, weak-bending, and photon-sphere evidence.

- [x] **Step 1: Add outgoing radial-null analytic test**

At `M=1`, `r0=10`, equator, use:

```text
p_t = -1
p_r = 1/f(r0) = 1.25
p_theta = p_phi = 0
```

Integrate to affine displacement `2`. Require:

```text
r = 12
t = 2 + 2*log((12-2)/(10-2))
```

within `2e-10`, and maximum normalized Hamiltonian error below `1e-10`.

- [x] **Step 2: Add photon-sphere RHS test**

At equatorial `r=3M` with:

```text
p_t=-1
p_r=0
p_theta=0
p_phi=3*sqrt(3)*M
```

Require `H`, `dr/dlambda`, and `dp_r/dlambda` below `2e-14`.

- [x] **Step 3: Add weak-field deflection test**

At `M=1`, impact parameter `b=100`, start/end radius `R=10000`, equator:

```cpp
const Mat4 inverse = metric.contravariant(initial.x);
const double base =
    inverse[0][0] +
    -2.0 * b * inverse[0][3] +
    b * b * inverse[3][3];
initial.p.v = Vec4{{
    -1.0,
    -std::sqrt(-base / inverse[1][1]),
    0.0,
    b,
}};
```

Use an `Increasing` event for `r-R`; the initial decreasing direction prevents
the start root from firing. Compute finite-distance deflection:

```text
deflection =
abs(phi_final-phi_initial) -
(pi - 2*asin(b/R))
```

Require it within 5% of `4M/b`, an ordinary-ray Hamiltonian error below
`1e-10`, and explicit `Escaped` termination.

- [x] **Step 4: Run and correct only demonstrated failures**

```bash
make tests/relativity/test_geodesics
./tests/relativity/test_geodesics
```

If a threshold fails, inspect convergence under halved `max_step` before
changing production. Do not loosen the v3 Hamiltonian gate.

- [x] **Step 5: Commit**

```bash
git add tests/relativity/test_geodesics_schwarzschild.cpp \
  src/relativity/geodesic_integrator.cpp \
  include/solar/numerics/dopri5.h
git commit -m "test: validate Schwarzschild geodesics"
```

Stage only production files that were actually changed for a demonstrated
failure.

### Task 8: Kerr invariants and transparent diagnostics

**Files:**
- Create: `tests/relativity/test_geodesics_kerr.cpp`
- Modify only for a demonstrated test failure:
  `src/relativity/geodesic_integrator.cpp`

**Interfaces:**
- Consumes: Kerr BL and explicit invariant-monitor flags.
- Produces: E/Lz drift evidence and unavailable-Carter semantics.

- [x] **Step 1: Add an ordinary Kerr null state**

At `M=1`, `chi=0.7`, `r=10`, `theta=1.2`, choose
`p_t=-1`, `p_phi=2`, `p_theta=0`. Solve only `p_r` from the independent
quadratic:

```cpp
const Mat4 inverse = metric.contravariant(initial.x);
const double nonradial =
    inverse[0][0] * p_t * p_t +
    2.0 * inverse[0][3] * p_t * p_phi +
    inverse[3][3] * p_phi * p_phi;
const double p_r =
    -std::sqrt(-nonradial / inverse[1][1]);
```

Enable `monitor_energy` and `monitor_lz`, integrate a short ordinary exterior
segment, and require:

```text
max_energy_rel_error = 0
max_lz_rel_error = 0
max_carter_rel_error = NaN
max_constraint_error < 1e-10
```

- [x] **Step 2: Add unmonitored invariant semantics**

For a Minkowski integration with both flags false, require energy, Lz, and
Carter diagnostic fields are all NaN. For a monitored exactly-zero Lz, require
the recorded field uses absolute drift and remains zero.

- [x] **Step 3: Run and verify**

```bash
make tests/relativity/test_geodesics
./tests/relativity/test_geodesics
```

Expected: all invariant/NaN semantics and the ordinary Hamiltonian gate pass.

- [x] **Step 4: Commit**

```bash
git add tests/relativity/test_geodesics_kerr.cpp \
  src/relativity/geodesic_integrator.cpp
git commit -m "test: validate Kerr geodesic invariants"
```

Stage the production file only if it changed for a demonstrated failure.

### Task 9: Phase 1 audit, validation, and gate

**Files:**
- Create: `docs/validation/relativity_01_hamiltonian_geodesics.md`
- Modify: `RELATIVITY_STATUS.md`
- Modify: `docs/superpowers/plans/2026-07-29-relativity-phase-1.md`

**Interfaces:**
- Consumes: all actual build/test/numerical results.
- Produces: the reproducible Phase 1 evidence and gate state.

- [ ] **Step 1: Run the complete release verification**

```bash
make clean
make
make test
./tests/relativity/test_dopri5
./tests/relativity/test_hamiltonian
./tests/relativity/test_geodesic_events
./tests/relativity/test_geodesics
git diff --check
```

Run an invalid-metric/constraint test separately and confirm the process exits
nonzero when its expected reason is deliberately changed.

- [ ] **Step 2: Run sanitizer verification**

Clean and rebuild the complete relativity test set with:

```bash
CXXFLAGS='-std=c++17 -O1 -g -Wall -Wextra -Iinclude \
  -fsanitize=address,undefined -fno-omit-frame-pointer'
```

Run every `tests/relativity/test_*` executable with:

```bash
ASAN_OPTIONS=detect_leaks=0:halt_on_error=1
UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1
```

Afterward, clean and rerun the normal release build/test so the worktree does
not retain sanitizer artifacts.

- [ ] **Step 3: Audit the complete Phase 1 diff**

Check line by line:

- Hamilton equations and normalized constraint denominator;
- canonical position/momentum variance and state order;
- DOPRI coefficients, error norm, sign preservation, and controller config;
- dense extension coefficients and endpoint behavior;
- bracket/direction/root-failure semantics;
- every termination reason and unavailable diagnostic;
- trial-domain rejection versus physical termination;
- no default projection;
- no unconditional E/Lz conservation;
- existing integrator API compatibility;
- no tetrad/observer/Phase 2 code;
- no files with mixed independent responsibilities that need an authorized
  split.

- [ ] **Step 4: Write the validation report**

Record:

- verified code commit and platform;
- every command actually run;
- test assertion totals and external skips;
- DOPRI convergence/dense/event errors;
- maximum Minkowski, Schwarzschild, and Kerr Hamiltonian errors;
- radial-null, weak-bending, reversibility, E/Lz results;
- model boundary and the v3 section 7.5 deferral;
- all remaining warnings/unverified platforms;
- at least three most likely bugs;
- fastest commands/parameter mutations to falsify the gate.

- [ ] **Step 5: Mark the phase gate**

Only after every required command succeeds, set:

```text
CURRENT_PHASE: 1
PHASE_STATE: PASSED
NEXT_ALLOWED_ACTION: Phase 2 only
```

Mark every plan checkbox complete.

- [ ] **Step 6: Commit the gate**

```bash
git add docs/validation/relativity_01_hamiltonian_geodesics.md \
  RELATIVITY_STATUS.md \
  docs/superpowers/plans/2026-07-29-relativity-phase-1.md
git commit -m "docs: pass relativity phase 1 gate"
```

- [ ] **Step 7: Publish a stacked draft PR**

Push `codex/relativity-phase-1` and create a draft PR targeting
`codex/relativity-phase-0b`, not `main`, while Phase 0B PR #3 remains open.
The PR body must state the Phase 0B dependency, model boundary, exact
verification totals, numerical maxima, skips, warnings, and next allowed
phase.
