# Solar Relativity Phase 4 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a validated Cartesian ingoing Kerr–Schild chart that evolves null and unit-mass timelike Hamiltonian geodesics through the Kerr outer horizon to an explicit finite interior cutoff.

**Architecture:** Keep coordinate-independent state and algebra in L0, place the scalar-generic Kerr–Schild fields, metric, chart transform, invariant evaluators, and standard events in L1, and reuse the existing generic L2 `GeodesicIntegrator`. Production inverse-metric derivatives use `Dual4`; independent analytic radius derivatives, five-point tensor differences, BL/KS common-event evolution, and horizon continuation provide separate validation paths.

**Tech Stack:** C++17, Solar `Metric`/Hamiltonian/DOPRI5/event APIs, `Dual4` automatic differentiation, CMake/Make, AddressSanitizer, UndefinedBehaviorSanitizer.

## Global Constraints

- Follow `docs/superpowers/specs/2026-07-30-relativity-phase-4-design.md` and `SOLAR_RELATIVITY_KERR_完整实施主提示词_v3.md`.
- Preserve signature `(-,+,+,+)`, coordinate order `(t,x,y,z)`, state order `(x^mu,p_mu)`, covariant canonical momentum, and geometrized units `G=c=1`.
- Accept only finite `M>0` and finite subextremal `|chi|<1`; do not add extremal transforms or a negative-`r` extension.
- Keep the dependency direction `L3 -> L2 -> L1 -> L0`; do not add a Kerr–Schild-specific integrator.
- Existing BL invariant monitoring and existing public calls must behave identically when new callbacks are absent.
- Metric inverse gates are `<5e-13` ordinarily and `<1e-10` at safe near-horizon points.
- Kerr–Schild inverse derivative versus five-point difference is relative `<3e-8`; analytic radius gradient versus independent `Dual4` is relative `<1e-12`.
- Transform Jacobian identity is `<5e-12`; canonical round trip and Hamiltonian/pairing invariance are `<1e-10`.
- Ordinary BL/KS common-event position and momentum P95 are `<1e-8`; safe near-horizon errors are `<1e-6`.
- Ordinary normalized Hamiltonian error is `<1e-10`; energy and axial angular momentum relative drift are `<1e-12`.
- Horizon freely falling tetrad orthonormality is `<1e-12`; interior termination occurs at `0.05M` with `InteriorCutoff`.
- Do not add radiative transfer, material models, rendering, WASM, GPU work, long-time symplectic integration, or Cartesian Carter evaluation in this phase.
- Use test-first changes, one independently reviewable commit per task, and do not weaken a physical gate to make a test pass.

## File Map

| File | Layer | Responsibility |
|---|---|---|
| `include/solar/relativity/geodesic_integrator.h` | L1/L2 boundary | Append optional chart-aware energy and Lz evaluators to integration configuration. |
| `src/relativity/geodesic_invariant_monitor.h/.cpp` | L1 | Evaluate enabled invariant callbacks with explicit failure semantics while preserving BL defaults. |
| `include/solar/relativity/kerr_schild_metric.h` | L1 public | Public Cartesian KS metric, radius, gradient, horizon accessors, and Cartesian conserved quantities. |
| `src/relativity/kerr_schild_fields.h` | L1 internal | Scalar-generic stable radius, null one-form, scalar `H`, metric, and inverse expressions. |
| `src/relativity/kerr_schild_metric.cpp` | L1 | Validate parameters/points and expose double/`Dual4` metric evaluation. |
| `include/solar/relativity/kerr_chart_transform.h` | L1 public | Safe-overlap BL↔KS position, Jacobian, and full canonical-state transforms. |
| `src/relativity/kerr_chart_fields.h` | L1 internal | Scalar-generic logarithmic offsets and forward coordinate expression. |
| `src/relativity/kerr_chart_transform.cpp` | L1 | Safe-overlap validation, position maps, and forward/inverse Jacobians. |
| `src/relativity/kerr_chart_state_transform.cpp` | L1 | Full canonical covector and affine-preserving state transforms. |
| `include/solar/relativity/kerr_schild_events.h` | L1 public | Standard outer-horizon and finite-interior event factories. |
| `src/relativity/kerr_schild_events.cpp` | L1 | Validate event configuration and capture immutable metric values safely. |
| `tests/relativity/test_geodesic_invariant_callbacks.cpp` | verification | Callback selection, drift, exceptions, and non-finite values. |
| `tests/relativity/test_kerr_schild.cpp` | verification | Metric algebra, stable radius, domain, Schwarzschild limit, and conserved quantities. |
| `tests/relativity/test_kerr_schild_derivatives.cpp` | verification | Independent radius and inverse-metric derivative checks. |
| `tests/relativity/test_kerr_chart_transform.cpp` | verification | Position/Jacobian/canonical round trips and rejection behavior. |
| `tests/relativity/test_kerr_schild_events.cpp` | verification | Event radius, direction, reason, ownership, and input validation. |
| `tests/relativity/test_kerr_bl_ks_crosscheck.cpp` | verification | Ordinary and near-horizon common-event BL/KS equivalence. |
| `tests/relativity/test_geodesics_kerr_schild.cpp` | verification | Timelike horizon crossing, tetrad validity, restart, and interior cutoff. |
| `tests/external_consumer/probe.cpp` | L3 verification | Prove the installed package exposes and links the new public metric. |
| `docs/validation/relativity_09_kerr_schild.md` | documentation | Record formulas, fixtures, commands, measured gates, scope, and residual risks. |
| `RELATIVITY_STATUS.md` | documentation | Mark Phase 4 complete only after all repository gates pass. |

Large existing files intentionally deferred: separated-geodesic sources, renderer/shadow tests, CLI routing, and legacy Solar dynamics are not on the Phase 4 call path.

---

### Task 1: Chart-aware invariant evaluators

**Files:**
- Create: `tests/relativity/test_geodesic_invariant_callbacks.cpp`
- Modify: `include/solar/relativity/geodesic_integrator.h`
- Modify: `src/relativity/geodesic_invariant_monitor.h`
- Modify: `src/relativity/geodesic_invariant_monitor.cpp`

**Interfaces:**
- Consumes: `using InvariantEvaluator = std::function<double(const PhaseSpaceState&)>`.
- Produces: appended `GeodesicIntegrationConfig::stationary_energy_evaluator` and `GeodesicIntegrationConfig::axial_angular_momentum_evaluator`.
- Preserves: empty evaluators select `-state.p.v[0]` and `state.p.v[3]`.

- [x] **Step 1: Write callback behavior tests**

Create a flat-metric integration fixture whose custom energy is
`state.x.v[1] + state.p.v[0]` and custom Lz is
`state.x.v[2] - 2.0 * state.p.v[3]`. Exercise:

```cpp
config.monitor_energy = true;
config.monitor_lz = true;
config.stationary_energy_evaluator =
    [](const PhaseSpaceState& state) {
        return state.x.v[1] + state.p.v[0];
    };
config.axial_angular_momentum_evaluator =
    [](const PhaseSpaceState& state) {
        return state.x.v[2] - 2.0 * state.p.v[3];
    };
```

Require the diagnostics to reflect those callbacks rather than the BL
defaults. Add separate integrations where each callback throws
`std::runtime_error("synthetic evaluator failure")` or returns NaN; require
`TerminationReason::NonFiniteState` and an explanatory diagnostic message.
Finally run with empty evaluators and require existing energy/Lz diagnostics
to retain their BL values.

- [x] **Step 2: Run the focused test and confirm the red state**

Run:

```bash
make -j4 tests/relativity/test_geodesic_invariant_callbacks
```

Expected: compilation fails because the two configuration fields do not
exist.

- [x] **Step 3: Add evaluator fields and one shared checked-evaluation path**

Append to `GeodesicIntegrationConfig`:

```cpp
InvariantEvaluator stationary_energy_evaluator;
InvariantEvaluator axial_angular_momentum_evaluator;
```

Store both callbacks in `GeodesicInvariantMonitor`. Add:

```cpp
bool evaluate(
    const char* label,
    const InvariantEvaluator& evaluator,
    const PhaseSpaceState& state,
    double& value,
    std::string& failure_message) const;
```

The method catches standard and non-standard exceptions and rejects
non-finite results. During initialization and update, evaluate a callback
only when its corresponding monitor flag is true; otherwise do no work.
When the callback is empty, evaluate the existing BL formula directly.
Use the same checked method for Carter so all callback failures have the same
termination semantics.

- [x] **Step 4: Run focused and existing invariant tests**

Run:

```bash
make -j4 tests/relativity/test_geodesic_invariant_callbacks \
  tests/relativity/test_geodesics_kerr \
  tests/relativity/test_geodesic_failures
./tests/relativity/test_geodesic_invariant_callbacks
./tests/relativity/test_geodesics_kerr
./tests/relativity/test_geodesic_failures
```

Expected: all executables report zero failures.

- [x] **Step 5: Commit**

```bash
git add include/solar/relativity/geodesic_integrator.h \
  src/relativity/geodesic_invariant_monitor.h \
  src/relativity/geodesic_invariant_monitor.cpp \
  tests/relativity/test_geodesic_invariant_callbacks.cpp
git commit -m "feat(relativity): support chart-aware invariant monitors"
```

### Task 2: Stable Cartesian Kerr–Schild metric and production derivatives

**Files:**
- Create: `include/solar/relativity/kerr_schild_metric.h`
- Create: `src/relativity/kerr_schild_fields.h`
- Create: `src/relativity/kerr_schild_metric.cpp`
- Create: `tests/relativity/test_kerr_schild.cpp`

**Interfaces:**
- Consumes: `Metric`, `Mat4`, `Vec3`, `Dual4`, and `PhaseSpaceState`.
- Produces: `KerrSchildCartesianMetric` and the two Cartesian invariant helper functions declared in the design.

- [x] **Step 1: Write metric value and domain tests**

For `M=2`, `chi=0.6`, sample generic exterior, axis, outer-horizon,
just-inside-horizon, and interior points. Reconstruct

```cpp
const double null_norm =
    -l[0] * l[0] + l[1] * l[1] +
    l[2] * l[2] + l[3] * l[3];
const Mat4 identity = multiply(
    metric.covariant(x), metric.contravariant(x));
```

Require exact symmetry, `abs(null_norm)<5e-13`, ordinary identity error
`<5e-13`, and near-horizon identity error `<1e-10`. Check the implicit
quartic residual
`r^4-(rho^2-a^2)r^2-a^2z^2` with scale-aware tolerance. For `chi=0`,
require the analytic Schwarzschild KS form
`g=eta+2(M/r)l⊗l`, `l=(1,x/r,y/r,z/r)`.

Require `chart()==Chart::KerrSchildCartesian`, finite horizon accessors,
axis/horizon/interior validity, and explicit rejection of the ring,
zero-radius Schwarzschild point, non-finite coordinates, invalid mass, and
extremal/superextremal spin. Verify:

```cpp
check_near("stationary energy",
           kerr_schild_stationary_energy(state),
           -state.p.v[0], 0.0);
check_near("axial angular momentum",
           kerr_schild_axial_angular_momentum(state),
           state.x.v[1] * state.p.v[2] -
               state.x.v[2] * state.p.v[1],
           0.0);
```

Also require the Schwarzschild radius gradient to equal
`(x/r,y/r,z/r)`, every inverse-metric derivative to be finite, symmetric,
and its stationary time derivative to be exactly zero. These basic
behavior assertions prevent an incomplete concrete `Metric` implementation;
Task 3 supplies independent precision validation.

- [x] **Step 2: Run the focused test and confirm the red state**

Run:

```bash
make -j4 tests/relativity/test_kerr_schild
```

Expected: compilation fails because
`solar/relativity/kerr_schild_metric.h` is absent.

- [x] **Step 3: Implement scalar-generic Kerr–Schild fields**

In `kerr_schild_fields.h`, define a focused internal result carrying `r`,
`H`, `l_covariant`, `g_covariant`, and `g_contravariant`. Evaluate:

```cpp
q = x*x + y*y + z*z - a*a;
s = sqrt(q*q + 4*a*a*z*z);
r2 = q >= 0 ? 0.5*(q+s) : 2*a*a*z*z/(s-q);
r = sqrt(r2);
H = M*r/(r2 + a*a*(z/r)*(z/r));
l = {1,
     (r*x+a*y)/(r2+a*a),
     (r*y-a*x)/(r2+a*a),
     z/r};
g_cov[mu][nu] = eta[mu][nu] + 2*H*l[mu]*l[nu];
l_up = {-l[0], l[1], l[2], l[3]};
g_inv[mu][nu] = eta[mu][nu] - 2*H*l_up[mu]*l_up[nu];
```

Keep branch selection based on the scalar value for both `double` and
`Dual4`.

- [x] **Step 4: Implement public validation, accessors, and derivatives**

Validate constructor parameters once. `valid_point` checks finite
coordinates, stable `r > 64*epsilon*M`, nonzero finite denominators, and
finite metric values without excluding the axis, horizon, or positive-radius
interior. Metric methods throw `std::domain_error` when invalid. Compute
`a=chi*M` and `r_\pm=M±sqrt(M^2-a^2)`.

Implement the prompt's analytic radius gradient:

```text
D = 2r^2-rho^2+a^2
dr/dx = xr/D
dr/dy = yr/D
dr/dz = z(r^2+a^2)/(rD)
```

For inverse-metric derivatives, seed all four input coordinates as `Dual4`,
evaluate the scalar-generic inverse once, and copy each derivative component.
This keeps `KerrSchildCartesianMetric` fully executable at the end of Task 2;
no placeholder exception or zero derivative is permitted.

- [x] **Step 5: Run focused and baseline metric tests**

Run:

```bash
make -j4 tests/relativity/test_kerr_schild \
  tests/relativity/test_metrics \
  tests/relativity/test_metric_derivatives
./tests/relativity/test_kerr_schild
./tests/relativity/test_metrics
./tests/relativity/test_metric_derivatives
```

Expected: all tests pass and the reported maximum inverse errors satisfy the
declared gates.

- [x] **Step 6: Commit**

```bash
git add include/solar/relativity/kerr_schild_metric.h \
  src/relativity/kerr_schild_fields.h \
  src/relativity/kerr_schild_metric.cpp \
  tests/relativity/test_kerr_schild.cpp
git commit -m "feat(relativity): add Cartesian Kerr-Schild metric"
```

### Task 3: Independent Kerr–Schild derivative validation

**Files:**
- Create: `tests/relativity/test_kerr_schild_derivatives.cpp`
- Modify: `src/relativity/kerr_schild_metric.cpp`

**Interfaces:**
- Consumes: the production analytic gradient and scalar-generic AD inverse
  derivatives from Task 2.
- Produces: independent precision and boundary evidence; production changes
  occur only when a measured defect identifies its owning expression.

- [x] **Step 1: Write independent derivative tests**

At exterior, axis-adjacent, safe near-horizon, and interior fixtures for
positive and negative spin, independently seed `Dual4` coordinates and
compare its radius derivatives with:

```cpp
D = 2*r*r - (x*x + y*y + z*z) + a*a;
drdx = x*r/D;
drdy = y*r/D;
drdz = z*(r*r+a*a)/(r*D);
```

Require componentwise relative error `<1e-12`. For every
`partial_alpha g^{mu nu}`, compute the independent five-point stencil

```cpp
(-g(x+2h) + 8*g(x+h) - 8*g(x-h) + g(x-2h)) / (12*h)
```

with `h=epsilon^(1/5)*max(M,abs(x_alpha),1)` and require scaled relative
error `<3e-8`. Require time derivatives to be exactly zero and invalid
points to throw.

- [x] **Step 2: Run the independent derivative acceptance test**

Run:

```bash
make -j4 tests/relativity/test_kerr_schild_derivatives
```

Expected: the new independent test either passes every gate immediately or
reports the first measured component/fixture mismatch. A first-run pass is
valid because Task 2 already test-drove the complete concrete metric.

- [x] **Step 3: Correct only a measured derivative defect**

If a gate fails, use the printed coordinate, tensor indices, analytic value,
AD value, and finite-difference value to locate the defect in stable radius,
`H`, the null one-form, inverse assembly, or analytic `D`. Correct that
single owning expression. Preserve exact tensor symmetry by assigning paired
components from the same scalar expression; do not relax a gate or introduce
a finite-difference production fallback.

- [x] **Step 4: Run derivative, Hamiltonian, and geodesic baseline tests**

Run:

```bash
make -j4 tests/relativity/test_kerr_schild_derivatives \
  tests/relativity/test_hamiltonian \
  tests/relativity/test_geodesics
./tests/relativity/test_kerr_schild_derivatives
./tests/relativity/test_hamiltonian
./tests/relativity/test_geodesics
```

Expected: all pass with derivative maxima below the documented gates.

- [x] **Step 5: Commit**

```bash
git add src/relativity/kerr_schild_metric.cpp \
  tests/relativity/test_kerr_schild_derivatives.cpp
git commit -m "feat(relativity): differentiate Kerr-Schild inverse metric"
```

### Task 4: Complete BL↔KS canonical chart transform

**Files:**
- Create: `include/solar/relativity/kerr_chart_transform.h`
- Create: `src/relativity/kerr_chart_fields.h`
- Create: `src/relativity/kerr_chart_transform.cpp`
- Create: `src/relativity/kerr_chart_state_transform.cpp`
- Create: `tests/relativity/test_kerr_chart_transform.cpp`

**Interfaces:**
- Consumes: `KerrSchildCartesianMetric`, `Dual4`, `Mat4::inverse`, BL
  coordinate/state conventions.
- Produces: `KerrChartTransform` public methods declared in the design.

- [x] **Step 1: Write position, Jacobian, and canonical-state tests**

Use masses `{1,3}`, spins `{-0.8,0,0.7}`, and safe exterior points away from
the axis. Require wrapped position round trips and:

```cpp
const Mat4 forward =
    transform.boyer_lindquist_to_kerr_schild_jacobian(bl.x);
const Mat4 reverse =
    transform.kerr_schild_to_boyer_lindquist_jacobian(ks.x);
check(max_identity_error(multiply(reverse, forward)) < 5e-12);
```

Transform a nontrivial covariant momentum, then require every momentum
component to round-trip within `<1e-10`. Test the covector-vector pairing
with `V_KS=J V_BL` and `p_KS=(J^-1)^T p_BL`, and compare BL/KS Hamiltonians
within `<1e-10`. Include negative spin to fix the azimuth orientation.
Require explicit failure for invalid constructor inputs, `r<=r_++margin`,
the polar axis, non-finite inputs, and an unresolvable Jacobian.

- [x] **Step 2: Run the transform test and confirm the red state**

Run:

```bash
make -j4 tests/relativity/test_kerr_chart_transform
```

Expected: compilation fails because the transform header is absent.

- [x] **Step 3: Implement safe logarithmic coordinate functions**

For subextremal horizons `r_+` and `r_-`, use:

```text
F_t(r) =
  [2 M r_+/(r_+-r_-)] log|r-r_+|
  -[2 M r_-/(r_+-r_-)] log|r-r_-|
F_phi(r) =
  [a/(r_+-r_-)] log|(r-r_+)/(r-r_-)|
```

The safe-overlap check makes both arguments positive. Evaluate forward
position with `Dual4` so its full Jacobian includes `dt_KS/dr` and
`dphi_tilde/dr`. Use the design’s Cartesian formulas exactly.

The master prompt's earlier differential contract and its later displayed
`F_phi` sign conflict. Use the positive expression above because it alone
differentiates to `+a/Delta`, matches the documented ingoing Kerr transform,
and is protected by a direct finite-difference sign regression.

- [x] **Step 4: Implement inverse position and canonical transforms**

Recover `r` from `KerrSchildCartesianMetric`, then:

```text
theta = acos(z/r)
phi_tilde = atan2(r*y-a*x, r*x+a*y)
t_BL = t_KS-F_t(r)
phi_BL = wrap(phi_tilde-F_phi(r))
```

Compute the inverse Jacobian by inverting the forward Jacobian at the
recovered BL point. Apply:

```cpp
V_ks = J * V_bl;
p_ks = transpose(inverse(J)) * p_bl;
p_bl = transpose(J) * p_ks;
```

Reject unsafe or non-finite transformations with `std::domain_error`; never
copy or reinterpret `p_r`.

- [x] **Step 5: Run focused transform and BL metric tests**

Run:

```bash
make -j4 tests/relativity/test_kerr_chart_transform \
  tests/relativity/test_kerr_bl \
  tests/relativity/test_spacetime_algebra
./tests/relativity/test_kerr_chart_transform
./tests/relativity/test_kerr_bl
./tests/relativity/test_spacetime_algebra
```

Expected: all pass with printed round-trip and Jacobian maxima below gates.

- [x] **Step 6: Commit**

```bash
git add include/solar/relativity/kerr_chart_transform.h \
  src/relativity/kerr_chart_fields.h \
  src/relativity/kerr_chart_transform.cpp \
  src/relativity/kerr_chart_state_transform.cpp \
  tests/relativity/test_kerr_chart_transform.cpp
git commit -m "feat(relativity): add canonical Kerr chart transforms"
```

### Task 5: Standard Kerr–Schild horizon and interior events

**Files:**
- Create: `include/solar/relativity/kerr_schild_events.h`
- Create: `src/relativity/kerr_schild_events.cpp`
- Create: `tests/relativity/test_kerr_schild_events.cpp`

**Interfaces:**
- Consumes: `KerrSchildCartesianMetric`, `GeodesicEvent`.
- Produces: the three event helpers declared in the design.

- [x] **Step 1: Write event contract tests**

Require:

```cpp
check_near("default cutoff",
           kerr_schild_interior_cutoff_radius(metric, 0.0),
           0.05 * metric.mass(), 0.0);
check_near("configured floor",
           kerr_schild_interior_cutoff_radius(metric, 0.01),
           0.05 * metric.mass(), 0.0);
```

and a larger configured radius to be retained. Verify the horizon event is
named, decreasing, `HorizonCrossing`, and evaluates to zero at `r_+`;
verify the cutoff event is decreasing, `InteriorCutoff`, and zero at its
radius. Create an event from a temporary metric and invoke it later to prove
owned lifetime. Require negative/non-finite configured radii and
non-positive/non-finite root tolerances to throw.

- [x] **Step 2: Run the event test and confirm the red state**

Run:

```bash
make -j4 tests/relativity/test_kerr_schild_events
```

Expected: compilation fails because the event header is absent.

- [x] **Step 3: Implement validated event factories**

Capture `KerrSchildCartesianMetric` by value in each event callback and
return:

```cpp
metric.radial_coordinate(state.x) - target_radius;
```

Use `std::max(0.05*metric.mass(), configured_radius_M)`, treating zero as
the default. Event construction validates all scalar inputs before creating
the callback.

- [x] **Step 4: Run event tests**

Run:

```bash
make -j4 tests/relativity/test_kerr_schild_events \
  tests/relativity/test_geodesic_events
./tests/relativity/test_kerr_schild_events
./tests/relativity/test_geodesic_events
```

Expected: both pass.

- [x] **Step 5: Commit**

```bash
git add include/solar/relativity/kerr_schild_events.h \
  src/relativity/kerr_schild_events.cpp \
  tests/relativity/test_kerr_schild_events.cpp
git commit -m "feat(relativity): add Kerr-Schild boundary events"
```

### Task 6: BL/KS common-event crosscheck

**Files:**
- Create: `tests/relativity/test_kerr_bl_ks_crosscheck.cpp`

**Interfaces:**
- Consumes: BL and KS metrics, chart transform, local initialization,
  generic integration, and user events.
- Produces: independent end-to-end coordinate-equivalence evidence.

- [x] **Step 1: Add moderate null and timelike fixtures**

Build at least six fixtures spanning null/timelike kind, both spin signs,
prograde/retrograde momenta, and distinct inclinations. Initialize constrained
BL states, transform the complete state to KS, and integrate both charts to
the same decreasing physical radius outside the overlap margin. Transform
the KS event state back to BL. Record scaled errors for
`(t,r,theta,wrapped phi)` and all four covariant momentum components, plus
termination reason, Hamiltonian constraint, E, and Lz drift.

- [x] **Step 2: Run the new acceptance crosscheck**

Run:

```bash
make -j4 tests/relativity/test_kerr_bl_ks_crosscheck
./tests/relativity/test_kerr_bl_ks_crosscheck
```

Expected: the executable either passes every declared gate immediately or
reports the first measured physical mismatch. A first-run pass is valid
because Tasks 1-5 already provide the complete production surface; do not
manufacture a failure or alter correct production code merely to force red.

- [x] **Step 3: Correct only identified implementation defects**

Use the printed worst fixture/component to isolate whether the defect is in
position mapping, covector Jacobian orientation, KS metric derivatives,
event localization, or invariant callback wiring. Correct the owning L1
module without changing a gate or adding fixture-specific behavior.

- [x] **Step 4: Add a safe near-horizon overlap fixture**

Place the target above `r_+ + overlap_margin`, close enough to exercise the
large logarithmic derivatives. Require matching termination and position/
momentum errors `<1e-6`. For ordinary fixtures, sort errors and compute the
nearest-rank P95; require position and momentum P95 `<1e-8`, Hamiltonian
error `<1e-10`, and E/Lz drift `<1e-12`.

- [x] **Step 5: Run the complete crosscheck twice**

Run:

```bash
./tests/relativity/test_kerr_bl_ks_crosscheck
./tests/relativity/test_kerr_bl_ks_crosscheck
```

Expected: both runs pass and report identical fixture classifications and
gate values to printed precision.

- [x] **Step 6: Commit**

```bash
git add include/solar/relativity src/relativity \
  tests/relativity/test_kerr_bl_ks_crosscheck.cpp
git commit -m "test(relativity): crosscheck BL and Kerr-Schild flows"
```

### Task 7: Timelike horizon crossing and interior continuation

**Files:**
- Create: `tests/relativity/test_geodesics_kerr_schild.cpp`

**Interfaces:**
- Consumes: `initialize_local_timelike`, `KerrChartTransform`,
  `GeodesicIntegrator`, standard KS events, `make_arbitrary_observer`, and
  tetrad validation.
- Produces: the Phase 4 horizon-to-interior acceptance proof.

- [x] **Step 1: Write the two-segment plunge**

Initialize a unit-mass ingoing timelike BL state from a safe exterior local
observer, transform it to KS, and configure:

```cpp
config.kind = GeodesicKind::TimelikeUnitMass;
config.monitor_energy = true;
config.monitor_lz = true;
config.stationary_energy_evaluator =
    kerr_schild_stationary_energy;
config.axial_angular_momentum_evaluator =
    kerr_schild_axial_angular_momentum;
```

First integrate with only the outer-horizon event. Require
`HorizonCrossing`, `abs(r-r_+)` within the event tolerance, finite position,
momentum, metric, inverse, affine, and Hamiltonian norm, constraint
`<1e-10`, and E/Lz drift `<1e-12`.

- [x] **Step 2: Run the new horizon acceptance test**

Run:

```bash
make -j4 tests/relativity/test_geodesics_kerr_schild
./tests/relativity/test_geodesics_kerr_schild
```

Expected: the executable either passes immediately or identifies the first
measured horizon, constraint, invariant, observer, or continuation failure.
This composes already test-driven units, so a first-run pass is acceptable.

- [x] **Step 3: Validate a freely falling horizon observer**

Raise the event momentum with the KS inverse metric, normalize it as the
observer four-velocity, and call `make_arbitrary_observer` with three
linearly independent spatial seeds. Require finite tetrad components and
maximum orthonormality error `<1e-12`.

- [x] **Step 4: Restart the exact event state to the interior cutoff**

Call `integrate` again using the exact first segment final state and only
`make_kerr_schild_interior_cutoff_event(metric, 0.0, tolerance)`. Require
`InteriorCutoff`, `r=0.05M` within tolerance, finite state, strictly larger
affine/proper time, and continuity of the initial second-segment state with
the first endpoint. Do not nudge or clamp the horizon state.

- [x] **Step 5: Correct owning code and run focused horizon tests**

If a gate fails, trace it to metric values, derivatives, callbacks, event
localization, or observer construction and correct only that owner. Run:

```bash
./tests/relativity/test_geodesics_kerr_schild
./tests/relativity/test_kerr_schild
./tests/relativity/test_kerr_schild_derivatives
./tests/relativity/test_kerr_schild_events
```

Expected: all pass with the measured horizon/interior gates printed.

- [x] **Step 6: Commit**

```bash
git add include/solar/relativity src/relativity \
  tests/relativity/test_geodesics_kerr_schild.cpp
git commit -m "test(relativity): validate Kerr horizon continuation"
```

### Task 8: Installed consumer and public surface

**Files:**
- Modify: `tests/external_consumer/probe.cpp`

**Interfaces:**
- Consumes: installed `Solar::Relativity` package.
- Produces: compile/link/runtime proof for public KS metric and short
  Hamiltonian integration.

- [x] **Step 1: Extend the external probe**

Include the public metric and event headers, instantiate
`KerrSchildCartesianMetric(1.0, 0.5)`, evaluate an exterior metric/inverse,
and run a short finite KS geodesic segment. Return a distinct nonzero code if
the chart, inverse identity, or integration result is invalid.

- [x] **Step 2: Verify the installed consumer**

Run:

```bash
make test-external-consumer
```

Expected: installation, standalone CMake configure/build, and probe execution
all pass without source-tree includes.

- [x] **Step 3: Verify no build-list edit is required**

Run:

```bash
cmake -S . -B build-phase4-release \
  -DCMAKE_BUILD_TYPE=Release -DSOLAR_BUILD_CLI=ON
cmake --build build-phase4-release -j4
```

Expected: CMake `CONFIGURE_DEPENDS` discovers every new source. Do not edit
`CMakeLists.txt` or `Makefile` unless this command demonstrates a discovery
or export defect.

- [x] **Step 4: Commit**

```bash
git add tests/external_consumer/probe.cpp
git commit -m "test(relativity): consume Kerr-Schild public API"
```

### Task 9: Validation, documentation, and release-candidate audit

**Files:**
- Create: `docs/validation/relativity_09_kerr_schild.md`
- Modify: `RELATIVITY_STATUS.md`
- Modify only if evidence requires it: Phase 4 implementation/tests above.

**Interfaces:**
- Consumes: all Phase 4 tests and measured outputs.
- Produces: reproducible validation record and a clean release candidate.

- [ ] **Step 1: Run clean Release verification**

Run:

```bash
make clean
make -j4
make -j4 test
make test-external-consumer
git diff --check
```

Record the exact passed/failed totals and all optional skips. A required
failure keeps Phase 4 incomplete.

- [ ] **Step 2: Run focused sanitizer verification**

Build the seven focused executables once with the combined sanitizers:

```bash
make clean
make CXXFLAGS='-std=c++17 -O1 -g -Wall -Wextra -Iinclude \
  -fsanitize=address,undefined -fno-omit-frame-pointer' \
  tests/relativity/test_geodesic_invariant_callbacks \
  tests/relativity/test_kerr_schild \
  tests/relativity/test_kerr_schild_derivatives \
  tests/relativity/test_kerr_chart_transform \
  tests/relativity/test_kerr_schild_events \
  tests/relativity/test_kerr_bl_ks_crosscheck \
  tests/relativity/test_geodesics_kerr_schild
```

Run each executable with:

```bash
ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 \
UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
  ./tests/relativity/test_kerr_schild
```

Repeat the same environment for the other six focused executables. Expected:
every assertion passes with no sanitizer diagnostic. Leak detection stays
disabled because the existing macOS sanitizer validation contract excludes
platform runtime allocations.

- [ ] **Step 3: Perform a targeted mutation check**

Temporarily transpose the KS covector transform in the wrong direction and
run `test_kerr_chart_transform` plus `test_kerr_bl_ks_crosscheck`; both must
fail. Restore the correct implementation using `apply_patch`, rerun both,
and require green. Do not commit the mutation.

- [ ] **Step 4: Self-audit architecture and test meaning**

Check that:

- no L0 type imports KS-specific code;
- no generic integrator imports the KS metric;
- callbacks are disabled by default and evaluate only when enabled;
- every public transform applies the complete Jacobian;
- the horizon is not rejected by `valid_point`;
- the interior cutoff is an event, not a metric-domain rule;
- tests would fail for a sign flip in `a`, `H`, `l_t`, or covector transform;
- new files over the 200-line review signal still own one cohesive
  responsibility or are split at a natural boundary;
- no dead code, unrelated edits, hidden fallback, or stale status item
  remains.

- [ ] **Step 5: Write measured validation documentation**

Create `docs/validation/relativity_09_kerr_schild.md` with formulas and
conventions, fixture families, exact commands, compiler/build type, measured
max/P95 values, sanitizer totals, installed-consumer result, explicit
non-goals, and any optional external-data skip. Update the existing
`RELATIVITY_STATUS.md` fields and lists (`CURRENT_PHASE`, `PHASE_STATE`,
verified commands, completed work, missing work, blockers, likely bugs,
falsification commands, and `NEXT_ALLOWED_ACTION`); move Phase 4 out of
current work only after every required gate is green.

- [ ] **Step 6: Run the final release-candidate command set**

Run:

```bash
make clean
make -j4 test
make test-external-consumer
git diff --check
git status --short
```

Expected: all required tests pass, sanitizer output is clean, diff checking
is clean, normal Release artifacts have been restored after sanitizer
testing, and only intentional Phase 4 files are modified.

- [ ] **Step 7: Commit**

```bash
git add docs/validation/relativity_09_kerr_schild.md \
  RELATIVITY_STATUS.md
git commit -m "docs(relativity): validate Kerr-Schild Phase 4"
```

- [ ] **Step 8: Publish one public pull request and validate Linux CI**

Push `codex/relativity-phase-4`, open a public PR against `main`, and let the
single intentional release-candidate CI run. If Linux/GCC exposes a
platform defect, reproduce or isolate it locally, add a behavior-focused
regression, and fix the cause without weakening the physics gate. Merge only
after CI passes and a final PR diff review finds no unrelated change.
