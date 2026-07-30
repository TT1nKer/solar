# Solar Relativity Phase 3 Kerr Separated/Mino Design

**Status:** Approved for implementation on 2026-07-30
**Source contract:** `SOLAR_RELATIVITY_KERR_完整实施主提示词_v3.md`
**Previous gate:** Phase 2 `PASSED`; Phase 3 starts from public `main` at
`699f6283e4874f14db79981e22cf64dd328ae762`

## Goal

Add a reusable, allocation-bounded Kerr Boyer-Lindquist fast path for null and
unit-mass timelike geodesics:

- separated radial and polar potentials in Mino time;
- explicit radial and polar turning-point handling;
- render-relevant diagnostics such as minimum radius and winding;
- the existing terminal geodesic event contract;
- worldline cross-checks against the authoritative Hamiltonian BL solver.

This is a Solar library capability. Gargantua Studio will consume only the
public installed API after the Phase 3 validation gate passes.

Phase 3 does not add Kerr-Schild coordinates, physical horizon crossing,
matter, radiative transfer, a renderer, CUDA, C ABI, WASM, or UI.

## Sources and authority

The implementation order of authority is:

1. the v3 project contract and Solar's documented sign/unit conventions;
2. the existing generic Hamiltonian BL integrator for numerical worldline
   cross-checks;
3. published separated Kerr equations and independent analytic fixtures.

Mino time decouples the radial and polar equations as described by
[Fujita and Hikida](https://arxiv.org/abs/0906.1420). Public fast Kerr ray
codes such as [GEOKERR](https://arxiv.org/abs/0903.0620) demonstrate the value
of separated equations, but Phase 3 does not copy their implementation or
introduce an analytic elliptic-integral backend. A future optimized backend
may use Carlson elliptic integrals behind the same validated public contract.

## Approaches considered

### Permanent second-order potential evolution

Evolving `r''=R'/2` and `mu''=U'/2` crosses simple turning points smoothly,
but numerical drift can move the solution away from the first-integral
surface. It also obscures the v3 requirement to reject negative-potential
trials, locate roots, and distinguish simple from double roots.

### Full analytic elliptic-integral solver

This can be faster and highly accurate, but its branch structure, root
classification, and licensing/provenance review are a much larger first
change. It is not the smallest Phase 3 implementation that unlocks a CPU
high-throughput reference path.

### Signed first-order equations with explicit turning roots

This is the selected approach. DOPRI5 evolves the signed first-order
separated equations. A trial that enters a forbidden potential region is
rejected. The solver brackets and locates the root, advances to it, flips the
appropriate direction only for a simple root, and performs a bounded local
release step. A double or numerically unresolved root terminates explicitly
as `NearCriticalOrbit`.

## Conventions and equations

- Signature: `(-,+,+,+)`.
- Coordinates: `(t,r,theta,phi)`.
- Canonical momentum: covariant `p_mu`.
- Geometrized units: `G=c=1`; `mass_M` is a length.
- `a=spin_chi*mass_M`.
- `mu_p^2` is `0` for null and `1` for unit-mass timelike motion.
- Mino parameter is `gamma`; affine/proper parameter is `lambda`.
- `d lambda = Sigma d gamma`.
- Backward tracing uses a negative Mino step and keeps physical momentum
  future-directed.

The internal fixed-size state is:

```text
y = [t, r, mu, phi, lambda],  mu = cos(theta)
```

`phi` remains unwrapped so winding is not lost.

```text
Delta = r^2 - 2 M r + a^2
Sigma = r^2 + a^2 mu^2
P     = E (r^2 + a^2) - a Lz
K     = mu_p^2 r^2 + (Lz-aE)^2 + Q
R     = P^2 - Delta K
A     = a^2 (mu_p^2-E^2)
U     = Q(1-mu^2) - mu^2[A(1-mu^2)+Lz^2]
```

```text
dr/dgamma      = s_r  sqrt(R)
dmu/dgamma     = s_mu sqrt(U)
dphi/dgamma    = a P/Delta + Lz/(1-mu^2) - a E
dt/dgamma      = (r^2+a^2) P/Delta
                 + a[Lz-aE(1-mu^2)]
dlambda/dgamma = Sigma
```

Analytic derivatives used only for root classification and local release are:

```text
R' = 4 E r P - Delta' K - Delta (2 mu_p^2 r)
Delta' = 2(r-M)
U' = -2 mu [Q + A + Lz^2 - 2 A mu^2]
```

The production RHS never uses `sqrt(max(V,0))` to hide a forbidden state.
A small negative value is admissible only inside root-location arithmetic,
never as a continuing trajectory state.

## Public API

`include/solar/relativity/kerr_separated.h` owns the stable L3-facing
contract:

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
    double min_mino_step;
    double max_mino_step;
    double min_radius_M;
    double azimuthal_advance;
    double winding;
    double max_radial_residual;
    double max_polar_residual;
    double max_constraint_error;
    double max_carter_rel_error;
    TerminationReason reason;
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
};
```

The exact default member values are implementation details except where the
existing aggregate compatibility tests require stable defaults. Invalid
enums, non-finite limits/tolerances, inconsistent step signs, invalid metric
points, and malformed events fail before trajectory work.

`potential_tolerance` is a dimensionless normalized first-integral tolerance.
`root_tolerance` bounds `abs(delta r)/M` or `abs(delta mu)`.
`critical_derivative_tolerance` is applied after scaling radial `R'` by
`M^3` and polar `U'` by `M^2`. `polar_axis_tolerance` is a dimensionless
lower bound on `1-mu^2`.

`TerminationReason` gains `NearCriticalOrbit`. Existing enum values keep
their order; the new value is appended to preserve source behavior for
ordinary callers.

## Internal boundaries

### L1 potential mathematics

`src/relativity/kerr_separated_potentials.{h,cpp}` owns `Delta`, `Sigma`,
`P`, `R`, `U`, analytic derivatives, normalized residuals, and finite-domain
validation. It has no integration loop or event knowledge.

### L1 state conversion

`src/relativity/kerr_separated_state.{h,cpp}` converts a valid Hamiltonian
state to constants, internal Mino state, and initial direction classification,
and reconstructs public `PhaseSpaceState`:

```text
p_t     = -E
p_phi   = Lz
p_r     = (dr/dgamma)/Delta
p_theta = -(dmu/dgamma)/sin(theta)
```

The signs are taken from the Hamiltonian tangent multiplied by `Sigma`, not
from canonical momentum signs. At a zero tangent, the potential derivative
and Hamiltonian trend determine a simple-root departure. Exact equatorial
`mu=0`, `p_theta=0`, `U=U'=0` is a locked symmetry plane, not an error.

Near `|mu|=1` with nonzero `Lz`, the BL azimuthal term is explicitly rejected
and deferred to Phase 4 Kerr-Schild work.

### L1 turning operation

`src/relativity/kerr_separated_turning.{h,cpp}` owns bracketed scalar
potential roots, simple/double-root classification, and the bounded release
state. It does not own the full trajectory loop.

### L2 trajectory flow

`src/relativity/kerr_separated.cpp` owns validation, adaptive stepping,
event selection, limits, diagnostics, and conversion to the public result.
It composes only the L1 blocks and existing L0/public types. No reverse
dependency from the generic Hamiltonian solver is introduced.

## Turning-point algorithm

For either `R(r)` or `U(mu)`:

1. require the accepted state to have nonnegative potential within the
   configured numerical tolerance;
2. attempt a signed DOPRI5 step;
3. if a stage or endpoint enters a forbidden potential region, reject the
   trial and shrink the step;
4. once positive and negative evaluations bracket a root, locate it with a
   safeguarded secant/bisection solve in the affected coordinate;
5. advance the complete state to the root using dense output, subject to the
   affine root tolerance converted by `d lambda/d gamma=Sigma`;
6. evaluate the analytic potential derivative;
7. if both potential and derivative meet the critical tolerances, terminate
   as `NearCriticalOrbit`;
8. otherwise flip only the affected direction, increment its turn count, and
   use `q''=V'/2` for one bounded local release displacement;
9. require the released point to lie in the allowed potential region before
   resuming first-order DOPRI5.

The local second-order relation is a root-release operation, not the
long-term evolution equation. Step underflow or an unresolvable bracket is
reported explicitly; the solver never sticks at a root or fabricates a
direction.

## Events and limits

The existing `GeodesicEvent` and `EventHit` types are reused. Event functions
receive reconstructed public Hamiltonian states. Terminal event ordering and
direction semantics remain identical to the generic solver.

An event's `root_tolerance` is an affine-length tolerance. The separated flow
converts it locally to a Mino interval using finite positive `Sigma`; it does
not reinterpret public event units.

Initial events are checked before any accepted or rejected step. The first
event in integration direction wins. Limits include affine, coordinate time,
total steps, per-step rejections, and representable Mino progress.

No full trajectory vector is returned. Cross-checks use terminal events at
common crossings, keeping the hot path allocation-bounded and avoiding a
premature sampling API.

## Diagnostics

The solver records:

- accepted and rejected trials;
- minimum and maximum accepted absolute Mino step;
- radial and polar turn counts;
- minimum Boyer-Lindquist radius;
- unwrapped azimuthal advance and `winding=advance/(2 pi)`;
- maximum normalized radial and polar first-integral residual;
- maximum normalized Hamiltonian constraint error and Carter relative drift
  after reconstructing public states;
- explicit termination reason and message.

Residual scales are dimensionally meaningful:

```text
radial scale = max(M^4, |(dr/dgamma)^2|, |R|)
polar scale  = max(M^2, |(dmu/dgamma)^2|, |U|)
```

Diagnostics use quiet NaN for values that were not evaluated. Non-finite
diagnostics are not silently replaced with zero.

## Validation

Behavior-focused tests are split by responsibility:

- `test_kerr_separated_potentials.cpp`: literal formulas, units, analytic
  derivatives versus finite differences, spin sign, null/timelike inputs,
  and invalid domains;
- `test_kerr_separated.cpp`: construction, reversal, equatorial lock,
  radial and polar simple turns, critical double roots, axis rejection,
  limits, malformed inputs, event ordering, and diagnostics;
- `test_kerr_separated_crosscheck.cpp`: common-event worldline comparisons
  against the generic Hamiltonian BL solver.

Required cross-check families are:

1. ordinary null and timelike motion at a common radial or coordinate-time
   event;
2. a scattering null ray with one radial turn and a return crossing;
3. an off-equatorial ray with a polar turn;
4. Schwarzschild and positive/negative Kerr spin cases;
5. a spherical critical photon orbit that returns `NearCriticalOrbit`;
6. convergence under a halved maximum step.

Comparisons are made at common physical events, never at equal Mino and
affine parameter values. Ordinary exterior position disagreement must meet
the v3 `p95 < 1e-8 M`, `max < 1e-7 M` target; the stricter test fixture
tolerances are recorded in the validation report. Accepted rays must also
retain the existing `1e-10` Hamiltonian/Carter reference gates when their
reconstructed states are evaluated.

The final gate runs focused Phase 3 tests, the full Solar suite, the external
installed-package consumer, sanitizer builds for the new tests, and
`git diff --check`. `RELATIVITY_STATUS.md` advances to Phase 3 only after all
required evidence passes.

## Compatibility and remaining risk

- Existing generic geodesic APIs and behavior remain unchanged except for the
  appended termination reason.
- The separated solver is BL-exterior-only. Invalid BL domain is not physical
  horizon capture.
- Simple turning roots are supported; double/near-critical roots are
  classified, not integrated for arbitrarily long dwell time.
- Near-polar nonzero-`Lz` motion remains unsupported until Phase 4.
- The first backend is numerical DOPRI5, not the final analytic or CUDA
  throughput ceiling.
- Public configuration and result types are source-level pre-1.0 contracts;
  serialized scene or render formats are not introduced in Solar Phase 3.
