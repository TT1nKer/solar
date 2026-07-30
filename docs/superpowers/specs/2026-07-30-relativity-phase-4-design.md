# Solar Relativity Phase 4: Cartesian Kerr–Schild Design

## Goal

Add an ingoing Cartesian Kerr–Schild metric and complete Boyer–Lindquist /
Kerr–Schild state transforms so Solar can evolve null and unit-mass timelike
Hamiltonian geodesics through the outer Kerr horizon to an explicit finite
interior cutoff.

Phase 4 is a coordinate and validation milestone. It does not add radiative
transfer, a long-time structure-preserving orbit integrator, a renderer, GPU
code, or singularity physics.

## Authority and conventions

The implementation follows the formulas and phase gate in
`SOLAR_RELATIVITY_KERR_完整实施主提示词_v3.md`.

Independent references:

- Chan, Medeiros, Özel, and Psaltis, “GRay2: A General Purpose Geodesic
  Integrator for Kerr Spacetimes,” <https://arxiv.org/abs/1706.07062>.
- SpECTRE `gr::KerrSchildCoords` documentation,
  <https://spectre-code.org/classgr_1_1KerrSchildCoords.html>.
- Einstein Toolkit `Exact` thorn Kerr–Schild documentation,
  <https://einsteintoolkit.org/thornguide/EinsteinInitialData/Exact/documentation.html>.

No external implementation code is copied.

Solar keeps:

- signature `(-,+,+,+)`;
- coordinate order `(t,x,y,z)` in Cartesian Kerr–Schild;
- canonical state order `(x^mu,p_mu)`;
- covariant canonical momentum;
- geometrized units `G=c=1`;
- subextremal `|chi|<1` first-version scope.

## Approaches considered

### Selected: analytic metric values plus forward AD derivatives

Evaluate the implicit radius, Kerr–Schild scalar, null one-form, covariant
metric, and analytic inverse directly. Instantiate the inverse expression
with `Dual4` to obtain all four coordinate derivatives. Independently compare
the AD radius derivative with the prompt’s analytic `dr/dx^i`, then compare
the inverse-metric derivatives with a five-point finite difference.

This follows the master prompt, preserves the generic Hamiltonian integrator,
and minimizes hand-derived tensor code.

### Rejected for Phase 4: fully handwritten metric derivatives

This can reduce arithmetic, but duplicates a large sign-sensitive derivative
surface before a measured performance need exists. It is a later optimization
only if profiles show the AD path dominates and exact equivalence tests exist.

### Rejected: finite-difference production derivatives

This is small to implement but makes Hamiltonian evolution depend on a
coordinate step heuristic precisely where the horizon and interior require
the most reliable derivatives. Five-point differences remain validation
evidence, not production physics.

## Layer and module boundaries

The dependency direction remains `L3 -> L2 -> L1 -> L0`.

### L0: existing mathematical and public state contracts

No coordinate-specific behavior enters `types.h`, `math.h`, or `dual4.h`.
The existing `Mat4`, `Contravariant4`, `Covariant4`, `PhaseSpaceState`, and
`Dual4` operations are sufficient. The existing `Chart::KerrSchildCartesian`
enumerator becomes executable.

### L1: Cartesian Kerr–Schild metric

New public header:

`include/solar/relativity/kerr_schild_metric.h`

Public class:

```cpp
class KerrSchildCartesianMetric final : public Metric {
public:
    KerrSchildCartesianMetric(double mass_M, double spin_chi);

    Chart chart() const noexcept override;
    std::string name() const override;
    Mat4 covariant(const Contravariant4& x) const override;
    Mat4 contravariant(const Contravariant4& x) const override;
    std::array<Mat4, 4>
    contravariant_derivatives(const Contravariant4& x) const override;
    bool valid_point(const Contravariant4& x) const noexcept override;

    double radial_coordinate(const Contravariant4& x) const;
    Vec3 radial_coordinate_gradient(const Contravariant4& x) const;
    double mass() const noexcept;
    double spin_chi() const noexcept;
    double spin_length() const noexcept;
    double outer_horizon_radius() const noexcept;
    double inner_horizon_radius() const noexcept;
};

double kerr_schild_stationary_energy(const PhaseSpaceState& state);
double kerr_schild_axial_angular_momentum(
    const PhaseSpaceState& state);
```

Implementation ownership:

- `src/relativity/kerr_schild_fields.h` owns the scalar-generic stable radius,
  `H`, null one-form, and metric expressions used by both `double` and
  `Dual4`.
- `src/relativity/kerr_schild_metric.cpp` owns validation, the public metric,
  the analytic radius gradient, and Cartesian invariant helpers.

The stable radius evaluation uses

```text
q = x^2 + y^2 + z^2 - a^2
s = sqrt(q^2 + 4 a^2 z^2)
r^2 = (q + s) / 2                         when q >= 0
r^2 = 2 a^2 z^2 / (s - q)                when q < 0
```

The second branch avoids cancellation inside the oblate disk. A point is
invalid when an input is non-finite, the derived radius is not finite and
strictly above `64 * epsilon * M`, a required denominator is zero/non-finite,
or a metric component is non-finite. The polar axis and both horizons are
valid. The ring/zero-radius branch is not.

Use the stable equivalent

```text
H = M r / [r^2 + a^2 (z/r)^2]
```

instead of forming fourth powers.

The analytic inverse is

```text
g^mu,nu = eta^mu,nu - 2 H l^mu l^nu
l^mu = eta^mu,nu l_nu
```

and is not obtained by numerical matrix inversion.

### L1: safe chart transform

New public header:

`include/solar/relativity/kerr_chart_transform.h`

Public class:

```cpp
class KerrChartTransform {
public:
    KerrChartTransform(
        double mass_M,
        double spin_chi,
        double overlap_margin_fraction = 1.0e-4);

    Contravariant4 position_to_kerr_schild(
        const Contravariant4& boyer_lindquist) const;
    Contravariant4 position_to_boyer_lindquist(
        const Contravariant4& kerr_schild) const;

    Mat4 boyer_lindquist_to_kerr_schild_jacobian(
        const Contravariant4& boyer_lindquist) const;
    Mat4 kerr_schild_to_boyer_lindquist_jacobian(
        const Contravariant4& kerr_schild) const;

    PhaseSpaceState state_to_kerr_schild(
        const PhaseSpaceState& boyer_lindquist) const;
    PhaseSpaceState state_to_boyer_lindquist(
        const PhaseSpaceState& kerr_schild) const;
};
```

Implementation ownership is split at the natural review boundary:

- `src/relativity/kerr_chart_fields.h` owns scalar-generic logarithmic
  offsets and the forward position expression;
- `src/relativity/kerr_chart_transform.cpp` owns parameter/domain checks,
  position maps, and Jacobians;
- `src/relativity/kerr_chart_state_transform.cpp` owns canonical covector
  and affine-preserving state transforms.

The transform uses the prompt’s zero-additive-constant convention:

```text
t_KS = t_BL + F_t(r)
phi_tilde = phi_BL + F_phi(r)
x = (r cos(phi_tilde) - a sin(phi_tilde)) sin(theta)
y = (r sin(phi_tilde) + a cos(phi_tilde)) sin(theta)
z = r cos(theta)
```

The inverse position uses:

```text
r = r(x,y,z)
theta = acos(z/r)
phi_tilde = atan2(r y - a x, r x + a y)
t_BL = t_KS - F_t(r)
phi_BL = phi_tilde - F_phi(r)
```

Azimuth comparisons wrap modulo `2*pi`.

The prompt's differential contract requires
`d(phi_tilde)=d(phi_BL)+a/Delta dr`, while its later displayed antiderivative
has a contradictory leading minus sign. Differentiating that displayed minus
would produce `-a/Delta`. Phase 4 follows the differential contract and the
same relation documented by SpECTRE:
`d(phi_BL)=d(phi_tilde)-a/Delta dr`. Therefore

```text
F_phi(r) =
  [a/(r_+-r_-)] log((r-r_+)/(r-r_-))
```

in the safe exterior overlap. A finite-difference regression checks
`dF_phi/dr=+a/Delta` directly.

The forward Jacobian is obtained from the full `Dual4` position expression.
For `J^alpha'_mu = partial x^alpha' / partial x^mu`:

```text
V^alpha' = J^alpha'_mu V^mu
p_alpha' = (J^-1)^mu_alpha' p_mu
```

The inverse-state path recovers the BL point, recomputes the forward Jacobian,
and applies `p_mu = J^alpha'_mu p_alpha'`. It never rotates a three-velocity
or copies `p_r` into a Cartesian component.

Both transform directions reject:

- non-finite inputs;
- `r <= r_+ + overlap_margin`;
- an unresolvable polar axis;
- a non-finite logarithm/Jacobian;
- a singular or ill-conditioned Jacobian.

This restriction applies only to the transform. It does not restrict the
Kerr–Schild metric or integration domain.

### L1: chart-correct invariant monitoring

Append two optional callbacks to `GeodesicIntegrationConfig`:

```cpp
InvariantEvaluator stationary_energy_evaluator;
InvariantEvaluator axial_angular_momentum_evaluator;
```

The existing booleans still decide whether each diagnostic is enabled.
When a callback is empty, the existing BL-compatible defaults remain
`E=-p_0` and `Lz=p_3`. When present, the monitor evaluates the callback at
the initial and every accepted/event state with the same explicit exception
and non-finite failure semantics as the Carter callback.

For Cartesian Kerr–Schild:

```text
E = -p_t
Lz = x p_y - y p_x
```

The generic monitor depends only on callbacks and state, never on the
Kerr–Schild class. Existing aggregate initialization and existing callers
remain source-compatible because the new fields are appended with empty
defaults.

### L1: standard horizon and interior events

New public header:

`include/solar/relativity/kerr_schild_events.h`

Public functions:

```cpp
double kerr_schild_interior_cutoff_radius(
    const KerrSchildCartesianMetric& metric,
    double configured_radius_M);

GeodesicEvent make_kerr_schild_horizon_event(
    const KerrSchildCartesianMetric& metric,
    double root_tolerance);

GeodesicEvent make_kerr_schild_interior_cutoff_event(
    const KerrSchildCartesianMetric& metric,
    double configured_radius_M,
    double root_tolerance);
```

The interior radius is
`max(0.05 * metric.mass(), configured_radius_M)`. A configured radius of zero
selects the `0.05M` default. Negative/non-finite configured radii and
non-positive/non-finite root tolerances are rejected. Each event owns a copy
of the immutable metric value, so its callback cannot dangle.

The horizon event is decreasing and reports `HorizonCrossing`. The interior
event is decreasing and reports `InteriorCutoff`. Both evaluate the same
implicit radial coordinate as the metric; neither classifies metric-domain
failure as an event.

### L2: generic Hamiltonian flow

No second Kerr–Schild integrator is introduced. `GeodesicIntegrator` already
consumes `Metric::contravariant_derivatives`, performs adaptive DOPRI5,
localizes public events, maintains unit-mass proper time through the affine
parameter, and rejects non-finite trial states.

Phase 4 supplies:

- a KS metric;
- chart-correct E/Lz evaluators;
- standard horizon/interior radial events;
- transformed safe-exterior initial states.

## Horizon and interior flow

The acceptance flow uses a unit-mass timelike plunge initialized with an
existing BL observer at a safe exterior matching radius. The complete
canonical state is transformed once to Cartesian KS.

1. Integrate to the decreasing event
   `r(x,y,z) - r_+ = 0` with reason `HorizonCrossing`.
2. Verify the event state, metric, inverse, momentum, Hamiltonian norm,
   affine/proper time, E, and Lz are finite.
3. Build an arbitrary freely falling observer from the horizon state’s
   raised timelike tangent and require the tetrad gate `<1e-12`.
4. Restart from the exact finite horizon event state without the horizon
   event.
5. Integrate to decreasing
   `r(x,y,z) - 0.05M = 0` with reason `InteriorCutoff`.
6. Require increasing continuous affine/proper time and finite state through
   the interior endpoint.

The two segments make the terminal event API explicit while proving that the
horizon state itself is a valid continuation state. No nudge, clamping, or
fabricated horizon crossing is allowed.

The interior event is a configured model boundary, not part of
`Metric::valid_point()` and not a singularity claim.

## BL / KS overlap flow

Use several moderate null/timelike fixtures and two common physical radius
events:

- ordinary exterior comparison;
- safe near-horizon exterior comparison.

For each fixture:

1. initialize a constrained BL state;
2. transform the full initial state to KS;
3. integrate the BL state with a BL radius event;
4. integrate the KS state with the same physical implicit-radius event;
5. transform the KS event state back to BL;
6. compare `(t,r,theta,wrapped phi)` and all four covariant momentum
   components;
7. compare event classifications and invariant gates.

Comparisons use common events, not equal integrator step counts or equal
coordinate parameter samples.

## Error handling

- Constructors reject non-positive/non-finite mass and non-finite or
  extremal/superextremal spin.
- Metric methods throw `std::domain_error` outside `valid_point()`.
- Transform methods throw `std::domain_error` outside the safe overlap or
  when a Jacobian cannot be inverted.
- Invariant callbacks that throw or return non-finite values terminate the
  geodesic as `NonFiniteState`; no diagnostic is silently dropped.
- BL transform failure is never classified as capture or a physical
  divergence.
- KS geodesic invalidity remains `InvalidMetricPoint` unless an explicit
  configured event was localized first.

## Validation gates

### Metric

- symmetry: exact component equality;
- Minkowski-null one-form:
  `|eta^mu,nu l_mu l_nu| < 5e-13`;
- ordinary inverse identity:
  `max|g_mu,alpha g^alpha,nu-delta_mu^nu| < 5e-13`;
- safe near-horizon inverse identity: `<1e-10`;
- KS AD inverse derivative versus five-point difference: relative `<3e-8`;
- analytic radius gradient versus independent `Dual4`: relative `<1e-12`;
- `a=0` agrees with the analytic Cartesian Schwarzschild KS expression;
- axis, outer horizon, and a point just inside the horizon remain valid;
- ring/zero-radius and non-finite points are rejected explicitly.

### Transform

- position round trip at multiple spins/masses;
- Jacobian products differ from identity by `<5e-12`;
- full canonical state round trip position/momentum error `<1e-10`;
- Hamiltonian and covector-vector pairing are invariant to `<1e-10`;
- negative-spin orientation is covered;
- unsafe horizon margin, axis, non-finite, and extremal inputs fail
  explicitly.

### Geodesic

- ordinary BL/KS common-event position P95 `<1e-8 M`;
- ordinary momentum P95 `<1e-8`;
- safe near-horizon overlap position and momentum `<1e-6`;
- termination classification agrees;
- ordinary normalized Hamiltonian error `<1e-10`;
- stationary energy relative drift `<1e-12`;
- Cartesian axial angular momentum relative drift `<1e-12`;
- timelike horizon and interior states are finite;
- affine/proper time is continuous and increases after crossing;
- horizon freely falling tetrad orthonormality `<1e-12`;
- interior termination is `InteriorCutoff` at `0.05M`.

### Repository

- focused Release tests pass;
- the complete existing Release suite passes;
- installed external CMake consumer still passes;
- new focused tests pass under AddressSanitizer and UndefinedBehaviorSanitizer;
- `git diff --check` is clean;
- Linux/GCC CI must pass before merge.

## Compatibility and migration

- Existing public signatures remain valid.
- Existing BL energy/Lz monitoring behavior is unchanged when evaluators are
  absent.
- CMake and Make source/test globs should discover the new `.cpp` files and
  tests; build files change only if a verified export/build defect exists.
- `Chart::KerrSchildCartesian` changes from representable-only to executable.
- The physics contract string remains unchanged unless an existing versioning
  rule explicitly requires a bump; Phase number and contract string are not
  assumed to be identical concepts.
- Gargantua remains pinned to the merged Phase 3 commit until Phase 4 is
  merged and independently consumed.

## Explicit non-goals

- no negative-`r` analytic extension;
- no extremal `|chi|=1` transform;
- no singularity or quantum-gravity model;
- no long-time symplectic/structure-preserving timelike solver;
- no Carter Killing-tensor evaluator in Cartesian KS;
- no selected-event history API;
- no radiative transfer, material model, renderer, WASM, or GPU work;
- no replacement of the separated solver for ordinary high-throughput
  exterior rays.
