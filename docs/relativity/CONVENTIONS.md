# Solar Relativity Conventions

This document is normative for `solar::relativity`. Code and later validation
reports must use these conventions unless a phase-gated API revision changes
both together.

## Model boundary

The relativity kernel models test photons and test particles on a fixed,
analytic background spacetime. The particle and radiation do not alter the
metric. Phase 0A provides only mathematical and unit foundations; it does not
yet implement a metric, geodesic, horizon, observer, radiative transfer, or
black-hole image.

It is not a solver for the Einstein equations, dynamical spacetime, black-hole
mergers, self-force, radiation reaction, GRMHD, or quantum gravity.

## Indices, signature, and coordinates

- Spacetime dimension: 4.
- Metric signature: `(-,+,+,+)`.
- Greek coordinate indices `mu,nu,alpha,beta` range over `{0,1,2,3}`.
- Local orthonormal-frame indices `(a),(b)` range over `{0,1,2,3}`.
- Repeated upper/lower indices imply Einstein summation in equations.
- Component `x^0` is always first.

Chart component order:

- Minkowski Cartesian: `(t,x,y,z)`.
- Boyer–Lindquist: `(t,r,theta,phi)`.
- Cartesian Kerr–Schild: `(t,x,y,z)`.

The Kerr spin axis is world `+z`. A frontend that displays `+y` as vertical
must convert axes only in its display adapter; physics data remains `+z`
aligned.

## Units and black-hole parameters

Integration uses geometrized units:

```text
G = c = 1
```

For SI mass `M_kg`:

```text
M_length = G M_kg / c^2  [metres]
M_time   = G M_kg / c^3  [seconds]
```

`GeometricUnits::from_solar_masses` uses
`1 solar mass = 1.98847e30 kg`. At one solar mass this implementation gives:

```text
M_length = 1476.6696910334392 m
M_time   = 4.925639893961039e-6 s
```

Solar's stored constants are in km-based units. `units.cpp` explicitly converts
`G` from km³/(kg·s²) to m³/(kg·s²) and `C_LIGHT` from km/s to m/s before using
the equations above.

Later Kerr parameters use:

- `M`: the geometrized mass length;
- `a = J/(M c)`: spin length;
- `chi = a/M`: dimensionless spin;
- routine first-version input requires `abs(chi) < 1`.

Phase 0A does not yet construct or validate Kerr parameters.

## Mathematical types

`Vec<N>` is an unlabelled fixed-size numeric container used inside explicit
contracts. `Vec3` and `Vec4` in `solar::relativity` are not implicitly
convertible to legacy `solar::Vec3`, whose components carry Newtonian km or
km/s semantics.

`Contravariant4` and `Covariant4` are distinct public wrappers. Code must not
raise or lower an index by copying components; a metric must perform that
operation in a later phase.

There is no default Euclidean `dot(Vec4,Vec4)`. Minkowski contraction is
available only through the deliberately named
`minkowski_dot_minus_plus_plus_plus`. General contraction must receive an
explicit metric.

`Mat4` is row-major. `inverse(Mat4)` uses partial-pivot Gauss–Jordan
elimination and throws `std::domain_error` for non-finite, singular, or
double-precision ill-conditioned input. It does not regularize a matrix or
silently replace a physical result.

## Phase-space state

The v3 public contract is:

```cpp
struct PhaseSpaceState {
    double affine;
    Contravariant4 x;
    Covariant4 p;
};
```

`x` lists coordinate components in the active chart order. Although a
coordinate tuple is mathematically a point rather than a vector, this wrapper
choice is retained to match the approved v3 API. Any later separation into a
coordinate-point type is an explicit migration.

`p` is covariant canonical momentum. It must not be replaced with a
three-velocity or a normalized spatial direction.

Geodesic kinds and Hamiltonian targets:

```text
Null:             H = 0
TimelikeUnitMass: H = -1/2
```

For a timelike unit-mass geodesic, the affine parameter may be proper time. A
null affine parameter has no unique absolute scale; observer initialization
will later set `nu_obs = 1`.

`GeodesicSample::proper_time` defaults to quiet NaN. Null paths must retain NaN
instead of reporting coordinate time as proper time.

## Dual4 automatic differentiation

`Dual4` stores one scalar value and derivatives with respect to the four
ordered coordinates. `Dual4::variable(value, mu)` seeds derivative `mu` with
one. Arithmetic and `sqrt`, `sin`, `cos`, `log`, and `atan2` apply forward
first-order differentiation.

Domain behavior:

- variable indices outside `[0,3]` are rejected;
- division by a zero or non-finite value is rejected;
- a quotient that produces non-finite value/derivative is rejected;
- square root rejects negative/non-finite values;
- `sqrt(0)` is valid for a constant zero but rejects a nonzero seeded
  derivative because the derivative is singular;
- logarithm requires a finite positive value;
- `atan2(0,0)` and non-finite radius are rejected.

No function clamps its input to make a result finite.

## Failure and finite-value policy

Construction or conversion input outside an API's mathematical domain throws
`std::invalid_argument` or `std::domain_error` as documented. Ordinary future
geodesic outcomes such as capture, escape, chart exit, or step-budget
exhaustion will use explicit result enums instead of exceptions.

`all_finite` is diagnostic, not permission to hide bad values. A caller must
surface non-finite state and stop the affected numerical path.
