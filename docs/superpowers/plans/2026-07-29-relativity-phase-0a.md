# Solar Relativity Phase −1/0A Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Establish an audited, recursively built, tested C++17 L0 foundation for later fixed-background general relativity in Solar.

**Architecture:** Preserve all legacy Newtonian APIs and add a one-way,
standard-library-only `solar::relativity` L0. Units, numeric containers,
tensor-variance wrappers, and forward automatic differentiation are separate
cohesive modules; tests consume these modules directly.

**Tech Stack:** C++17, GNU Make, Apple Clang/GCC-compatible standard library,
standalone executable tests.

## Global Constraints

- C++17 `double` CPU results are authoritative.
- Geometrized integration units use `G=c=1`.
- Spacetime signature is `(-,+,+,+)` and component zero is time.
- No runtime third-party dependencies and no mandatory CMake migration.
- Existing `libsolar.a`, `solar`, CLI commands, and old tests remain available.
- This plan stops before metrics, geodesics, Kerr, rendering, adapters, and UI.

---

### Task 1: Phase −1 evidence and API audit

**Files:**
- Create: `RELATIVITY_STATUS.md`
- Create: `docs/relativity/SOLAR_API_AUDIT.md`
- Create: `docs/validation/relativity_00_baseline.md`
- Modify: `Makefile`
- Modify: `tests/test_de.cpp`
- Modify: `tests/test_validation.cpp`

**Interfaces:**
- Consumes: existing Solar headers, Make targets, CLI commands, and standalone tests.
- Produces: recursive `SRCS`/`TEST_SRCS`, truthful aggregate test exit status,
  portable baseline tests, and phase evidence.

- [ ] **Step 1: Record the untouched baseline**

Run:

```bash
git status --short
git rev-parse HEAD
make clean
make
make test
```

Expected: build succeeds with existing warnings; the test build fails because
`std::ostringstream` is incomplete in `test_validation.cpp`.

- [ ] **Step 2: Apply the minimum portability correction**

Add:

```cpp
#include <sstream>
```

to `tests/test_validation.cpp`. In `test_de.cpp`, return success with an
explicit `SKIP` only when the default optional fixture is absent; preserve a
nonzero result for an explicitly supplied invalid path.

- [ ] **Step 3: Close recursive build and test failure propagation**

Use:

```make
SRCS := $(shell find src -name '*.cpp' -type f | sort)
TEST_SRCS := $(shell find tests -name 'test_*.cpp' -type f | sort)
```

and accumulate a nonzero child-test status in the `test` recipe before exiting.

- [ ] **Step 4: Verify the repaired old baseline**

Run:

```bash
make clean
make
make test
./tests/test_de /tmp/solar-definitely-missing-de440.asc
./solar bodies
./solar transfer Earth Mars 2026-04-12
./solar lagrange Sun Earth
./solar halo Sun Earth 1 50000
./solar energy 365 3600 --planets-only
./solar lambert Earth Mars 2026-10-26 2027-09-03
```

Expected: old suite passes or explicitly skips default DE440; the explicit
missing DE440 path exits 1; smoke outputs match the recorded validation values.

- [ ] **Step 5: Write the audit, baseline, and initial status**

Document actual interfaces, reusability decisions, observed outputs, warnings,
model boundary, and fastest falsification commands. Set Phase −1 passed only
after every Phase −1 gate has evidence.

- [ ] **Step 6: Commit Phase −1 artifacts**

```bash
git add Makefile tests/test_de.cpp tests/test_validation.cpp \
  RELATIVITY_STATUS.md docs/relativity/SOLAR_API_AUDIT.md \
  docs/validation/relativity_00_baseline.md
git commit -m "docs: record solar relativity phase minus one"
```

### Task 2: Geometric unit conversion

**Files:**
- Create: `include/solar/relativity/units.h`
- Create: `src/relativity/units.cpp`
- Create: `tests/relativity/test_units.cpp`

**Interfaces:**
- Consumes: `solar::constants::G` in km³/(kg·s²) and
  `solar::constants::C_LIGHT` in km/s.
- Produces: `GeometricUnits::from_mass_kg`, `from_solar_masses`,
  `length_si_to_M`, `length_M_to_si`, `time_si_to_M`, `time_M_to_si`,
  and `velocity_si_to_c`.

- [ ] **Step 1: Write the failing unit test**

Exercise one-solar-mass references:

```cpp
const auto units = GeometricUnits::from_solar_masses(1.0);
check_near("solar M length", units.M_length_m, 1476.669691, 1e-6);
check_near("solar M time", units.M_time_s, 4.925639893961039e-6, 1e-18);
check_near("length round trip",
           units.length_M_to_si(units.length_si_to_M(12345.0)),
           12345.0, 1e-10);
```

Also require invalid mass and non-finite conversion input to throw
`std::invalid_argument`.

- [ ] **Step 2: Verify the unit test fails for the missing API**

Run:

```bash
make tests/relativity/test_units
```

Expected: compilation fails because `solar/relativity/units.h` is absent.

- [ ] **Step 3: Implement the minimum checked conversion API**

Compute SI constants from existing kilometre-based Solar constants:

```cpp
constexpr double km3_to_m3 = 1.0e9;
constexpr double km_to_m = 1.0e3;
const double G_si = constants::G * km3_to_m3;
const double c_si = constants::C_LIGHT * km_to_m;
```

Reject non-finite/non-positive mass and non-finite conversion inputs.

- [ ] **Step 4: Verify units**

Run:

```bash
make tests/relativity/test_units
./tests/relativity/test_units
```

Expected: all unit reference, round-trip, and invalid-input checks pass.

- [ ] **Step 5: Commit units**

```bash
git add include/solar/relativity/units.h src/relativity/units.cpp \
  tests/relativity/test_units.cpp
git commit -m "feat: add checked geometric unit conversions"
```

### Task 3: Strong spacetime math and state types

**Files:**
- Create: `include/solar/relativity/math.h`
- Create: `include/solar/relativity/types.h`
- Create: `tests/relativity/test_math.cpp`
- Create: `tests/relativity/test_types.cpp`

**Interfaces:**
- Consumes: standard arrays and numeric functions only.
- Produces: `Vec<N>`, `Vec3`, `Vec4`, `Mat4`, arithmetic, `max_norm`,
  `multiply`, `inverse`, `minkowski_dot_minus_plus_plus_plus`,
  `Contravariant4`, `Covariant4`, `GeodesicKind`, `PhaseSpaceState`, and
  `GeodesicSample`.

- [ ] **Step 1: Write failing vector and tensor-type tests**

Require:

```cpp
const Vec4 a{{1.0, -2.0, 3.0, -4.0}};
const Vec4 b{{0.5, 2.0, -1.0, 8.0}};
check_near("max norm", max_norm(a), 4.0, 0.0);
check_near("Minkowski contraction",
           minkowski_dot_minus_plus_plus_plus(a, b), -39.5, 0.0);
static_assert(!std::is_convertible_v<Contravariant4, Covariant4>);
```

Require identity recovery from a hand-derived nonsingular matrix and
`std::domain_error` for a matrix with duplicate rows.

- [ ] **Step 2: Verify math tests fail for missing APIs**

Run:

```bash
make tests/relativity/test_math tests/relativity/test_types
```

Expected: compilation fails because the new headers are absent.

- [ ] **Step 3: Implement vectors and matrices**

Use partial-pivot Gauss–Jordan inversion with an epsilon-scaled pivot test.
Do not define an unnamed/default four-dimensional Euclidean dot product.
Reject zero/non-finite scalar division and non-finite/singular inversion.

- [ ] **Step 4: Implement strongly wrapped states**

Use the v3 contract:

```cpp
struct Contravariant4 { Vec4 v; };
struct Covariant4 { Vec4 v; };

struct PhaseSpaceState {
    double affine = 0.0;
    Contravariant4 x;
    Covariant4 p;
};
```

Give `GeodesicSample::proper_time` a quiet-NaN default so null paths cannot
accidentally report coordinate time as proper time.

- [ ] **Step 5: Verify math and types**

Run:

```bash
make tests/relativity/test_math tests/relativity/test_types
./tests/relativity/test_math
./tests/relativity/test_types
```

Expected: arithmetic, contraction, inverse, singular-domain, variance, and
proper-time convention checks pass.

- [ ] **Step 6: Commit math and types**

```bash
git add include/solar/relativity/math.h include/solar/relativity/types.h \
  tests/relativity/test_math.cpp tests/relativity/test_types.cpp
git commit -m "feat: add strong spacetime math types"
```

### Task 4: Forward-mode Dual4

**Files:**
- Create: `include/solar/relativity/dual4.h`
- Create: `tests/relativity/test_dual4.cpp`

**Interfaces:**
- Consumes: scalar `double` math.
- Produces: `Dual4::variable`, arithmetic with Dual4/scalars, unary negation,
  `sqrt`, `sin`, `cos`, `log`, `atan2`, and finite checks.

- [ ] **Step 1: Write the failing derivative tests**

For independent variables `x=2`, `y=3`, require the hand-derived result for
`f=x*x*y + sin(y)`:

```cpp
check_near("value", f.value, 12.0 + std::sin(3.0), 1e-14);
check_near("df/dx", f.derivative[0], 12.0, 1e-14);
check_near("df/dy", f.derivative[1], 4.0 + std::cos(3.0), 1e-14);
```

Also test quotient, square-root/log chains, `atan2`, invalid variable index,
zero division, and invalid elementary-function domains.

- [ ] **Step 2: Verify the Dual4 test fails for the missing API**

Run:

```bash
make tests/relativity/test_dual4
```

Expected: compilation fails because `solar/relativity/dual4.h` is absent.

- [ ] **Step 3: Implement the minimum forward-mode algebra**

Represent:

```cpp
struct Dual4 {
    double value = 0.0;
    std::array<double, 4> derivative{};
};
```

Apply product, quotient, and chain rules directly. Reject mathematical domain
violations without clamping.

- [ ] **Step 4: Verify Dual4**

Run:

```bash
make tests/relativity/test_dual4
./tests/relativity/test_dual4
```

Expected: all analytic derivative and failure-domain checks pass.

- [ ] **Step 5: Commit Dual4**

```bash
git add include/solar/relativity/dual4.h tests/relativity/test_dual4.cpp
git commit -m "feat: add four-variable forward autodiff"
```

### Task 5: Conventions, phase gate, and complete verification

**Files:**
- Create: `docs/relativity/CONVENTIONS.md`
- Modify: `RELATIVITY_STATUS.md`
- Modify: `docs/validation/relativity_00_baseline.md`

**Interfaces:**
- Consumes: tested Phase 0A code and actual command results.
- Produces: human-readable contract identical to code and a truthful phase gate.

- [ ] **Step 1: Write conventions matching code**

Document signature, coordinate order, spin axis, geometric units, tensor
variance, state normalization, null proper time, input failures, and the
fixed-background model boundary.

- [ ] **Step 2: Run focused and full verification**

Run:

```bash
make clean
make
make test
./tests/relativity/test_units
./tests/relativity/test_math
./tests/relativity/test_types
./tests/relativity/test_dual4
git diff --check
```

Expected: exit 0 for every command; known legacy compile warnings remain
recorded rather than hidden.

- [ ] **Step 3: Confirm recursive clean behavior**

Run:

```bash
make clean
test ! -e src/relativity/units.o
test ! -e tests/relativity/test_units
make
make test
```

Expected: nested objects and executables are removed, rediscovered, rebuilt,
and executed.

- [ ] **Step 4: Update status from fresh evidence**

Set `CURRENT_PHASE: 0A` and `PHASE_STATE: PASSED` only if every Phase −1 and
0A requirement is present and all required commands pass. Record the verified
commit and platform; list at least three plausible failure modes.

- [ ] **Step 5: Commit phase gate documentation**

```bash
git add docs/relativity/CONVENTIONS.md docs/validation/relativity_00_baseline.md \
  RELATIVITY_STATUS.md
git commit -m "docs: pass relativity phase 0a gate"
```
