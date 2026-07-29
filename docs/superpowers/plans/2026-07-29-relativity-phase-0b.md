# Solar Relativity Phase 0B Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add independently validated Minkowski, Schwarzschild, and Kerr
Boyer–Lindquist metrics plus a real metric-inspection CLI.

**Architecture:** Metric implementations are L1 building blocks over Phase 0A
math/types. The CLI is an L3 adapter in a separate translation unit. Analytic
metric formulas, inverse formulas, and derivatives are cross-checked through
independent identities, Dual4 expressions, and five-point finite differences.

**Tech Stack:** C++17, GNU Make, standalone executable tests, existing
`solar::relativity` L0 types.

## Global Constraints

- Signature is `(-,+,+,+)` and coordinate zero is time.
- CPU C++17 `double` is authoritative.
- Geometrized units use `G=c=1`.
- Kerr uses `a=M*chi` with finite `M>0` and `abs(chi)<1`.
- Ordinary-point inverse error must be below `5e-13`.
- BL derivative comparison normalized error must be below `1e-8`.
- BL horizon and polar singularities are rejected, never clamped.
- No geodesic, observer, Kerr–Schild, renderer, adapter, or UI code is added.

---

### Task 1: Metric interface and Minkowski implementation

**Files:**
- Create: `include/solar/relativity/metric.h`
- Create: `include/solar/relativity/minkowski_metric.h`
- Create: `src/relativity/minkowski_metric.cpp`
- Create: `tests/relativity/test_metrics.cpp`

**Interfaces:**
- Consumes: `Contravariant4`, `Mat4`.
- Produces: `Chart`, `chart_name`, abstract `Metric`, and `MinkowskiMetric`.

- [ ] **Step 1: Write the missing-interface test**

Require exact diagonal and zero derivatives:

```cpp
MinkowskiMetric metric;
const Contravariant4 x{Vec4{{2.0, -3.0, 4.0, 5.0}}};
check_near("g_tt", metric.covariant(x)[0][0], -1.0, 0.0);
check_near("g_xx", metric.contravariant(x)[1][1], 1.0, 0.0);
check("all derivatives zero",
      max_abs(metric.contravariant_derivatives(x)) == 0.0);
```

Also require a non-finite coordinate to make `valid_point` false and matrix
evaluation throw `std::domain_error`.

- [ ] **Step 2: Verify RED**

Run:

```bash
make tests/relativity/test_metrics
```

Expected: compilation fails because the metric headers do not exist.

- [ ] **Step 3: Implement the exact v3 interface and Minkowski**

Use:

```cpp
class Metric {
public:
    virtual ~Metric() = default;
    virtual Chart chart() const noexcept = 0;
    virtual std::string name() const = 0;
    virtual Mat4 covariant(const Contravariant4& x) const = 0;
    virtual Mat4 contravariant(const Contravariant4& x) const = 0;
    virtual std::array<Mat4, 4>
        contravariant_derivatives(const Contravariant4& x) const = 0;
    virtual bool valid_point(const Contravariant4& x) const noexcept = 0;
};
```

Return literal diagonal matrices; do not call matrix inversion.

- [ ] **Step 4: Verify GREEN and commit**

```bash
make tests/relativity/test_metrics
./tests/relativity/test_metrics
git add include/solar/relativity/metric.h \
  include/solar/relativity/minkowski_metric.h \
  src/relativity/minkowski_metric.cpp tests/relativity/test_metrics.cpp
git commit -m "feat: add metric interface and Minkowski spacetime"
```

### Task 2: Schwarzschild Boyer–Lindquist metric

**Files:**
- Create: `include/solar/relativity/schwarzschild_metric.h`
- Create: `src/relativity/schwarzschild_metric.cpp`
- Modify: `tests/relativity/test_metrics.cpp`

**Interfaces:**
- Consumes: `Metric`, matrix multiplication, geometrized `M`.
- Produces: `SchwarzschildBoyerLindquistMetric`, analytic covariant and inverse
  derivatives, horizon radius, and BL validity checks.

- [ ] **Step 1: Add literal Schwarzschild tests**

At `M=1,r=10,theta=pi/2`, require:

```cpp
g_tt = -0.8
g_rr = 1.25
g_thetatheta = 100
g_phiphi = 100
g_inverse_tt = -1.25
g_inverse_rr = 0.8
g_inverse_thetatheta = 0.01
g_inverse_phiphi = 0.01
```

Require inverse identity error below `5e-13`; reject `M<=0`, `r<=2M+margin`,
the polar axis, theta outside `(0,pi)`, and non-finite input.

- [ ] **Step 2: Verify RED**

```bash
make tests/relativity/test_metrics
```

Expected: compilation fails because `schwarzschild_metric.h` is absent.

- [ ] **Step 3: Implement explicit metric and analytic derivative**

Use `f=1-2M/r`, explicit covariant/inverse diagonal matrices, analytic
`partial_r` and `partial_theta` covariant entries, then:

```cpp
partial_inverse[k] =
    -multiply(multiply(g_inverse, partial_covariant[k]), g_inverse);
```

- [ ] **Step 4: Verify GREEN and commit**

```bash
make tests/relativity/test_metrics
./tests/relativity/test_metrics
git add include/solar/relativity/schwarzschild_metric.h \
  src/relativity/schwarzschild_metric.cpp tests/relativity/test_metrics.cpp
git commit -m "feat: add Schwarzschild Boyer-Lindquist metric"
```

### Task 3: Kerr Boyer–Lindquist metric and surfaces

**Files:**
- Create: `include/solar/relativity/kerr_bl_metric.h`
- Create: `src/relativity/kerr_bl_metric.cpp`
- Create: `tests/relativity/test_kerr_bl.cpp`

**Interfaces:**
- Consumes: `Metric`, explicit matrix operations, `M`, `spin_chi`.
- Produces: Kerr covariant/inverse matrices, analytic derivatives, horizon and
  stationary-limit surfaces, and checked BL domain.

- [ ] **Step 1: Write literal Kerr and limit tests**

For `M=1,chi=0.9,r=10,theta=pi/2`, compare every nonzero matrix component to
precomputed literals independent of production code. Require:

```cpp
max_abs(g_covariant * g_contravariant - identity) < 5e-13
outer_horizon_radius() == 1.4358898943540672 within 1e-15
inner_horizon_radius() == 0.5641101056459328 within 1e-15
outer_stationary_limit_radius(pi/2) == 2 within 1e-15
```

At several exterior points compare `Kerr(M,0)` component by component with
`Schwarzschild(M)`.

- [ ] **Step 2: Verify RED**

```bash
make tests/relativity/test_kerr_bl
```

Expected: compilation fails because `kerr_bl_metric.h` is absent.

- [ ] **Step 3: Implement explicit Kerr matrices**

Implement exact v3 `Sigma`, `Delta`, `A`, covariant entries, inverse entries,
and analytic covariant derivatives. Obtain inverse derivatives only from the
matrix identity.

- [ ] **Step 4: Implement checked surfaces and domain**

Reject non-finite/non-positive mass, `abs(chi)>=1`, non-finite coordinates,
axis points, `Sigma<=0`, radii at/below outer horizon plus margin, and Delta
below its double-precision scale floor.

- [ ] **Step 5: Verify GREEN and commit**

```bash
make tests/relativity/test_kerr_bl
./tests/relativity/test_kerr_bl
git add include/solar/relativity/kerr_bl_metric.h \
  src/relativity/kerr_bl_metric.cpp tests/relativity/test_kerr_bl.cpp
git commit -m "feat: add Kerr Boyer-Lindquist metric"
```

### Task 4: Independent derivative validation

**Files:**
- Create: `tests/relativity/test_metric_derivatives.cpp`
- Modify when a failing test proves necessary:
  `src/relativity/schwarzschild_metric.cpp`
- Modify when a failing test proves necessary:
  `src/relativity/kerr_bl_metric.cpp`

**Interfaces:**
- Consumes: analytic metric derivatives, Dual4, metric evaluation.
- Produces: independent derivative evidence at multiple exterior points.

- [ ] **Step 1: Write independent Dual4 covariant expressions**

In the test only, evaluate the v3 Kerr covariant equations with coordinate
`r` and `theta` seeded as Dual4 variables. Compare every analytic covariant
partial with the resulting derivative using normalized error:

```cpp
error = abs(actual - expected) /
        max({1.0, abs(actual), abs(expected)})
```

Require maximum below `1e-12` at ordinary points.

- [ ] **Step 2: Write five-point inverse finite differences**

For coordinate `k` use:

```text
[-g(x+2h)+8g(x+h)-8g(x-h)+g(x-2h)] / (12h)
```

Compare against `contravariant_derivatives` below `1e-8`, including exact zero
stationary/axisymmetric derivatives.

- [ ] **Step 3: Verify RED against deliberate mutation**

Temporarily negate one copied expected analytic derivative in the test and
confirm the executable exits nonzero; restore it before implementation fixes.

- [ ] **Step 4: Run and correct production only for demonstrated failures**

```bash
make tests/relativity/test_metric_derivatives
./tests/relativity/test_metric_derivatives
```

Expected final result: Dual4 and five-point comparisons pass at all sampled
ordinary exterior points.

- [ ] **Step 5: Commit validation**

```bash
git add tests/relativity/test_metric_derivatives.cpp \
  src/relativity/schwarzschild_metric.cpp src/relativity/kerr_bl_metric.cpp
git commit -m "test: cross-check Boyer-Lindquist metric derivatives"
```

### Task 5: Real metric CLI

**Files:**
- Create: `cli/relativity_metric.h`
- Create: `cli/relativity_metric.cpp`
- Modify: `cli/main.cpp`
- Modify: `Makefile`
- Create: `tests/relativity/test_metric_cli.cpp`

**Interfaces:**
- Consumes: `Metric` implementations and external argv strings.
- Produces: `dispatch_relativity`, human/JSON metric inspection, truthful exit
  codes, recursively compiled CLI translation units.

- [ ] **Step 1: Write real-process CLI tests**

Run the actual binary from the test:

```bash
./solar relativity metric --metric kerr-bl --M 1 --spin 0.9 \
  --x 0,10,1.5707963267948966,0 --json
```

Require exit 0 and output fields `"metric":"kerr-bl"` and
`"inverse_error"`. Require nonzero exits for missing `--x`, `--spin 1`,
unknown metric, malformed coordinate list, and unknown option.

- [ ] **Step 2: Verify RED**

```bash
make tests/relativity/test_metric_cli
./tests/relativity/test_metric_cli
```

Expected: test fails because the command is unknown.

- [ ] **Step 3: Make CLI translation units recursive**

Replace the single `CLI_SRC` link with sorted recursive `CLI_SRCS`,
`CLI_OBJS`, dependency files, and a `cli/%.o` rule. Make `test` depend on
`solar` so real-process tests have the executable.

- [ ] **Step 4: Implement strict parsing and output**

Use `std::stod` with full-string consumption and an exact four-value coordinate
parser. Construct the selected metric, verify `valid_point`, compute both
matrices and inverse-identity error, then print JSON or human output. Catching
at top-level `main` preserves nonzero failures.

- [ ] **Step 5: Verify GREEN and commit**

```bash
make clean
make tests/relativity/test_metric_cli
./tests/relativity/test_metric_cli
git add Makefile cli/main.cpp cli/relativity_metric.h \
  cli/relativity_metric.cpp tests/relativity/test_metric_cli.cpp
git commit -m "feat: add relativity metric CLI"
```

### Task 6: Phase 0B documentation and gate

**Files:**
- Create: `docs/validation/relativity_0b_metrics.md`
- Modify: `RELATIVITY_STATUS.md`
- Modify: `docs/superpowers/plans/2026-07-29-relativity-phase-0b.md`

**Interfaces:**
- Consumes: actual test/build/CLI results.
- Produces: reproducible metric evidence and Phase 0B state.

- [ ] **Step 1: Run complete verification**

```bash
make clean
make
make test
./tests/relativity/test_metrics
./tests/relativity/test_kerr_bl
./tests/relativity/test_metric_derivatives
./tests/relativity/test_metric_cli
./solar relativity metric --metric kerr-bl --M 1 --spin 0.9 \
  --x 0,10,1.5707963267948966,0 --json
git diff --check
```

- [ ] **Step 2: Audit the final diff**

Check formulas against v3 line by line, dependency direction, every invalid
domain, absence of clamping, exact zero derivative dimensions, test mutation
coverage, unchanged legacy CLI behavior, and no Phase 1 code.

- [ ] **Step 3: Write validation and status**

Record literal matrix values, maximum inverse error, maximum Dual4 error,
maximum five-point error, commands, platform, skipped external tests, model
boundary, at least three likely bugs, and fastest falsification inputs.

- [ ] **Step 4: Mark and commit the gate**

Set `CURRENT_PHASE: 0B` and `PHASE_STATE: PASSED` only if every command exits
as expected.

```bash
git add docs/validation/relativity_0b_metrics.md RELATIVITY_STATUS.md \
  docs/superpowers/plans/2026-07-29-relativity-phase-0b.md
git commit -m "docs: pass relativity phase 0b gate"
```
