# Solar Relativity Phase 5 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an installable, chart-aware, unpolarized invariant
radiative-transfer kernel with analytic Kerr disk/torus matter and bounded
thin-disk surface composition.

**Architecture:** Pure constant-segment transfer and redshift value types form
the L0/L1 base. Trusted evaluation composes untrusted fluid and emission
interfaces without coupling them to geodesic integration. Kerr analytic
matter owns chart-aware kinematics, while thin-disk surface composition remains
separate from volume coefficients so Phase 6 can feed both from accepted rays.

**Tech Stack:** C++17, existing fixed-size relativity math/types,
`Metric`, `KerrChartTransform`, observer/orbit utilities, standard-library
`expm1`, repository Make/CMake discovery, behavior-focused executable tests.

## Global Constraints

- Follow
  `docs/superpowers/specs/2026-07-30-relativity-phase-5-design.md` and
  `SOLAR_RELATIVITY_KERR_完整实施主提示词_v3.md`.
- Preserve signature `(-,+,+,+)`, geometrized spacetime units, coordinate
  order from `Metric::chart()`, and canonical covariant photon momentum.
- Observer-initialized photon momentum is normalized to
  `-p_mu u_obs^mu=1`; require positive caller-supplied `nu_observer`.
- Backward accumulation accepts only `ds=-d lambda>=0`; never reuse the
  future-affine differential sign directly.
- Thin disk is a surface model. Only the analytic torus and explicit volume
  fluids produce volume coefficients.
- Keep metric, Hamiltonian, DOPRI, and event-root layers independent of Phase 5
  types. Do not duplicate the geodesic integrator.
- No camera, renderer, RGB/film grading, CLI transfer command, Planck
  bandpass, polarization, scattering, GRMHD, GPU, WASM, or radiation feedback.
- Constructor/configuration errors throw `std::invalid_argument`; untrusted
  runtime model/sample failures return explicit `TransferError`.
- Keep focused production files near the 200-line review signal where a
  natural ownership boundary exists; do not create one-line wrapper modules.
- Use test-first red/green cycles, exact analytic gates, deliberate sign/chart
  mutations, clean Release, installed-consumer, ASan/UBSan, and Linux/GCC CI.

## File map

| File | Layer | Responsibility |
|---|---|---|
| `include/solar/relativity/radiative_transfer.h` | L0/L1 public | Transfer values, redshift/evaluation results, exact backward advance. |
| `src/relativity/radiative_transfer.cpp` | L1 | Stable constant-segment formal solution and observer intensity conversion. |
| `src/relativity/transfer_evaluation.cpp` | L2 | Validate fluid/frequency and invoke untrusted emission models. |
| `include/solar/relativity/fluid_model.h` | L1 public | Fluid interface, analytic model configs/classes. |
| `include/solar/relativity/emission_model.h` | L1 public | Emission interface plus vacuum/grey/debug models. |
| `src/relativity/emission_model.cpp` | L1 | Invariant coefficient models. |
| `src/relativity/kerr_fluid_kinematics.h` | L1 internal | Common BL/KS point and circular-velocity conversion contract. |
| `src/relativity/kerr_fluid_kinematics.cpp` | L1 internal | Parameter checks, circular normalization, chart transform. |
| `src/relativity/analytic_disk_fluid.cpp` | L1 | Surface-bounded circular disk profile. |
| `src/relativity/analytic_torus_fluid.cpp` | L1 | Compact Gaussian torus profile. |
| `include/solar/relativity/thin_disk.h` | L1/L2 public | Surface emission, crossing value, bounded recorder. |
| `src/relativity/thin_disk_surface.cpp` | L1 | Surface source validation and intensity law. |
| `src/relativity/thin_disk_geometry.{h,cpp}` | L1 internal | BL/KS surface-normal construction and validation. |
| `src/relativity/thin_disk.cpp` | L2 | `g^3/g^4`, opacity composition, and crossing bounds. |
| focused tests below | validation | Independent analytic and chart/failure contracts. |

---

### Task 1: Exact observer-to-past formal solution

**Files:**
- Create: `include/solar/relativity/radiative_transfer.h`
- Create: `src/relativity/radiative_transfer.cpp`
- Create: `tests/relativity/test_radiative_transfer.cpp`

**Interfaces:**
- Consumes: standard floating-point math only.
- Produces: `TransferCoefficients`, `TransferError`,
  `BackwardTransferState`, `TransferAdvanceResult`,
  `advance_backward_transfer`, and `specific_intensity_at_observer`.

- [x] **Step 1: Write the failing analytic test**

Create a local harness and exercise the public API:

```cpp
const BackwardTransferState vacuum{};
const auto emission = advance_backward_transfer(
    vacuum, TransferCoefficients{2.5, 0.0}, 4.0);
check_near("constant emission", emission.state.invariant_intensity,
           10.0, 2.0e-14);

const auto source = advance_backward_transfer(
    vacuum, TransferCoefficients{3.0, 0.4}, 5.0);
check_near("constant source exact",
           source.state.invariant_intensity,
           (3.0 / 0.4) * (1.0 - std::exp(-2.0)), 5.0e-14);
check_near("constant absorption transmission",
           source.state.transmission, std::exp(-2.0), 5.0e-15);
```

Split the same interval into 1, 2, 7, and 100 equal steps; require intensity,
transmission, and optical depth within `5e-14` of the one-step solution.
Include `A*ds=1e-12`, `A*ds=1000`, zero-length, zero coefficients, foreground
attenuation, observer-frequency cubic conversion, and invalid negative/NaN
inputs.

- [x] **Step 2: Run the test and confirm the red state**

Run:

```bash
make -j4 tests/relativity/test_radiative_transfer
```

Expected: compilation fails because `radiative_transfer.h` does not exist.

- [x] **Step 3: Declare exact transfer value types**

Use the design names and keep the state scalar:

```cpp
struct TransferCoefficients {
    double invariant_emissivity = 0.0;
    double invariant_absorption = 0.0;
};

enum class TransferError {
    None,
    NonFiniteInput,
    InvalidObserverFrequency,
    InvalidMetricPoint,
    InvalidFluidSample,
    FourVelocityNotUnitTimelike,
    NonFutureDirectedPhoton,
    EmissionModelFailure,
    InvalidCoefficients,
    InvalidStep,
    NonFiniteResult,
    CrossingLimitReached,
};

struct BackwardTransferState {
    double invariant_intensity = 0.0;
    double transmission = 1.0;
    double optical_depth = 0.0;
};
```

`TransferAdvanceResult::operator bool()` is true only for `None`.

- [x] **Step 4: Implement the stable formal solution**

Validate before calculation. For finite `tau=A*ds` use:

```cpp
double attenuation_integral(double tau, double ds) {
    if (tau == 0.0) {
        return ds;
    }
    if (std::fabs(tau) < 1.0e-6) {
        const double tau2 = tau * tau;
        return ds * (1.0 - 0.5 * tau + tau2 / 6.0);
    }
    return ds * (-std::expm1(-tau) / tau);
}
```

For positive overflowed `tau`, use the exact thick limit `1/A`, set
transmission to zero, and optical depth to positive infinity. Compute:

```text
I_new = I + T * J * attenuation_integral
T_new = T * exp(-tau)
tau_new = tau_old + tau
```

Reject any non-finite intensity result. Accept prior infinite optical depth
only when prior transmission is exactly zero.

- [x] **Step 5: Run focused and regression tests**

Run:

```bash
make -j4 tests/relativity/test_radiative_transfer
./tests/relativity/test_radiative_transfer
make -j4 test
```

Expected: new analytic assertions and all baseline assertions pass.

- [x] **Step 6: Commit**

```bash
git add include/solar/relativity/radiative_transfer.h \
  src/relativity/radiative_transfer.cpp \
  tests/relativity/test_radiative_transfer.cpp
git commit -m "feat(relativity): add invariant transfer kernel"
```

---

### Task 2: Trusted fluid/emission evaluation boundary

**Files:**
- Modify: `include/solar/relativity/radiative_transfer.h`
- Create: `include/solar/relativity/fluid_model.h`
- Create: `include/solar/relativity/emission_model.h`
- Create: `src/relativity/transfer_evaluation.cpp`
- Create: `src/relativity/emission_model.cpp`
- Create: `tests/relativity/test_emission_models.cpp`
- Create: `tests/relativity/test_transfer_failures.cpp`

**Interfaces:**
- Consumes: Task 1 values, `Metric`, spacetime pairing, and
  `observer_measured_frequency`.
- Produces: `FluidSample`, `FluidModel`, `EmissionModel`,
  `VacuumFluid`, `VacuumEmission`, `GreyEmission`, `DebugPaintEmission`,
  `RedshiftSample`, `RedshiftResult`, `TransferEvaluationResult`,
  `evaluate_redshift`, and `evaluate_transfer_coefficients`.

- [x] **Step 1: Write failing frequency and emission tests**

Use Minkowski metric at the origin, photon
`p_mu=(-1,1,0,0)`, and boosted emitter
`u^mu=gamma(1,v,0,0)`:

```text
nu_emit_normalized = gamma (1-v)
g = 1 / [gamma (1-v)]
nu_emit = nu_observer gamma (1-v)
```

For `density=2`, grey `j/rho=3`, `alpha/rho=0.4`, and
`nu_observer=230e9`, require:

```text
J = 6 / nu_emit^2
A = 0.8 nu_emit
```

within `2e-14` relative error. Verify vacuum fluid and vacuum emission both
return exact zero coefficients, while debug emission returns its configured
invariant values.

- [x] **Step 2: Write failing untrusted-boundary tests**

Add local fake models that return or throw:

- valid=false with NaNs in ignored fields: exact vacuum success;
- valid=true with NaN/negative density or temperature;
- spacelike, zero, non-unit, or non-finite velocity;
- thrown fluid/emission exceptions;
- negative/NaN coefficients;
- past-directed photon frequency;
- invalid metric point;
- non-positive/non-finite observer frequency.

Require the exact `TransferError` and a non-empty message. No exception may
escape the evaluation boundary.

- [x] **Step 3: Run both tests and confirm red**

Run:

```bash
make -j4 \
  tests/relativity/test_emission_models \
  tests/relativity/test_transfer_failures
```

Expected: compilation fails because the model headers and evaluation API are
absent.

- [x] **Step 4: Implement the public interfaces**

Declare `FluidSample` and the two virtual interfaces exactly as the design.
Add:

```cpp
struct RedshiftSample {
    double normalized_emitter_frequency;
    double emitter_frequency;
    double redshift_g;
};

struct RedshiftResult {
    TransferError error = TransferError::None;
    RedshiftSample sample;
    std::string message;
    explicit operator bool() const noexcept;
};

struct TransferEvaluationResult {
    TransferError error = TransferError::None;
    TransferCoefficients coefficients;
    RedshiftSample redshift;
    std::string message;
    explicit operator bool() const noexcept;
};

RedshiftResult evaluate_redshift(
    const Metric& metric,
    const PhaseSpaceState& photon,
    const Contravariant4& emitter_four_velocity,
    double observer_frequency);

TransferEvaluationResult evaluate_transfer_coefficients(
    const Metric& metric,
    const PhaseSpaceState& photon,
    double observer_frequency,
    const FluidModel& fluid_model,
    const EmissionModel& emission_model);
```

Forward-declare the model interfaces in `radiative_transfer.h` so L0 values do
not include model implementation headers.

- [x] **Step 5: Implement vacuum, grey, and diagnostic emission**

Constructors reject negative/non-finite configured coefficients.
`GreyEmission::coefficients` computes `nu_emit` from
`-p_mu u^mu * nu_observer` and then the invariant formulas. Vacuum returns
zeros. Debug returns fixed invariant coefficients for a valid sample and is
documented as nonphysical.

- [x] **Step 6: Implement the trusted evaluation sequence**

In `transfer_evaluation.cpp`:

1. validate point, state, and observer frequency;
2. call `fluid_model.sample` inside `try/catch`;
3. return vacuum immediately for `valid=false`;
4. evaluate metric and require velocity norm error `<1e-10`;
5. require positive finite normalized and physical emitter frequency;
6. call emission inside `try/catch`;
7. require finite nonnegative coefficients.

Use `metric_inner_product` and the existing covector/vector pairing. Never
renormalize malformed model velocity.

- [x] **Step 7: Run focused and baseline tests**

Run:

```bash
make -j4 \
  tests/relativity/test_emission_models \
  tests/relativity/test_transfer_failures
./tests/relativity/test_emission_models
./tests/relativity/test_transfer_failures
make -j4 test
```

Expected: all required tests pass.

- [x] **Step 8: Commit**

```bash
git add include/solar/relativity/radiative_transfer.h \
  include/solar/relativity/fluid_model.h \
  include/solar/relativity/emission_model.h \
  src/relativity/transfer_evaluation.cpp \
  src/relativity/emission_model.cpp \
  tests/relativity/test_emission_models.cpp \
  tests/relativity/test_transfer_failures.cpp
git commit -m "feat(relativity): add transfer model boundary"
```

---

### Task 3: Chart-aware Kerr kinematics and circular disk fluid

**Files:**
- Modify: `include/solar/relativity/fluid_model.h`
- Create: `src/relativity/kerr_fluid_kinematics.h`
- Create: `src/relativity/kerr_fluid_kinematics.cpp`
- Create: `src/relativity/analytic_disk_fluid.cpp`
- Create: `tests/relativity/test_fluid_models.cpp`

**Interfaces:**
- Consumes: Task 2 fluid interface, Kerr BL/KS metrics, chart transform,
  `OrbitSense`, and ISCO/orbit formulas.
- Produces: `AnalyticCircularDiskConfig`,
  `AnalyticCircularDiskFluid`, plus internal
  `detail::KerrFluidPoint` and `detail::evaluate_kerr_circular_fluid_point`.

- [x] **Step 1: Write failing disk-profile tests**

Configure `M=2`, `chi=0.5`, prograde, `r_out=30M`,
`density_scale=4`, `temperature_scale=10`, `density_power=1.5`.
Require:

- default `r_in` equals `kerr_isco_radius`;
- samples at `r_in`, `2*r_in`, and `r_out` follow the literal design profile;
- inside/above radial bounds and outside surface tolerance return vacuum;
- temperature is zero at the zero-torque inner boundary;
- invalid mass/spin/radii/scales/power/tolerance throw;
- both spin signs and orbit senses produce finite future unit-timelike
  velocities.

- [x] **Step 2: Add failing BL/KS common-event tests**

Transform at least six equatorial safe-overlap points spanning
`chi={-0.7,0,0.6}` and both senses. Sample the same event in BL and KS.
Require density and temperature relative error `<1e-12`, norm error `<1e-10`,
and:

```cpp
const Mat4 J =
    transform.boyer_lindquist_to_kerr_schild_jacobian(bl_point);
const Contravariant4 expected_ks_u = multiply_vector(J, bl_sample.four_velocity);
```

component maximum error `<1e-10`. Mismatched model/metric Kerr parameters and
unsupported metrics throw `std::domain_error` directly and become
`InvalidFluidSample` through `evaluate_transfer_coefficients`.

- [x] **Step 3: Run the fluid test and confirm red**

Run:

```bash
make -j4 tests/relativity/test_fluid_models
```

Expected: compilation fails because the disk config/class are absent.

- [x] **Step 4: Implement shared Kerr fluid kinematics**

Declare the public disk configuration exactly:

```cpp
struct AnalyticCircularDiskConfig {
    double mass_M;
    double spin_chi;
    OrbitSense sense = OrbitSense::Prograde;
    std::optional<double> inner_radius_M;
    double outer_radius_M;
    double density_scale;
    double temperature_scale;
    double density_power = 0.0;
    double surface_height_tolerance = 1.0e-8;
};

class AnalyticCircularDiskFluid final : public FluidModel {
public:
    explicit AnalyticCircularDiskFluid(
        AnalyticCircularDiskConfig config);
    FluidSample sample(
        const Metric& metric,
        const Contravariant4& x) const override;
    double inner_radius() const noexcept;
    double outer_radius() const noexcept;
};
```

`detail::KerrFluidPoint` stores BL position, caller-chart position,
caller-chart four-velocity, radius, theta, and equatorial height.

For BL, use the point directly. For KS, use `KerrChartTransform` to recover BL.
Evaluate:

```text
Omega = rotation_sign /
        [M ((r/M)^(3/2) + rotation_sign chi)]
u^t = 1/sqrt[-(g_tt+2 Omega g_tphi+Omega^2 g_phiphi)]
u^phi = Omega u^t
```

Then transform `u_KS=J u_BL`. Reject non-finite normalization, non-timelike
worldlines and unsafe overlap with an invalid internal result. Throw
`std::domain_error` for unsupported metric classes or model/metric parameter
mismatch; do not spatially copy velocity.

- [x] **Step 5: Implement the disk profile**

Use `std::optional<double> inner_radius_M` in the config. Resolve an absent
value to ISCO at construction. Sampling checks:

```text
r_in <= r <= r_out
abs(cos(theta)) <= surface_height_tolerance
```

then evaluates the exact profile from the design. Return `valid=true` even at
the zero-temperature inner edge because it is a legitimate disk boundary.

- [x] **Step 6: Run focused and baseline tests**

Run:

```bash
make -j4 tests/relativity/test_fluid_models
./tests/relativity/test_fluid_models
make -j4 test
```

Expected: all assertions pass.

- [x] **Step 7: Commit**

```bash
git add include/solar/relativity/fluid_model.h \
  src/relativity/kerr_fluid_kinematics.h \
  src/relativity/kerr_fluid_kinematics.cpp \
  src/relativity/analytic_disk_fluid.cpp \
  tests/relativity/test_fluid_models.cpp
git commit -m "feat(relativity): add analytic Kerr disk fluid"
```

---

### Task 4: Compact analytic optically thin torus

**Files:**
- Modify: `include/solar/relativity/fluid_model.h`
- Create: `src/relativity/analytic_torus_fluid.cpp`
- Modify: `tests/relativity/test_fluid_models.cpp`

**Interfaces:**
- Consumes: Task 3 internal Kerr kinematics.
- Produces: `AnalyticOpticallyThinTorusConfig` and
  `AnalyticOpticallyThinTorus`.

- [x] **Step 1: Add failing torus tests**

Configure:

```text
M=1, chi=0.6, prograde
r_center=8, radial_width=2, angular_width=0.2
density_scale=5, temperature_scale=7
temperature_power=0.4, density_cutoff_fraction=1e-4
```

Require exact center density/temperature, north/south reflection, the literal
Gaussian value one radial width away, compact vacuum below cutoff, finite
future unit-timelike velocity, BL/KS common-event agreement, negative spin,
and explicit unsupported/mismatched metric errors. Reject every non-finite or
out-of-domain config field. Model/metric errors must throw directly and become
`InvalidFluidSample` at the trusted evaluation boundary.

- [x] **Step 2: Run and confirm the new red state**

Run:

```bash
make -j4 tests/relativity/test_fluid_models
```

Expected: compilation fails because the torus config/class are absent.

- [x] **Step 3: Implement the exact compact Gaussian**

Declare:

```cpp
struct AnalyticOpticallyThinTorusConfig {
    double mass_M;
    double spin_chi;
    OrbitSense sense = OrbitSense::Prograde;
    double center_radius_M;
    double radial_width_M;
    double angular_width;
    double density_scale;
    double temperature_scale;
    double temperature_power = 0.0;
    double density_cutoff_fraction = 1.0e-4;
};

class AnalyticOpticallyThinTorus final : public FluidModel {
public:
    explicit AnalyticOpticallyThinTorus(
        AnalyticOpticallyThinTorusConfig config);
    FluidSample sample(
        const Metric& metric,
        const Contravariant4& x) const override;
};
```

Use:

```cpp
const double radial_offset =
    (radius - center_radius_M) / radial_width_M;
const double angular_offset =
    std::cos(theta) / angular_width;
const double shape = std::exp(
    -0.5 * (radial_offset * radial_offset +
            angular_offset * angular_offset));
```

Return vacuum when `shape<density_cutoff_fraction`; otherwise use the design
density/temperature powers and the shared caller-chart circular velocity.

- [x] **Step 4: Run focused and baseline tests**

Run:

```bash
make -j4 tests/relativity/test_fluid_models
./tests/relativity/test_fluid_models
make -j4 test
```

Expected: all assertions pass.

- [x] **Step 5: Commit**

```bash
git add include/solar/relativity/fluid_model.h \
  src/relativity/analytic_torus_fluid.cpp \
  tests/relativity/test_fluid_models.cpp
git commit -m "feat(relativity): add analytic Kerr torus"
```

---

### Task 5: Thin-disk surface redshift and bounded composition

**Files:**
- Create: `include/solar/relativity/thin_disk.h`
- Create: `src/relativity/thin_disk.cpp`
- Create: `src/relativity/thin_disk_surface.cpp`
- Create: `src/relativity/thin_disk_geometry.h`
- Create: `src/relativity/thin_disk_geometry.cpp`
- Create: `tests/relativity/test_thin_disk.cpp`

**Interfaces:**
- Consumes: Tasks 1–3 transfer, disk fluid, redshift, metric, and chart
  transform contracts.
- Produces: `DiskOpacityMode`, `ThinDiskSurfaceEmission`,
  `ThinDiskCrossing`, `ThinDiskObservedState`, `ThinDiskRecorderConfig`,
  `ThinDiskRecordResult`, and `ThinDiskCrossingRecorder`.

- [x] **Step 1: Write failing one-crossing/redshift tests**

Use an exact disk surface state and controlled photon/emitter pairing with
known `g=0.5`. For temperature `8`, emitted specific intensity `6`, emitted
bolometric intensity `10`, require:

```text
T_obs = 4
I_specific_obs = 0.125 * 6
I_bolometric_obs = 0.0625 * 10
```

Check affine, position, radius, velocity, emitter frequency, `g`, unit
spacelike normal, `u dot n=0`, front/back, and image order zero.

- [x] **Step 2: Write failing multiple/opacity tests**

For two sheet sources `S1,S2`, depths `tau1,tau2` in observer-to-past order,
require:

```text
I = S1(1-exp(-tau1))
  + exp(-tau1) S2(1-exp(-tau2))
T = exp(-(tau1+tau2))
```

Require opaque mode closes after the first valid crossing with full source
intensity and zero transmission. Record exactly eight semi-transparent
crossings, reject the ninth as `CrossingLimitReached`, and preserve the first
eight. Invalid/out-of-radius disk samples do not consume image order.

- [x] **Step 3: Add failing BL/KS surface agreement**

At common disk events require:

- radius and redshift agreement `<1e-10`;
- transformed normal component error `<1e-10`;
- each normal norm equals `+1` within `1e-10`;
- each normal is orthogonal to emitter velocity within `1e-10`;
- identical front/back classification.

- [x] **Step 4: Run and confirm red**

Run:

```bash
make -j4 tests/relativity/test_thin_disk
```

Expected: compilation fails because `thin_disk.h` is absent.

- [x] **Step 5: Implement surface emission and normal**

Use the exact recorder surface:

```cpp
struct ThinDiskObservedState {
    double specific_intensity = 0.0;
    double bolometric_intensity = 0.0;
    double transmission = 1.0;
};

struct ThinDiskRecorderConfig {
    DiskOpacityMode opacity_mode = DiskOpacityMode::Opaque;
    std::size_t max_crossings = 8;
};

class ThinDiskCrossingRecorder {
public:
    ThinDiskCrossingRecorder(
        ThinDiskRecorderConfig config,
        AnalyticCircularDiskFluid disk,
        ThinDiskSurfaceEmission emission);
    ThinDiskRecordResult record(
        const Metric& metric,
        const PhaseSpaceState& photon,
        double observer_frequency);
    const std::vector<ThinDiskCrossing>& crossings() const noexcept;
    const ThinDiskObservedState& observed() const noexcept;
    bool closed() const noexcept;
};
```

`ThinDiskSurfaceEmission` validates nonnegative finite intensity scales and
surface depth. Its emitted values are:

```text
specific = specific_scale * temperature
bolometric = bolometric_scale * temperature^4
```

Build the north normal in BL with negative theta direction, normalize it
against the BL metric, and transform it contravariantly for KS. Refuse any
non-unit/non-orthogonal result instead of repairing it.

- [x] **Step 6: Implement the bounded recorder**

Reserve exactly `max_crossings` in the constructor. Do not grow beyond it.
For each valid disk sample:

1. compute redshift through `evaluate_redshift`;
2. evaluate surface source;
3. apply `g`, `g^3`, and `g^4`;
4. classify face with `p_mu n^mu`;
5. append crossing with current vector size as image order;
6. update both linear observed intensities and their shared transmission.

Opaque mode applies the full observed source and closes. Semi-transparent mode
uses `-expm1(-surface_tau)`. Surface depth positive infinity is allowed only
for opaque saturation.

- [x] **Step 7: Run focused and baseline tests**

Run:

```bash
make -j4 tests/relativity/test_thin_disk
./tests/relativity/test_thin_disk
make -j4 test
```

Expected: all assertions pass.

- [x] **Step 8: Commit**

```bash
git add include/solar/relativity/thin_disk.h \
  src/relativity/thin_disk.cpp \
  tests/relativity/test_thin_disk.cpp
git commit -m "feat(relativity): add thin-disk surface transfer"
```

---

### Task 6: Acceptance cross-checks and mutation sensitivity

**Files:**
- Modify: `tests/relativity/test_radiative_transfer.cpp`
- Modify: `tests/relativity/test_fluid_models.cpp`
- Modify: `tests/relativity/test_thin_disk.cpp`
- Modify if a real defect is found: focused Phase 5 implementation above.

**Interfaces:**
- Consumes: complete Phase 5 public API.
- Produces: measured acceptance output and tests sensitive to the two highest
  risk mutations.

- [ ] **Step 1: Add deterministic measured summaries**

Print only final maximum values:

```text
max_constant_solution_error
max_subdivision_error
max_redshift_error
max_bl_ks_velocity_error
max_velocity_norm_error
max_surface_normal_error
max_surface_composition_error
```

Use `std::setprecision(17)` so the validation report records real values.

- [ ] **Step 2: Run the five focused executables**

Run:

```bash
./tests/relativity/test_radiative_transfer
./tests/relativity/test_emission_models
./tests/relativity/test_transfer_failures
./tests/relativity/test_fluid_models
./tests/relativity/test_thin_disk
```

Expected: every assertion passes and each maximum is below its design gate.

- [ ] **Step 3: Perform the backward-sign mutation**

Temporarily change `exp(-tau)`/`expm1(-tau)` to the wrong positive exponent.
Rebuild and run `test_radiative_transfer` and `test_thin_disk`; both must fail.
Restore with `apply_patch`, rebuild, and require green.

- [ ] **Step 4: Perform the chart-velocity mutation**

Temporarily replace `u_KS=J u_BL` with a component copy. Rebuild and run
`test_fluid_models` and `test_thin_disk`; both must fail. Restore with
`apply_patch`, rebuild, and require green.

- [ ] **Step 5: Self-audit architecture and failure meaning**

Verify with focused searches and diff review:

- generic metric/geodesic/numerics modules do not import transfer headers;
- disk never returns volume opacity;
- torus is explicitly kinematic and compact;
- all model callbacks are caught at the trusted boundary;
- invalid samples never become vacuum; only `valid=false` material support is
  vacuum, while unsupported metric/chart combinations fail explicitly;
- every chart vector uses the full Jacobian;
- no sRGB/RGB, renderer, CLI, GPU, or hidden unit conversion entered Phase 5;
- files over the review signal still own one cohesive responsibility;
- no unrelated edit, stale status item, dead code, or fallback remains.

- [ ] **Step 6: Commit measured test improvements if changed**

```bash
git add tests/relativity/test_radiative_transfer.cpp \
  tests/relativity/test_fluid_models.cpp \
  tests/relativity/test_thin_disk.cpp
git commit -m "test(relativity): crosscheck Phase 5 transfer"
```

Do not create an empty commit and never commit either mutation.

---

### Task 7: Installed public consumer

**Files:**
- Modify: `tests/external_consumer/probe.cpp`

**Interfaces:**
- Consumes: installed Task 1–5 public headers and symbols only.
- Produces: clean external construction and execution evidence.

- [ ] **Step 1: Extend the installed consumer**

Include:

```cpp
#include "solar/relativity/emission_model.h"
#include "solar/relativity/fluid_model.h"
#include "solar/relativity/radiative_transfer.h"
#include "solar/relativity/thin_disk.h"
```

Construct a Minkowski constant transfer with a known analytic result and a
Kerr disk/torus sample. Record JSON fields:

```text
transfer_intensity
transfer_transmission
disk_temperature
torus_density
surface_specific_intensity
surface_crossings
```

Return nonzero unless all are finite, the analytic transfer error is
`<5e-14`, and one surface crossing is recorded.

- [ ] **Step 2: Run the installed consumer**

Run:

```bash
make test-external-consumer
```

Expected: a clean temporary install configures, builds, links, runs, prints the
new fields, and exits zero.

- [ ] **Step 3: Commit**

```bash
git add tests/external_consumer/probe.cpp
git commit -m "test(relativity): consume Phase 5 transfer API"
```

---

### Task 8: Validation, sanitizer, documentation, and release candidate

**Files:**
- Create: `docs/validation/relativity_10_radiative_transfer.md`
- Modify: `RELATIVITY_STATUS.md`
- Modify only if evidence requires it: focused Phase 5 implementation/tests.

**Interfaces:**
- Consumes: all Phase 5 measured gates.
- Produces: reproducible evidence and a clean release candidate.

- [ ] **Step 1: Run clean Release verification**

Run:

```bash
make clean
make -j4 test
make test-external-consumer
git diff --check
```

Record exact relativity, legacy, total, and Phase 5 focused assertion counts.
Record the optional DE440 skip separately. Any required failure keeps Phase 5
incomplete.

- [ ] **Step 2: Run combined ASan/UBSan**

After `make clean`, build the five focused executables once with:

```bash
make CXXFLAGS='-std=c++17 -O1 -g -Wall -Wextra -Iinclude \
  -fsanitize=address,undefined -fno-omit-frame-pointer' \
  tests/relativity/test_radiative_transfer \
  tests/relativity/test_emission_models \
  tests/relativity/test_transfer_failures \
  tests/relativity/test_fluid_models \
  tests/relativity/test_thin_disk
```

Run each with:

```bash
ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 \
UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
  ./tests/relativity/test_radiative_transfer
```

Repeat for the other four and require no diagnostic.

- [ ] **Step 3: Write measured validation evidence**

`docs/validation/relativity_10_radiative_transfer.md` must include:

- claim and precise non-claims;
- invariant equations and backward sign;
- constant-emission/absorption/source analytic values;
- measured maxima and thresholds;
- grey/redshift formula and boosted fixture;
- BL/KS disk/torus and surface-normal gates;
- one/multiple/opaque/semi-transparent results;
- exact commands, compiler, build type, assertion totals, skips;
- installed consumer JSON;
- both mutation failures;
- limitations and fastest falsification commands;
- RAPTOR, ipole, and EHT verification references.

Update every `RELATIVITY_STATUS.md` field. Before Linux CI use
`LOCAL_PASSED_AWAITING_LINUX_CI`; move to `PASSED` only after the PR head CI
succeeds.

- [ ] **Step 4: Run the final release-candidate set**

Restore normal Release objects, then:

```bash
make clean
make -j4 test
make test-external-consumer
git diff --check
git status --short
```

Expected: only intentional Phase 5 files are changed and every required gate
passes.

- [ ] **Step 5: Commit**

```bash
git add docs/validation/relativity_10_radiative_transfer.md \
  RELATIVITY_STATUS.md \
  docs/superpowers/plans/2026-07-30-relativity-phase-5.md
git commit -m "docs(relativity): validate Phase 5 transfer"
```

- [ ] **Step 6: Publish, validate Linux, and merge**

Push `codex/relativity-phase-5`, open one public PR against `main`, and run the
single intentional Linux/GCC release-candidate CI. If CI exposes a real
platform defect, add a focused regression and fix the cause. After the PR head
is green, record the Linux evidence in the status/report, revalidate the final
head, merge, and require the automatic `main` CI to pass.
