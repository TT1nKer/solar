# Solar Relativity Phase 5: Invariant Radiative Transfer Design

## 1. Goal and gate

Phase 5 adds a public, installable, unpolarized covariant radiative-transfer
kernel and controlled analytic matter models for later CPU rendering. It
implements:

- observer-to-past formal-solution accumulation;
- local emission frequency and redshift validation;
- vacuum, grey, and diagnostic volume emission;
- an analytic circular thin-disk surface model;
- an analytic optically thin torus volume model;
- bounded thin-disk crossing records and linear-intensity composition;
- flat-space analytic, redshift, surface, chart, and failure regressions.

Phase 5 passes only when the clean Release suite, focused sanitizers, installed
consumer, analytic transfer gates, BL/KS fluid agreement, and deliberate sign
mutations pass locally and Linux/GCC CI is green.

This phase does not add a camera, pixels, RGB/film grading, a reference image,
CLI transfer command, GPU, WASM, polarized transport, scattering, GRMHD,
returning radiation, disk self-gravity, or radiation feedback. Those remain
separate later milestones.

## 2. Fixed conventions

All existing Solar relativity conventions remain unchanged:

- signature `(-,+,+,+)`;
- geometrized spacetime units `G=c=1`;
- coordinate order reported by `Metric::chart()`;
- canonical state order `(x^mu,p_mu)`;
- photon momentum is the covector `p_mu`;
- observer-initialized photon momentum is normalized so
  `-p_mu u_obs^mu=1`;
- a renderer traces toward the past with decreasing affine parameter and
  supplies the positive segment length `ds=-d lambda`.

For a valid emitter four-velocity:

```text
nu_emit_normalized = -p_mu u_emit^mu
nu_emit = nu_observer * nu_emit_normalized
g = nu_observer / nu_emit
  = 1 / nu_emit_normalized
```

`nu_observer` is a positive caller-selected physical or normalized frequency.
The model does not infer physical density, temperature, or frequency units.
Model documentation must name the units chosen by its caller.

The invariant quantities are:

```text
I = I_nu / nu^3
J = j_nu / nu^2
A = nu alpha_nu
```

Along future affine parameter:

```text
dI/dlambda = J - A I
```

The observer-to-past accumulator does not integrate this equation with a
silently reversed sign. It accumulates foreground-to-background contribution
with:

```text
tau_step = A ds
e = exp(-tau_step)
delta_I = T J ds phi(tau_step)
phi(tau) = (1-exp(-tau))/tau
I <- I + delta_I
T <- T e
```

The continuous limit is `phi(0)=1`. The implementation uses `expm1` and a
small-series branch so optically thin steps do not lose digits. Large optical
depth may set `T` exactly to zero; that is physical saturation, not an error.

## 3. Chosen architecture

The design keeps four independently testable responsibilities:

1. `radiative_transfer` owns invariant scalar algebra, validation results,
   redshift, coefficient evaluation, and exact constant-segment accumulation.
2. `fluid_model` owns material samples and analytic Kerr disk/torus
   kinematics. It never accumulates intensity.
3. `emission_model` converts a valid fluid sample into invariant volume
   coefficients. It never integrates a geodesic or applies display color.
4. `thin_disk` owns surface-only emission, crossing metadata, bounded
   foreground-to-background composition, and face classification.

No Phase 5 type is imported by `Metric`, Hamiltonian, DOPRI, event-root, or
generic geodesic layers. Phase 6 will sample accepted ray steps and pass those
samples into this kernel. Phase 5 therefore does not duplicate the geodesic
integrator or add an unbounded trajectory allocation.

Layering is:

```text
L0:
  radiative transfer value types, enums, pure scalar formulas

L1:
  fluid/emission interfaces, analytic models, Kerr chart kinematics,
  thin-disk surface algebra

L2:
  coefficient evaluation and bounded crossing recorder/compositor

L3:
  installed public headers and external-consumer proof
```

## 4. Public transfer contract

`include/solar/relativity/radiative_transfer.h` declares:

```cpp
struct TransferCoefficients {
    double invariant_emissivity = 0.0;
    double invariant_absorption = 0.0;
};

struct RedshiftSample {
    double normalized_emitter_frequency;
    double emitter_frequency;
    double redshift_g;
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

struct RedshiftResult {
    TransferError error = TransferError::None;
    RedshiftSample sample;
    std::string message;

    explicit operator bool() const noexcept;
};

struct BackwardTransferState {
    double invariant_intensity = 0.0;
    double transmission = 1.0;
    double optical_depth = 0.0;
};

struct TransferAdvanceResult {
    TransferError error = TransferError::None;
    BackwardTransferState state;
    std::string message;

    explicit operator bool() const noexcept;
};

TransferAdvanceResult advance_backward_transfer(
    const BackwardTransferState& current,
    const TransferCoefficients& coefficients,
    double positive_ds);

RedshiftResult evaluate_redshift(
    const Metric& metric,
    const PhaseSpaceState& photon,
    const Contravariant4& emitter_four_velocity,
    double observer_frequency);

double specific_intensity_at_observer(
    double invariant_intensity,
    double observer_frequency);
```

`advance_backward_transfer` is functional: an invalid input returns an error
and never partially mutates caller state. It requires:

- finite nonnegative `I`, `A`, `J`, and `ds`;
- finite `T` in `[0,1]`;
- nonnegative optical depth, where positive infinity is accepted only after
  physical attenuation saturation;
- finite output intensity.

Zero-length and vacuum steps are identities. Negative `ds` is rejected so the
backward sign cannot be silently applied twice.

## 5. Fluid and emission contracts

`include/solar/relativity/fluid_model.h` declares the master-prompt interface:

```cpp
struct FluidSample {
    bool valid = false;
    double density = 0.0;
    double temperature = 0.0;
    Contravariant4 four_velocity;
};

class FluidModel {
public:
    virtual ~FluidModel() = default;
    virtual FluidSample sample(
        const Metric& metric,
        const Contravariant4& x) const = 0;
};
```

`valid=false` means vacuum. Other fields are ignored in that case. A valid
sample must have finite nonnegative density and temperature and a finite unit
timelike four-velocity at the sampled point:

```text
g_mu_nu u^mu u^nu = -1
```

within `1e-10`. Model exceptions and malformed samples become explicit
`TransferError` results at the evaluation boundary.

The public model set is:

```cpp
class VacuumFluid final : public FluidModel;
class AnalyticCircularDiskFluid final : public FluidModel;
class AnalyticOpticallyThinTorus final : public FluidModel;
```

`include/solar/relativity/emission_model.h` declares:

```cpp
class EmissionModel {
public:
    virtual ~EmissionModel() = default;
    virtual TransferCoefficients coefficients(
        const FluidSample& fluid,
        const Covariant4& photon_p,
        double nu_observer) const = 0;
};

class VacuumEmission final : public EmissionModel;
class GreyEmission final : public EmissionModel;
class DebugPaintEmission final : public EmissionModel;
```

`evaluate_transfer_coefficients(...)` in `radiative_transfer.h` is the trusted
boundary that:

1. validates metric point and photon/sample finiteness;
2. treats `valid=false` as vacuum without reading unused sample fields;
3. validates the emitter four-velocity;
4. computes positive emitter frequency and `g`;
5. calls the emission model;
6. rejects negative, non-finite, or thrown coefficients.

`GreyEmission` takes nonnegative comoving coefficients per unit density:

```text
j_nu = emissivity_per_density * density
alpha_nu = absorption_per_density * density
J = j_nu / nu_emit^2
A = nu_emit alpha_nu
```

It is spectrally grey and ignores temperature by definition.

`VacuumEmission` always returns zero coefficients.

`DebugPaintEmission` returns configured invariant coefficients for every valid
sample. It is an explicitly diagnostic nonphysical model used to expose
sampling, sign, and composition bugs. It is never a scientific default and
does not define RGB.

## 6. Analytic Kerr matter models

Both analytic matter models are configured with finite subextremal
`mass_M`, `spin_chi`, and an `OrbitSense`. They support
`KerrBoyerLindquistMetric` and `KerrSchildCartesianMetric` whose parameters
match the model. An unsupported metric or parameter mismatch throws
`std::domain_error` rather than returning vacuum and hiding configuration
drift. The trusted evaluation boundary converts that exception to
`InvalidFluidSample`. A supported point outside the material's geometric
support returns vacuum.

An internal `kerr_fluid_kinematics` L1 module owns:

- BL/KS point conversion in the safe exterior overlap;
- the circular angular velocity
  `Omega=s/[M((r/M)^(3/2)+s chi)]`;
- normalization of `u^mu=u^t(1,0,0,Omega)` against the metric at the actual
  `(r,theta)`;
- contravariant velocity transformation `u_KS=J u_BL`;
- unit-timelike and parameter-consistency checks.

Here `s` is the existing coordinate rotation sign from `OrbitSense`. Away
from the equator this is a prescribed circular kinematic flow, not a claim of
geodesic or hydrostatic equilibrium.

### 6.1 Circular thin-disk fluid

`AnalyticCircularDiskConfig` contains:

- `mass_M`, `spin_chi`, and `OrbitSense`;
- optional inner radius, defaulting to the matching Kerr ISCO;
- finite outer radius greater than the inner radius;
- positive density and temperature scales;
- nonnegative density power;
- a positive dimensionless surface-height tolerance.

The sample is valid only in the configured radial interval and within the
equatorial surface tolerance. The controlled profile is:

```text
x = r / r_in
density = density_scale x^(-density_power)
flux_shape = x^(-3) max(0, 1-sqrt(1/x))
temperature = temperature_scale flux_shape^(1/4)
```

This is a Novikov–Thorne-like zero-torque shape, not a full Page–Thorne
spectral disk. It excludes magnetic stresses, returning radiation, disk
self-gravity, pressure support, and GRMHD evolution.

### 6.2 Optically thin torus

`AnalyticOpticallyThinTorusConfig` contains:

- matching Kerr parameters and `OrbitSense`;
- center radius, radial width, and angular width;
- positive density and temperature scales;
- nonnegative temperature power;
- a density cutoff fraction strictly between zero and one.

For BL coordinates:

```text
q_r = (r-r_center)/radial_width
q_theta = cos(theta)/angular_width
shape = exp[-0.5(q_r^2+q_theta^2)]
density = density_scale shape
temperature = temperature_scale shape^(temperature_power)
```

Samples below the cutoff are vacuum. The compact cutoff prevents an
effectively infinite atmosphere and makes integration bounds testable.
This is a kinematic Gaussian torus for transfer and renderer validation, not a
Fishbone–Moncrief equilibrium or a GRMHD snapshot.

BL and KS samples of the same event must agree in density and temperature
within `1e-12` relative error. Their four-velocities must agree after the
canonical chart Jacobian within `1e-10`, and each chart norm must be `-1`
within `1e-10`.

## 7. Thin-disk surface contract

Thin-disk surface emission is not returned as a volume
`TransferCoefficients` sample.

`include/solar/relativity/thin_disk.h` declares:

```cpp
enum class DiskOpacityMode {
    Opaque,
    SemiTransparent,
};

struct ThinDiskSurfaceEmissionSample {
    double emitted_specific_intensity;
    double emitted_bolometric_intensity;
    double surface_optical_depth;
};

class ThinDiskSurfaceEmission {
public:
    ThinDiskSurfaceEmission(
        double specific_intensity_scale,
        double bolometric_intensity_scale,
        double surface_optical_depth);

    ThinDiskSurfaceEmissionSample evaluate(
        const FluidSample& fluid) const;
};

struct ThinDiskCrossing {
    double affine;
    Contravariant4 position;
    double disk_radius;
    Contravariant4 emitter_four_velocity;
    double emitter_frequency;
    double redshift_g;
    Contravariant4 surface_normal;
    bool front_facing;
    std::size_t image_order;
    double observed_temperature;
    double observed_specific_intensity;
    double observed_bolometric_intensity;
};

struct ThinDiskObservedState {
    double specific_intensity = 0.0;
    double bolometric_intensity = 0.0;
    double transmission = 1.0;
};

struct ThinDiskRecorderConfig {
    DiskOpacityMode opacity_mode = DiskOpacityMode::Opaque;
    std::size_t max_crossings = 8;
};

struct ThinDiskRecordResult {
    TransferError error = TransferError::None;
    bool recorded = false;
    bool closed = false;
    std::string message;

    explicit operator bool() const noexcept;
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

The normalized surface model uses:

```text
I_emit_specific = specific_intensity_scale * temperature
I_emit_bolometric = bolometric_intensity_scale * temperature^4
T_obs = g T_emit
I_obs_specific = g^3 I_emit_specific
I_obs_bolometric = g^4 I_emit_bolometric
```

The scales define caller-selected intensity units. This validates the
relativistic transformations without falsely claiming an absolute spectrum.
A physical Planck spectrum can be added later behind the same surface
boundary once unit and bandpass contracts exist.

The north-facing unit surface normal is built in BL and transformed as a
contravariant vector to KS. `front_facing` is determined by the sign of
`p_mu n^mu` for the future photon momentum. The normal must be finite,
unit-spacelike, and orthogonal to the emitter velocity.

`ThinDiskCrossingRecorder` owns a bounded vector with a required maximum
between 1 and a conservative implementation ceiling. Default maximum is 8.
It receives already localized surface states in observer-to-past order.

- invalid/out-of-radius surface samples are ignored explicitly;
- valid crossings record every field listed by the master prompt;
- image order starts at zero and is monotonic;
- opaque mode accepts the first valid crossing, applies its full source
  intensity, sets transmission to zero, and closes the recorder;
- semi-transparent mode applies the exact sheet solution
  `T source (1-exp(-tau_surface))`, then
  `T<-T exp(-tau_surface)`;
- intensity composition is scalar linear intensity, never sRGB alpha blend;
- exceeding the configured crossing count returns
  `CrossingLimitReached` without reallocating beyond the bound.

Geodesic event localization remains owned by the existing event subsystem.
Phase 6 supplies disk-surface events and accepted crossings to this recorder.
This avoids a second geodesic integrator and avoids pretending an arbitrary
surface root is a volume coefficient.

## 8. Failure and numerical semantics

Constructor/configuration errors throw `std::invalid_argument`, matching
existing Solar configuration conventions. Runtime samples from metrics,
fluids, and emission models are untrusted and return explicit
`TransferError` with a message.

The following never receive a fallback value:

- non-finite/negative density, temperature, emissivity, absorption, or step;
- invalid metric point;
- non-unit or non-timelike emitter velocity;
- non-positive emitter frequency;
- mismatched Kerr model/metric parameters;
- non-finite redshift or transformed velocity;
- non-finite accumulated intensity;
- malformed surface optical depth or crossing state.

Vacuum is not an error. It is represented only by `FluidSample.valid=false`
or zero coefficients.

## 9. Validation strategy

### 9.1 Exact flat-space transfer

Use constant coefficients in Minkowski space and compare:

```text
J=2.5, A=0, ds=4       => I=10, T=1
J=0, A=0.7, ds=3       => background attenuation exp(-2.1)
J=3, A=0.4, ds=5       => source contribution (J/A)(1-exp(-2))
```

Split each interval into 1, 2, 7, and 100 constant substeps. Results must
agree with the one-step analytic solution within `5e-14`. Optically thin
`A ds=1e-12` and thick `A ds=1000` cases lock the `expm1` and saturation
branches. A deliberate `exp(+tau)` mutation must fail.

### 9.2 Frequency and coefficient gates

Minkowski boosted emitters have analytic Doppler frequency. Require
`nu_emit`, `g`, grey `J`, and grey `A` within `2e-14`. Reject past-directed
photons, spacelike/non-unit emitters, invalid fluid fields, thrown emission
models, negative coefficients, non-finite observer frequency, and negative
`ds`.

### 9.3 Fluid models

Validate disk ISCO default, inner/outer/surface refusal, profile values,
negative spin, both orbit senses, torus center/cutoff/symmetry, finite
positive samples, and timelike normalization. Cross-chart fixtures require
the BL/KS gates from Section 6.

### 9.4 Surface transfer

Synthetic exact crossings validate:

- one crossing and all recorded fields;
- two and eight ordered crossings;
- maximum-count refusal;
- opaque first-hit closure;
- semi-transparent two-sheet analytic composition;
- front/back classification;
- `g^3` specific intensity;
- `g^4` bolometric intensity;
- `g` temperature shift;
- BL/KS normal and redshift agreement.

Mutation checks flip the backward attenuation sign and replace the canonical
emitter velocity transform with a spatial copy; the focused tests must fail.

### 9.5 Release gate

Run:

```bash
make clean
make -j4 test
make test-external-consumer
git diff --check
```

Then build all Phase 5 focused executables once with combined ASan/UBSan and
run every assertion. Record exact totals, compiler, skips, maximum errors,
installed-consumer output, and mutation evidence in
`docs/validation/relativity_10_radiative_transfer.md`.

The optional missing `data/de440.asc` fixture may remain a visible skip.
No required transfer, fluid, surface, installation, sanitizer, or Linux CI
gate may be skipped.

## 10. External references

The invariant equation and constant-coefficient formal solution follow:

- Bronzwaer et al., RAPTOR I, arbitrary-spacetime GR radiative transfer:
  <https://arxiv.org/abs/1801.10452>
- Mościbrodzka and Gammie, `ipole`, analytic constant-coefficient covariant
  transport:
  <https://arxiv.org/abs/1712.03057>
- Event Horizon Telescope Collaboration, radiative-transfer code
  verification:
  <https://doi.org/10.3847/1538-4357/ab96c6>

External implementations are references only. No code is copied. Solar keeps
its existing signature, chart, affine, momentum-variance, and unit contracts.

## 11. Expected file ownership

Public:

```text
include/solar/relativity/radiative_transfer.h
include/solar/relativity/fluid_model.h
include/solar/relativity/emission_model.h
include/solar/relativity/thin_disk.h
```

Implementation:

```text
src/relativity/radiative_transfer.cpp
src/relativity/transfer_evaluation.cpp
src/relativity/emission_model.cpp
src/relativity/kerr_fluid_kinematics.h
src/relativity/kerr_fluid_kinematics.cpp
src/relativity/analytic_disk_fluid.cpp
src/relativity/analytic_torus_fluid.cpp
src/relativity/thin_disk.cpp
src/relativity/thin_disk_surface.cpp
src/relativity/thin_disk_geometry.h
src/relativity/thin_disk_geometry.cpp
```

Focused tests:

```text
tests/relativity/test_radiative_transfer.cpp
tests/relativity/test_emission_models.cpp
tests/relativity/test_fluid_models.cpp
tests/relativity/test_thin_disk.cpp
tests/relativity/test_transfer_failures.cpp
```

The split follows scalar transfer, runtime evaluation, matter kinematics, and
surface composition boundaries. It does not introduce one-line wrappers or a
new framework.
