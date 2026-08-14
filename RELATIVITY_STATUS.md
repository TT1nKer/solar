# Solar Relativity Status

CURRENT_PHASE: 5
PHASE_STATE: PASSED
LAST_VERIFIED_COMMIT: a65befbc1ae1462b52297b2a2e5c9e82e040eb3c
LAST_VERIFIED_PLATFORM: Darwin 23.6.0 arm64 / Apple Clang 16.0.0;
  Ubuntu 24.04 x86_64 / GCC 13.3.0 via GitHub C++ CI
LAST_VERIFIED_COMMANDS:
- `make clean`
- `make -j4 test`
- `make test-external-consumer`
- all 39 `tests/relativity/test_*` executables
- combined AddressSanitizer/UndefinedBehaviorSanitizer build and all five
  Phase 5 focused executables
- wrong-sign observer-to-past attenuation mutation
- BL-to-KS fluid four-velocity component-copy mutation
- GitHub C++ CI build, complete test suite, installed external consumer, CLI,
  and sample-command gates
- `git diff --check`

ACTUALLY_COMPLETED:
- Preserved the Phase -1/0A/0B convention, tensor, unit, derivative, and metric foundation.
- Added canonical Hamilton equations for general fixed-background metrics with coordinate order `x^mu,p_mu`.
- Added the exact v3 normalized null/timelike Hamiltonian constraint and finite-scale rejection.
- Added a shared fixed/vector container DOPRI5 stage engine with per-component tolerances, stable RMS or legacy maximum norm, adaptive rejection, and sign-preserving negative steps.
- Added fourth-order dense output with closed-interval evaluation and no extrapolation.
- Added directed bracketed event roots with safeguarded secant/bisection, endpoint roots, iteration failure, and first-event selection.
- Added explicit v3 diagnostics and termination reasons, including `InteriorCutoff`, without claiming it is executable before Kerr–Schild work.
- Added an allocation-bounded L2 geodesic flow with affine/proper/coordinate/step limits, invalid-trial shrink, rejection recovery, constraint gates, and no default projection.
- Rejects DOPRI5 steps that cannot advance the floating-point independent variable, so canonical state cannot advance while affine remains frozen.
- Validates geodesic/event enums and every event contract before trajectory work instead of silently choosing fallback physics or hiding malformed input behind a later metric failure.
- Restricts Kerr BL `valid_point()` to the v3 `1e-10` near-horizon inverse-identity precision gate; ill-conditioned exterior points now terminate explicitly and remain distinct from horizon capture.
- Kept E, Lz, and Carter monitoring opt-in; the generic integrator receives
  Carter through an explicit invariant evaluator and reports relative and
  absolute drift without depending on Kerr.
- Preserved the existing dynamic-vector generic integrator output and the specialized N-body interface.
- Validated Minkowski null/timelike lines, reversal, Schwarzschild radial null motion, photon sphere, weak bending, and ordinary Kerr E/Lz behavior.
- Added metric contraction, index raising/lowering, and covector/vector pairing
  primitives with fixed variance.
- Added static, arbitrary, look-at, ZAMO, and equatorial circular observers
  with Lorentzian Gram–Schmidt, right-handed tetrads, full 16-component frame
  validation, and explicit nonexistence/failure outcomes.
- Added local future-directed photon and subluminal timelike initialization,
  observer-frequency measurement, photon normalization, and Hamiltonian gates.
- Added Kerr `E`, `Lz`, and Carter evaluation plus generic opt-in Carter drift
  monitoring and explicit callback failure semantics.
- Added analytic Schwarzschild/Kerr ISCO, equatorial photon, marginally bound,
  and circular timelike quantities with spin-relative orbit sense.
- Added the asymptotic Bardeen critical curve, inclination-dependent visible
  tip solves, cancellation-safe screen-alpha sampling, a stable small-spin
  Schwarzschild branch, and finite-domain filtering without fabricated roots.
- Added an independent distant-ZAMO CPU backward-ray benchmark. It keeps
  momentum future-directed, uses negative affine integration, distinguishes BL
  chart failure from capture, enforces `alpha=-Lz/E`, recovers both horizontal
  Kerr shadow edges at `r=1000M` and `2000M`, and verifies the Schwarzschild
  critical radius.
- Added literal Kerr separated radial and `mu=cos(theta)` potentials with
  analytic derivatives, normalized residual scales, null/timelike support,
  and Hamiltonian/canonical state round trips.
- Added the public `KerrSeparatedIntegrator` Mino-time CPU path with fixed-size
  DOPRI5 state, affine/proper limits, public events, negative integration,
  explicit invalid-domain failures, and no per-step allocation.
- Added rejected forbidden-potential trials, bracketed machine-refined simple
  roots, controlled radial/polar release, exact turn counters, and explicit
  `NearCriticalOrbit` refusal for double/critical roots.
- Added a synchronized smooth turning-phase state that preserves
  `dr/dgamma,dmu/dgamma` through the square-root endpoint, localizes later
  velocity crossings, and refuses normalized phase drift above the `1e-10`
  CPU constraint gate.
- Localizes public termination events inside the nonzero Mino interval used
  to approach and release a turning root, so affine and coordinate limits
  cannot be skipped by the turn transition.
- Added render-relevant minimum radius, azimuth/winding, radial/polar residual,
  Hamiltonian-constraint, Carter-drift, accepted/rejected-step, and
  termination diagnostics.
- Cross-validated six common-event families against the generic Hamiltonian
  solver: null/timelike, both spin signs, Schwarzschild, radial return, and
  polar return. Release P95/maximum worldline error is `2.27726e-11`.
- Proved the separated API through a clean installed external consumer.
- Added Cartesian ingoing Kerr–Schild covariant/inverse metrics from the
  null-form identity, a cancellation-safe positive implicit radius, analytic
  radius gradient, and `Dual4` inverse-metric derivatives.
- Added a finite exterior Boyer–Lindquist/Kerr–Schild overlap transform with
  full position Jacobians, inverse maps, angle wrapping, explicit axis/margin
  refusal, and inverse-transpose canonical momentum conversion.
- Added chart-aware stationary-energy and axial-angular-momentum callbacks to
  the generic invariant monitor while preserving the existing Boyer–Lindquist
  defaults and aggregate initialization compatibility.
- Added explicit decreasing-radius outer-horizon and interior-cutoff events.
  A localized timelike plunge restarts at the exact horizon event state and
  evolves a positive affine interval to `max(0.05M, configured)`.
- Cross-validated ordinary and near-horizon Kerr–Schild Hamiltonian flows
  against Boyer–Lindquist at common exterior events. Ordinary position and
  momentum P95 errors are `6.64417e-12` and `6.24937e-12`; near-horizon errors
  are `6.13639e-10` and `2.47107e-10`.
- Proved the installed public API through a clean consumer that constructs the
  Kerr–Schild metric, transforms canonical state, integrates, and localizes an
  event.
- Passed 2267/2267 relativity assertions across 34 executables and 67/67
  fixture-independent legacy assertions in release mode.
- Passed all 189 assertions in the seven Phase 4 focused executables under
  AddressSanitizer and UndefinedBehaviorSanitizer.
- Added a stable observer-to-past invariant formal solution with
  `expm1`/small-optical-depth handling, exact constant-segment absorption, and
  explicit finite-state failures.
- Added validated local emitter frequency and redshift plus vacuum, grey, and
  diagnostic emission models behind public fluid/emission boundaries.
- Added analytic Kerr circular-disk and compact kinematic Gaussian-torus
  material models with full Boyer–Lindquist/Kerr–Schild four-vector
  transformation and explicit metric-parameter mismatch failures.
- Added bounded thin-disk surface records with unit north normal, `g`, `g^3`,
  and `g^4` transformations, front/back classification, opaque closure,
  semi-transparent linear-intensity composition, and an explicit crossing
  limit.
- Proved Phase 5 public headers and symbols through a clean installed external
  consumer that executes constant transfer, disk/torus samples, and a surface
  crossing.
- Passed 2462/2462 relativity assertions across 39 executables and 67/67
  fixture-independent legacy assertions in Release mode.
- Passed all 195 assertions in the five Phase 5 focused executables under
  AddressSanitizer and UndefinedBehaviorSanitizer.

NOT_COMPLETED:
- Analytic elliptic Kerr geodesics, fundamental frequencies, or long-time
  structure-preserving timelike integration.
- Extremal/negative-radius Kerr extension, ring-singularity treatment, or
  quantum-gravity interior model.
- CPU reference renderer, thin-disk event integration, reference images,
  movie pipeline, Solar adapter, WASM/GPU, UI, or visual regression.
- Absolute spectral calibration, detector bandpass, polarization, scattering,
  returning radiation, self-gravity, magnetic evolution, or GRMHD.
- DE440's eight external-data assertions on this machine.
- Linux sanitizer or non-x86_64/non-arm64 verification.
- The one-off Carter-Q, Kerr inverse-domain, convergence, and multiprecision
  audit probes are supplementary evidence, not committed CI regression
  executables.

CURRENT_BLOCKERS:
- None.

MOST_LIKELY_BUGS:
- Endpoint-bracket event detection can miss tangencies, multiple roots, or an even number of roots inside one accepted step.
- Boyer–Lindquist invalid-domain termination remains chart failure, not
  physical capture; callers requiring horizon crossing must use the validated
  Cartesian Kerr–Schild path.
- The Kerr BL precision boundary is intentionally stricter than the geometric
  exterior and can vary slightly with floating-point platform; callers must
  treat rejection as chart/numerical invalidity, not capture.
- Long bound timelike trajectories can accumulate secular error because Phase 1 has no structure-preserving integrator.
- Near-degenerate arbitrary-observer seeds can cross the fixed numerical
  rejection threshold differently on other floating-point platforms.
- Kerr–Schild supports positive-radius subextremal Kerr only. The BL/KS
  transform intentionally refuses the polar axis and `r<=r_++margin`; callers
  must initialize in a safe exterior overlap.
- Interior evolution stops at the explicit numerical/model boundary
  `max(0.05M, configured)` and does not represent the ring singularity.
- The universal tolerance factory follows v3 component defaults but is not
  chart-aware; non-unit mass scales with angular BL coordinates need a
  dedicated convergence sweep before scientific use.
- Near-extremal, near-margin, extreme-momentum, and long-duration cases lack a multiprecision reference sweep.
- Carter monitoring is callback-based and the validated generic trajectory is
  short; long bound-orbit secular drift is not characterized.
- Sub-ULP near-axis Bardeen intervals may be rejected explicitly when the two
  physical tips cannot be distinguished in the platform floating-point type.
- The CPU shadow benchmark checks only two equatorial horizontal edges at
  `r=1000M` and `2000M`, `chi=0.5`, plus the Schwarzschild critical radius; it
  does not validate a full 2D image or near-extremal convergence.
- Circular-orbit formulas are double precision and have not received a dense
  near-extremal sweep against an independent package.
- The separated cross-check covers six moderate fixtures and one simple turn
  of each type; long bound, repeated-turn, and near-extremal scientific sweeps
  remain uncharacterized.
- Near-axis separated trajectories with nonzero `Lz` are rejected because the
  Boyer–Lindquist azimuth term is coordinate-singular.
- The smooth turn-phase invariant ceiling is `1e-10`; a trajectory exceeding
  it fails explicitly, but denser multiprecision characterization is pending.
- Sanitizer validation is limited to Apple Clang 16 on macOS arm64; the
  Ubuntu/GCC CI run covers Release behavior only.
- Missing default DE440 data remains a visible skip.
- The disk temperature profile is a controlled zero-torque approximation and
  the torus is a kinematic Gaussian; neither is a calibrated plasma model.
- Material, frequency, and intensity units remain caller-selected. No absolute
  spectrum or detector response has been validated.
- Phase 5 consumes caller-supplied constant-coefficient segments and localized
  surface states; Phase 6 must validate adaptive accepted-step sampling and
  repeated disk-event localization.
- Near-extremal, long-path, and dense multi-crossing transfer lack a broader
  independent precision sweep.

FASTEST_WAY_TO_FALSIFY:
- `./tests/relativity/test_dopri5`: tableau, RMS, per-component scale, controller, dense output, finite-value gates, and representable independent-variable progress.
- `./tests/relativity/test_hamiltonian`: Hamilton equations, variance/order, constraint denominator, and no-false-symmetry path.
- `./tests/relativity/test_geodesic_events`: direction, unknown-enum rejection, endpoint, negative-step, bracket, and iteration semantics.
- `./tests/relativity/test_geodesic_failures`: malformed event contracts fail before integration; invalid trial points must shrink/reject and never become `HorizonCrossing`.
- `./tests/relativity/test_geodesics`: analytic Minkowski lines, limits, first event, reversal, enum validation, and affine-resolution underflow.
- `./tests/relativity/test_kerr_bl`: accepted near-horizon BL points must satisfy the `1e-10` inverse-identity gate.
- `./tests/relativity/test_geodesics_schwarzschild`: radial null, photon sphere, weak bending, and `1e-10` constraint gate.
- `./tests/relativity/test_observers`: static/arbitrary/look-at/ZAMO
  normalization, handedness, round trips, ergosphere and invalid-domain
  failures.
- `./tests/relativity/test_local_initialization`: photon/timelike constraints,
  observer frequency, future direction, negative-affine semantics, and invalid
  local inputs.
- `./tests/relativity/test_kerr_constants`: literal null/timelike Carter
  values and every evaluator failure boundary.
- `./tests/relativity/test_kerr_orbits`: Schwarzschild/Kerr special radii,
  stability/existence, signed spin, and lowered circular-observer invariants.
- `./tests/relativity/test_kerr_shadow`: Schwarzschild limit, Kerr endpoints,
  off-equatorial/near-axis interiors, reflection, mass scaling, small-spin
  stability, overflow, and invalid inputs.
- `./tests/relativity/test_kerr_shadow_raytrace`: future-directed
  negative-affine rays, exact conserved screen mapping, explicit
  capture/escape classification, two-radius Kerr and Schwarzschild boundary
  comparisons, and Hamiltonian/Carter gates.
- `./tests/relativity/test_geodesics_kerr`: ordinary Kerr constraint, exact
  monitored E/Lz, Carter drift, denominator semantics, and callback failures.
- `./tests/relativity/test_geodesic_invariant_callbacks`: chart-aware E/Lz
  evaluators, disabled-monitor behavior, legacy defaults, and callback failure
  semantics.
- `./tests/relativity/test_kerr_schild`: radius quartic, null one-form,
  analytic metric inverse, horizon/interior domain, Cartesian invariants, and
  invalid parameters.
- `./tests/relativity/test_kerr_schild_derivatives`: independent radius
  gradient and five-point inverse-derivative comparisons.
- `./tests/relativity/test_kerr_chart_transform`: position/Jacobian/canonical
  momentum round trips, pairing/Hamiltonian invariance, differential signs,
  overlap/axis refusal, and both spin signs.
- `./tests/relativity/test_kerr_bl_ks_crosscheck`: ordinary and near-horizon
  common-event worldlines, constraints, E/Lz, and convergence.
- `./tests/relativity/test_geodesics_kerr_schild`: exact horizon localization,
  tetrad regularity, restart continuity, positive interior evolution, cutoff
  localization, and timelike constraints.
- `./tests/relativity/test_kerr_separated_potentials`: literal `R/U`,
  derivatives, scaling, and invalid domains.
- `./tests/relativity/test_kerr_separated_state`: Hamilton/Mino direction,
  locked-root, critical-root, and canonical round trips.
- `./tests/relativity/test_kerr_separated`: limits, events, negative Mino
  direction, radial/polar turns, critical refusal, and diagnostics.
- `./tests/relativity/test_kerr_separated_crosscheck`: six common-event
  Hamiltonian comparisons, convergence, constraints, Carter drift, and
  P95/maximum worldline gates.
- `./tests/relativity/test_kerr_separated_turning_phase`: phase crossing,
  root-time accuracy, duplicate-count prevention, and critical classification.
- `make test-external-consumer`: installed headers, symbols, observer/local
  initialization, separated integration, invariant transfer, Kerr matter, and
  thin-disk surface composition.
- `./tests/relativity/test_radiative_transfer`: constant analytic transfer,
  stable thin/thick limits, subdivision, sign, and failure gates.
- `./tests/relativity/test_emission_models`: boosted redshift and invariant
  grey/vacuum/debug coefficient contracts.
- `./tests/relativity/test_transfer_failures`: malformed runtime samples,
  callback exceptions, non-timelike emitters, and finite-range failures.
- `./tests/relativity/test_fluid_models`: disk/torus profiles, support,
  timelike normalization, parameter mismatch, and BL/KS vector agreement.
- `./tests/relativity/test_thin_disk`: one/two/eight crossings, opacity,
  `g` powers, BL/KS normal, bounds, and failure atomicity.
- Rebuild all relativity tests with ASan/UBSan; any runtime diagnostic invalidates the gate.

NEXT_ALLOWED_ACTION:
- Phase 6 CPU reference rendering backed by the installed Phase 5 transfer
  API. Do not begin GPU, WASM, or UI before the CPU reference and regression
  image gates pass.

## Dynamic Collapse Track (branch codex/dynamic-collapse-bh)

Parallel track for the star-formation / black-hole-formation collapse
pipeline: Newtonian free-fall -> post-Newtonian corrections -> exact GR
dust collapse, with an analytic anchor at every stage.

- Milestone 01 (3D post-Newtonian collapse) — DONE:
  - Regime A: Barnes-Hut octree gravity (theta-ladder verified), free-fall
    cycloid anchor (<1% surface error), turbulent-cloud initial conditions
    (Kritsuk P(k) ~ k^-2, log-normal density, rotation; FFT-verified).
  - Regime B: field-based 1PN force (a = g[1 + 4 Phi/c^2 - v^2/c^2] -
    4(g.v)v/c^2 via the shared monopole walk) with c->infinity reduction
    at 1e-15, radial spot check at machine precision, perihelion
    precession 0.0100596 vs 0.0103569 rad/orbit; spherical-limit collapse
    regression (tracks the radial 1PN shell model to 0.94%, signed
    coordinate-time lag 0.135% vs 0.167% analytic); per-particle
    blending driver with compactness diagnostics and hand-off tagging.
  - Executables: test_barnes_hut_gravity, test_collapse_freefall,
    test_turbulent_cloud, test_turbulent_spectrum, test_pn_gravity,
    test_collapse_pn_spherical, test_collapse_blend (all passing).
- Milestone 02 (LTB compact stage) — DONE:
  - LTBCollapse: per-shell LTB dust trajectories, singularity/horizon
    times, shell-crossing monitor; OS surface observables (t_obs,
    redshift, (1+z)^-2 luminosity tail).
  - OS anchors: cycloid vs independent RK4 at 6.9e-12; simultaneous
    crunch; horizon at the analytic collapse angle; t_obs log-divergence
    coefficient 3 G M / c^3 to 0.3%; luminosity tail slope to 0.15%.
  - N-body -> LTB hand-off: initial-radius recovery 1.7%, exact passage
    through the hand-off state, horizon-then-singularity completion,
    per-shell profile recovery 2.1%.
  - Executables: test_ltb_os, test_ltb_handoff (passing).
- Next: milestone 03 (nebula-scale multi-mass/multi-component collapse
  on the turbulent ICs with the full blending + hand-off driver), then
  the traceable collapse video pipeline.
