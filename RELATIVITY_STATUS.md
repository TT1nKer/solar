# Solar Relativity Status

CURRENT_PHASE: 2
PHASE_STATE: PASSED
LAST_VERIFIED_COMMIT: 8a6a2533972685f31dea8eae5f361f948546f285
LAST_VERIFIED_PLATFORM: Darwin 23.6.0 arm64 / Apple Clang 16.0.0
LAST_VERIFIED_COMMANDS:
- `make clean`
- `make`
- `make test`
- all 22 `tests/relativity/test_*` executables
- AddressSanitizer and UndefinedBehaviorSanitizer build plus all 22 relativity
  tests and `tests/test_integrator`
- DOPRI, Hamiltonian, RMS, dense-output, invalid-domain, observer, local-state,
  Carter, special-orbit, and analytic/numerical-shadow checks
- horizontal-screen-sign, Bardeen-`xi`-sign, and
  invalid-metric-as-capture mutation checks
- independent 80-decimal Kerr special-radius and Bardeen-edge probe
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
- Added the asymptotic Bardeen critical curve, a stable small-spin
  Schwarzschild branch, and finite-domain filtering without fabricated roots.
- Added an independent distant-ZAMO CPU backward-ray benchmark. It keeps
  momentum future-directed, uses negative affine integration, distinguishes BL
  chart failure from capture, and recovers both horizontal Kerr shadow edges.
- Passed 1880/1880 relativity assertions across 22 executables and 67/67
  fixture-independent legacy assertions in release mode.
- Passed all 1880 relativity assertions plus the 11-assertion legacy DOPRI
  adapter test under AddressSanitizer and UndefinedBehaviorSanitizer.

NOT_COMPLETED:
- Separated Kerr radial/polar potentials, Mino-time solver, fundamental
  frequencies, turning-point handling, or long-time structure-preserving
  timelike integration.
- Kerr–Schild coordinates, reliable physical horizon crossing, interior evolution, or singularity treatment.
- Disk/material models, radiative transfer, reference renderer, image/movie pipeline, Solar adapter, WASM/GPU, UI, or visual regression.
- DE440's eight external-data assertions on this machine.
- GCC/Linux or non-arm64 verification.
- The one-off Carter-Q, Kerr inverse-domain, convergence, and multiprecision
  audit probes are supplementary evidence, not committed CI regression
  executables.

CURRENT_BLOCKERS:
- None.

MOST_LIKELY_BUGS:
- Endpoint-bracket event detection can miss tangencies, multiple roots, or an even number of roots inside one accepted step.
- Boyer–Lindquist invalid-domain termination is not physical horizon capture; Phase 4 Kerr–Schild validation is still required.
- The Kerr BL precision boundary is intentionally stricter than the geometric
  exterior and can vary slightly with floating-point platform; callers must
  treat rejection as chart/numerical invalidity, not capture.
- Long bound timelike trajectories can accumulate secular error because Phase 1 has no structure-preserving integrator.
- Near-degenerate arbitrary-observer seeds can cross the fixed numerical
  rejection threshold differently on other floating-point platforms.
- The universal tolerance factory follows v3 component defaults but is not
  chart-aware; non-unit mass scales with angular BL coordinates need a
  dedicated convergence sweep before scientific use.
- Near-extremal, near-margin, extreme-momentum, and long-duration cases lack a multiprecision reference sweep.
- Carter monitoring is callback-based and the validated generic trajectory is
  short; long bound-orbit secular drift is not characterized.
- The analytic off-equatorial Bardeen curve samples the spherical-photon
  interval rather than solving visible-branch endpoint roots adaptively, so
  coarse sample counts can under-resolve the tips.
- The CPU shadow benchmark checks only two equatorial horizontal edges at
  `r=1000M`, `chi=0.5`; it does not validate a full 2D image or near-extremal
  convergence.
- Circular-orbit formulas are double precision and have not received a dense
  near-extremal sweep against an independent package.
- Only Apple Clang 16 on macOS arm64 was verified locally.
- Missing default DE440 data remains a visible skip.

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
  reflection, mass scaling, small-spin stability, and invalid inputs.
- `./tests/relativity/test_kerr_shadow_raytrace`: future-directed
  negative-affine rays, explicit capture/escape classification, both Bardeen
  edge comparisons, and Hamiltonian/Carter gates.
- `./tests/relativity/test_geodesics_kerr`: ordinary Kerr constraint, exact
  monitored E/Lz, Carter drift, denominator semantics, and callback failures.
- Rebuild all relativity tests with ASan/UBSan; any runtime diagnostic invalidates the gate.

NEXT_ALLOWED_ACTION:
- Phase 3 only: implement separated Kerr radial/polar potentials, turning
  points, and Mino-time validation. Do not enter Phase 4 Kerr–Schild, matter,
  transfer, renderer, GPU, or UI work before the Phase 3 gate passes.
