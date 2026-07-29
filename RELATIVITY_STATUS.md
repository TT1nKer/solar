# Solar Relativity Status

CURRENT_PHASE: 1
PHASE_STATE: PASSED
LAST_VERIFIED_COMMIT: 1bf7442c141e5440642af18677dc8680f96bfc9d
LAST_VERIFIED_PLATFORM: Darwin 23.6.0 arm64 / Apple Clang 16.0.0
LAST_VERIFIED_COMMANDS:
- `make clean`
- `make`
- `make test`
- all 15 `tests/relativity/test_*` executables
- AddressSanitizer and UndefinedBehaviorSanitizer build plus all relativity tests
- DOPRI, Hamiltonian, RMS, dense-output, and invalid-domain mutation checks
- independent Kerr inverse-domain, Carter-Q, Schwarzschild convergence, and
  high-precision bending probes
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
- Made E and Lz monitoring opt-in and left unavailable Carter diagnostics as NaN.
- Preserved the existing dynamic-vector generic integrator output and the specialized N-body interface.
- Validated Minkowski null/timelike lines, reversal, Schwarzschild radial null motion, photon sphere, weak bending, and ordinary Kerr E/Lz behavior.
- Passed 403/403 relativity assertions and 67/67 fixture-independent legacy assertions in release mode.
- Passed all 403 relativity assertions plus the 11-assertion legacy DOPRI adapter test under AddressSanitizer and UndefinedBehaviorSanitizer.

NOT_COMPLETED:
- Observer/tetrad construction, local photon/timelike initialization, measured-frequency normalization, or general future-direction classification.
- Carter constant, separated Kerr/Mino-time solver, or long-time structure-preserving timelike integration.
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
- Caller-supplied coordinate states can have an unintended local direction until Phase 2 tetrad initialization exists.
- The universal tolerance factory follows v3 component defaults but is not
  chart-aware; non-unit mass scales with angular BL coordinates need a
  dedicated convergence sweep before scientific use.
- Near-extremal, near-margin, extreme-momentum, and long-duration cases lack a multiprecision reference sweep.
- E/Lz monitoring is caller-enabled and symmetry-dependent; Carter is intentionally unavailable.
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
- `./tests/relativity/test_geodesics_kerr`: ordinary Kerr constraint, exact monitored E/Lz, and unavailable Carter.
- Rebuild all relativity tests with ASan/UBSan; any runtime diagnostic invalidates the gate.

NEXT_ALLOWED_ACTION:
- Phase 2 only: implement observer and tetrad construction, local initialization, frequency normalization, and round-trip validation. Do not enter Phase 3.
