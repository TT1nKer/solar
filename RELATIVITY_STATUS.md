# Solar Relativity Status

CURRENT_PHASE: 1
PHASE_STATE: PASSED
LAST_VERIFIED_COMMIT: 57b65256094e6d6b41c32ffd2906310b0ed4b134
LAST_VERIFIED_PLATFORM: Darwin 23.6.0 arm64 / Apple Clang 16.0.0
LAST_VERIFIED_COMMANDS:
- `make clean`
- `make`
- `make test`
- all 15 `tests/relativity/test_*` executables
- AddressSanitizer and UndefinedBehaviorSanitizer build plus all relativity tests
- DOPRI, Hamiltonian, RMS, dense-output, and invalid-domain mutation checks
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
- Made E and Lz monitoring opt-in and left unavailable Carter diagnostics as NaN.
- Preserved the existing dynamic-vector generic integrator output and the specialized N-body interface.
- Validated Minkowski null/timelike lines, reversal, Schwarzschild radial null motion, photon sphere, weak bending, and ordinary Kerr E/Lz behavior.
- Passed 386/386 relativity assertions and 66/66 fixture-independent legacy assertions in release mode.
- Passed all 386 relativity assertions under AddressSanitizer and UndefinedBehaviorSanitizer.

NOT_COMPLETED:
- Observer/tetrad construction, local photon/timelike initialization, measured-frequency normalization, or general future-direction classification.
- Carter constant, separated Kerr/Mino-time solver, or long-time structure-preserving timelike integration.
- Kerr–Schild coordinates, reliable physical horizon crossing, interior evolution, or singularity treatment.
- Disk/material models, radiative transfer, reference renderer, image/movie pipeline, Solar adapter, WASM/GPU, UI, or visual regression.
- DE440's eight external-data assertions on this machine.
- GCC/Linux or non-arm64 verification.

CURRENT_BLOCKERS:
- None.

MOST_LIKELY_BUGS:
- Endpoint-bracket event detection can miss tangencies, multiple roots, or an even number of roots inside one accepted step.
- Boyer–Lindquist invalid-domain termination is not physical horizon capture; Phase 3 Kerr–Schild validation is still required.
- Long bound timelike trajectories can accumulate secular error because Phase 1 has no structure-preserving integrator.
- Caller-supplied coordinate states can have an unintended local direction until Phase 2 tetrad initialization exists.
- Near-extremal, near-margin, extreme-momentum, and long-duration cases lack a multiprecision reference sweep.
- E/Lz monitoring is caller-enabled and symmetry-dependent; Carter is intentionally unavailable.
- Only Apple Clang 16 on macOS arm64 was verified locally.
- Missing default DE440 data remains a visible skip.

FASTEST_WAY_TO_FALSIFY:
- `./tests/relativity/test_dopri5`: tableau, RMS, per-component scale, controller, dense output, and finite-value gates.
- `./tests/relativity/test_hamiltonian`: Hamilton equations, variance/order, constraint denominator, and no-false-symmetry path.
- `./tests/relativity/test_geodesic_events`: direction, endpoint, negative-step, bracket, and iteration semantics.
- `./tests/relativity/test_geodesic_failures`: invalid trial points must shrink/reject and never become `HorizonCrossing`.
- `./tests/relativity/test_geodesics`: analytic Minkowski lines, limits, first event, and reversal.
- `./tests/relativity/test_geodesics_schwarzschild`: radial null, photon sphere, weak bending, and `1e-10` constraint gate.
- `./tests/relativity/test_geodesics_kerr`: ordinary Kerr constraint, exact monitored E/Lz, and unavailable Carter.
- Rebuild all relativity tests with ASan/UBSan; any runtime diagnostic invalidates the gate.

NEXT_ALLOWED_ACTION:
- Phase 2 only: implement observer and tetrad construction, local initialization, frequency normalization, and round-trip validation. Do not enter Phase 3.
