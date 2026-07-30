# Solar Relativity Status

CURRENT_PHASE: 0B
PHASE_STATE: PASSED
LAST_VERIFIED_COMMIT: a234b160890d127cc546bbecc4d2addb97410a62
LAST_VERIFIED_PLATFORM: macOS 14.8.7 arm64 / Apple Clang 16.0.0
LAST_VERIFIED_COMMANDS:
- `make clean`
- `make`
- `make test`
- all eight `tests/relativity/test_*` executables
- AddressSanitizer and UndefinedBehaviorSanitizer build plus all relativity tests
- real Kerr metric CLI with JSON parsing
- `git diff --check`

ACTUALLY_COMPLETED:
- Preserved the Phase -1/0A build, convention, math, unit, tensor, and Dual4 foundation.
- Added the exact v3 `Metric` contract and chart identity.
- Added independently implemented Minkowski Cartesian and Schwarzschild BL metrics.
- Added subextremal Kerr BL covariant/inverse matrices, analytic covariant derivatives, inverse derivatives, horizons, and stationary limit.
- Stabilized Kerr surface calculations for large finite masses.
- Rejected horizon, axis, non-finite, overflowing, underflowing, and otherwise unrepresentable metric domains without clamping.
- Cross-checked Kerr analytic derivatives with Dual4 and inverse derivatives with five-point finite differences.
- Cross-checked Kerr `chi=0` against Schwarzschild, inverse identities, symmetry, signature, and literal surface values.
- Added the real `solar relativity metric` human/JSON CLI with strict parsing and truthful nonzero errors.
- Made CLI translation-unit discovery recursive while preserving legacy command routing.
- Passed 202/202 relativity assertions and 56/56 fixture-independent legacy assertions in release mode.
- Passed all 202 relativity assertions under AddressSanitizer and UndefinedBehaviorSanitizer.

NOT_COMPLETED:
- Hamiltonian/geodesic integration, dense output, events, or horizon handling.
- Observer/tetrad, conserved quantities, Kerr separated solver, or proper-time evolution.
- Kerr-Schild coordinates or reliable horizon crossing.
- Radiative transfer, rendering, Solar adapter, WASM/GPU, UI, or movie pipeline.
- DE440's eight external-data assertions on this machine.
- GCC/Linux verification.

CURRENT_BLOCKERS:
- None.

MOST_LIKELY_BUGS:
- Near-extremal spin and points close to the BL margin lack an arbitrary-precision reference sweep.
- The fixed polar-axis and epsilon-scaled Delta floors may need integrator-aware policy in Phase 1.
- Deterministic derivative samples could miss a correlated transcription error outside the tested domain.
- Extreme representable magnitudes can lose precision before domain checks reject an evaluation.
- Only Apple Clang 16 on macOS arm64 was verified locally; GCC/Linux remains dependent on hosted CI or another environment.
- The CLI JSON has no versioned downstream schema yet.
- Missing default DE440 data is a visible skip, so its external validation can regress without a fixture-provisioned job.

FASTEST_WAY_TO_FALSIFY:
- `./tests/relativity/test_kerr_bl`: literals, inverse, signature, surfaces, limits, and invalid domains must all pass.
- `./tests/relativity/test_metric_derivatives`: Dual4 and five-point errors must remain below their gates.
- `./tests/relativity/test_metric_cli`: all three metrics must succeed and malformed/invalid inputs must fail.
- The documented Kerr JSON CLI must report `inverse_error < 5e-13`.
- Kerr CLI at the outer horizon, polar axis, `spin=1`, or a NaN coordinate must exit nonzero.
- Rebuild the relativity tests with ASan/UBSan; any diagnostic invalidates the gate.

NEXT_ALLOWED_ACTION:
- Phase 1 only: implement and independently validate the general Hamiltonian geodesic foundation. Do not enter Phase 2.
