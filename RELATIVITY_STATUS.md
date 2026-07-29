# Solar Relativity Status

CURRENT_PHASE: 0A
PHASE_STATE: PASSED
LAST_VERIFIED_COMMIT: 210450b1db3a9903a608db74a67c75c6b847b28b
LAST_VERIFIED_PLATFORM: macOS 14.8.7 arm64 / Apple Clang 16.0.0
LAST_VERIFIED_COMMANDS:
- `make clean`
- nested object, dependency, and test absence checks
- `make`
- `make test`
- `./tests/relativity/test_units`
- `./tests/relativity/test_math`
- `./tests/relativity/test_types`
- `./tests/relativity/test_dual4`
- `./tests/test_de /tmp/solar-definitely-missing-de440.asc` (expected exit 1)
- README smoke commands listed in `docs/validation/relativity_00_baseline.md`
- `git diff --check`

ACTUALLY_COMPLETED:
- Audited the existing Solar API, build, CLI, tests, and model boundary.
- Reproduced and recorded the untouched `main` baseline at commit `1268d89`.
- Fixed the macOS direct-include failure in `test_validation`.
- Made absent default DE440 data an explicit skip while preserving explicit-path failure.
- Made Make discover nested sources/tests, track header dependencies, clean nested artifacts, and propagate child test failures.
- Fixed signature, coordinate, spin-axis, geometric-unit, and null proper-time conventions.
- Added checked SI/geometrized unit conversions.
- Added fixed-size vectors, explicit Minkowski contraction, 4x4 multiplication/inversion, finite checks, and tensor-variance wrappers.
- Added PhaseSpaceState and GeodesicSample value contracts.
- Added four-coordinate forward automatic differentiation with checked domain failures.
- Passed 62/62 Phase 0A assertions and 56/56 fixture-independent legacy assertions.

NOT_COMPLETED:
- Minkowski, Schwarzschild, or Kerr metric implementations.
- Hamiltonian/geodesic integration, dense output, events, or horizon handling.
- Observer/tetrad, conserved quantities, radiative transfer, rendering, Solar adapter, WASM/GPU, or UI.
- DE440's eight external-data assertions on this machine.

CURRENT_BLOCKERS:
- None. Phase 0B still requires an explicit continue instruction.

MOST_LIKELY_BUGS:
- The fixed `1.98847e30 kg` solar-mass convention may differ from a downstream caller that expects the IAU nominal solar mass parameter.
- The 4x4 inverse intentionally rejects pivots below an epsilon-scaled threshold and may reject a mathematically invertible, badly scaled matrix.
- Dual4 has analytic unit tests but no finite-difference sweep over randomized expressions or extreme magnitudes.
- `PhaseSpaceState::x` follows v3 as `Contravariant4`; treating a coordinate point as a vector could enable a semantically invalid future operation.
- Only Apple Clang 16 on macOS arm64 was verified locally; GCC/Linux remains dependent on hosted CI or another environment.
- Missing default DE440 data is a visible skip, so its external validation can regress without a fixture-provisioned job.

FASTEST_WAY_TO_FALSIFY:
- `./tests/relativity/test_units`: one solar mass must yield `1476.6696910334392 m` and round trips within recorded tolerances.
- `./tests/relativity/test_math`: the hand-set matrix times its inverse must differ from identity by less than `1e-14`.
- `./tests/relativity/test_dual4`: `x*x*y+sin(y)` at `(2,3)` must give derivatives `(12, 4+cos(3), 0, 0)`.
- Touch `include/solar/relativity/dual4.h`; `make tests/relativity/test_dual4` must recompile the test.
- Run `./tests/test_de /tmp/solar-definitely-missing-de440.asc`; any zero exit is wrong.

NEXT_ALLOWED_ACTION:
- Phase 0B only: metric interface plus independently validated Minkowski, Schwarzschild, and Kerr Boyer-Lindquist metrics. Stop until explicitly instructed.
