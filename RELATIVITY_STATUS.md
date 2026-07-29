# Solar Relativity Status

CURRENT_PHASE: -1
PHASE_STATE: IN_PROGRESS
LAST_VERIFIED_COMMIT: e33f1d7
LAST_VERIFIED_PLATFORM: macOS 14.8.7 arm64 / Apple Clang 16.0.0
LAST_VERIFIED_COMMANDS:
- `make clean && make && make test`
- `./tests/test_de /tmp/solar-definitely-missing-de440.asc` (expected exit 1)
- README smoke commands listed in `docs/validation/relativity_00_baseline.md`

ACTUALLY_COMPLETED:
- Audited the existing Solar API, build, CLI, tests, and model boundary.
- Reproduced and recorded the untouched `main` baseline at commit `1268d89`.
- Fixed the macOS direct-include failure in `test_validation`.
- Made absent default DE440 data an explicit skip while preserving explicit-path failure.
- Made Make discover nested sources/tests and propagate child test failures.

NOT_COMPLETED:
- Phase 0A unit, math, type, and Dual4 implementation.
- Any metric, Hamiltonian, geodesic, Kerr, observer, rendering, adapter, or UI code.

CURRENT_BLOCKERS:
- None for Phase 0A.

MOST_LIKELY_BUGS:
- Recursive Make patterns have not yet been exercised by real nested sources.
- The DE440 skip path can hide loss of optional external validation unless CI provisions the fixture.
- Existing generic DOPRI5 and NBodySim time handling are not suitable for Hamiltonian geodesics without redesign.

FASTEST_WAY_TO_FALSIFY:
- Add `tests/relativity/test_probe.cpp`, run `make test`, and verify it builds and runs.
- Run `./tests/test_de /tmp/solar-definitely-missing-de440.asc`; any zero exit is wrong.
- Replace one existing test command with `false`; `make test` must exit nonzero.

NEXT_ALLOWED_ACTION:
- Finish Phase −1 evidence documents, then implement Phase 0A by failing tests first.
