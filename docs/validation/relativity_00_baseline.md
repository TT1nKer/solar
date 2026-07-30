# Validation: Solar baseline before relativity Phase 0A

## Claim

The unmodified `main` library and CLI build on the audit platform, but the
test baseline is red because of one missing direct include. After the bounded
Phase −1 repairs, all fixture-independent legacy tests and README smoke
commands run, missing optional DE440 data is visible as a skip, and an explicit
bad DE440 path remains a failure.

## Model boundary

This validates existing Solar build, Newtonian/1PN educational models, CLI
examples, and test wiring. It validates no metric, geodesic, Kerr, horizon, or
radiative-transfer physics.

## Platform

```text
macOS 14.8.7 arm64 (Build 23J520)
Apple clang version 16.0.0 (clang-1600.0.26.6)
Target: arm64-apple-darwin23.6.0
```

## Untouched reference

```text
git status --short
(no output after quarantining four generated PR #2 artifacts)

git rev-parse HEAD
1268d89f5f84779d2fc66a1487280af586d1a245
```

The quarantined files were generated Mach-O executables
`tests/test_dynamics`, `tests/test_relativity` and their `.dSYM` directories.
They contained no source and were moved to
`/tmp/solar-spike-artifacts.rtcBkM`.

## Original commands and output

`make clean` removed the source objects, `libsolar.a`, `solar`, and the six
top-level test binaries.

`make` exited 0. It compiled 24 source files, archived `libsolar.a`, and linked
`solar`. The complete warning set was:

```text
src/cr3bp.cpp:388:9: warning: variable 'sub_iter' set but not used
src/integrator.cpp:100:25: warning: unused variable 'dp_a71'
src/integrator.cpp:100:51: warning: unused variable 'dp_a73'
src/integrator.cpp:100:76: warning: unused variable 'dp_a74'
src/integrator.cpp:100:102: warning: unused variable 'dp_a75'
src/integrator.cpp:100:127: warning: unused variable 'dp_a76'
src/integrator.cpp:114:55: warning: unused variable 'dp_c7'
```

`make test` exited 2 while compiling `tests/test_validation.cpp`:

```text
tests/test_validation.cpp:54:28: error: implicit instantiation of undefined
template 'std::basic_ostringstream<char>'
/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1/
__fwd/sstream.h:27:28: note: template is declared here
```

The same error occurred at lines 60, 94, 99, 115, 143, 158, 174, 184, and
215. Root cause: the file used `std::ostringstream` without including
`<sstream>` and compiled only where transitive includes happened to expose it.

After adding the direct include, `make test` exposed a second baseline defect:
`test_de` returned 1 for absent `data/de440.asc`, all later executables passed,
and the Make recipe still exited 0 because its shell loop returned the last
test's status.

## Phase −1 repairs

- `tests/test_validation.cpp`: include `<sstream>` directly.
- `tests/test_de.cpp`: missing default optional data prints `SKIP` and exits 0;
  an explicit invalid path exits 1.
- `Makefile`: recursively discover nested sources/tests and aggregate any child
  failure.

These changes do not alter library physics or public APIs.

## Repaired test output

Command:

```bash
make clean && make && make test
```

Result: exit 0.

```text
test_de:         SKIP (default data/de440.asc absent)
test_horizons:   print-only comparison, exit 0
test_kepler:     16 passed, 0 failed
test_montecarlo: 13 passed, 0 failed
test_network:    13 passed, 0 failed
test_validation: 14 passed, 0 failed
```

Fixture-independent assertion total: 56 passed, 0 failed. The same seven
legacy compiler warnings remain. DE440's eight data-backed assertions were not
run and are not counted as passes.

Explicit failure check:

```bash
./tests/test_de /tmp/solar-definitely-missing-de440.asc
```

```text
DE: cannot open /tmp/solar-definitely-missing-de440.asc
Cannot load DE file: /tmp/solar-definitely-missing-de440.asc
```

Actual exit: 1, expected exit: 1.

## README smoke commands

All commands exited 0:

```bash
./solar bodies
./solar transfer Earth Mars 2026-04-12
./solar lagrange Sun Earth
./solar halo Sun Earth 1 50000
./solar energy 365 3600 --planets-only
./solar lambert Earth Mars 2026-10-26 2027-09-03
```

Observed values:

```text
bodies: 17 rows headed by Sun, Mercury, Venus, Earth, Mars
transfer total delta-v: 5.594 km/s
transfer time of flight: 258.871 days
Sun-Earth L1 x: 0.990027 normalized
halo convergence: 140 iterations
halo period: 177.89 days
energy relative drift: 7.656552e-06
angular-momentum relative drift: 1.785352e-08
Lambert total delta-v: 5.642 km/s
Lambert time of flight: 312.000 days
```

## Result

Phase −1 baseline behavior is repaired and reproducible on this platform.
Nested Phase 0A sources subsequently proved recursive build, header dependency,
test execution, and clean behavior.

## Phase 0A final verification

Verified commit:

```text
210450b1db3a9903a608db74a67c75c6b847b28b
```

Commands:

```bash
make clean
test ! -e src/relativity/units.o
test ! -e src/relativity/units.d
test ! -e tests/relativity/test_units
test ! -e tests/relativity/test_units.d
make
make test
./tests/relativity/test_units
./tests/relativity/test_math
./tests/relativity/test_types
./tests/relativity/test_dual4
git diff --check
```

Actual result:

```text
Phase 0A test_units: 14 passed, 0 failed
Phase 0A test_math:  15 passed, 0 failed
Phase 0A test_types: 8 passed, 0 failed
Phase 0A test_dual4: 25 passed, 0 failed
Legacy test_kepler: 16 passed, 0 failed
Legacy test_montecarlo: 13 passed, 0 failed
Legacy test_network: 13 passed, 0 failed
Legacy test_validation: 14 passed, 0 failed
Legacy test_de: SKIP because default external fixture is absent
Legacy test_horizons: print-only, exit 0
```

Unique assertion result: 118 passed, 0 failed, eight DE440 assertions skipped.
All required commands exited 0. The seven pre-existing compiler warnings remain.

Header dependency check:

```bash
touch include/solar/relativity/dual4.h
make tests/relativity/test_dual4
```

Actual output included a new compile command for
`tests/relativity/test_dual4.cpp`; the stale-binary defect is closed.

## Limitations

- DE440 external data was absent; its eight assertions were skipped.
- `test_horizons` prints differences but has no pass/fail thresholds.
- The seven existing compiler warnings remain.
- Smoke commands reproduce documented examples; they do not independently
  validate the physical models.
- Phase 0A validates foundational arithmetic and derivatives, not any spacetime
  metric or geodesic.

## Fastest falsification

1. Remove `<sstream>` and run `make test` with Apple Clang; compilation must
   fail.
2. Run `./tests/test_de /tmp/solar-definitely-missing-de440.asc`; a zero exit
   means explicit input failure is being hidden.
3. Add a nested failing `tests/relativity/test_probe.cpp`; `make test` must
   build it and return nonzero.
4. Change a nested source and run `make`; its object must be rebuilt into
   `libsolar.a`.
