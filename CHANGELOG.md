# Changelog

All notable changes to this project will be documented in this file.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Changed
- **Major rescope**: project repositioned from "NASA-grade simulation engine" to "experimental educational orbital mechanics sandbox". README rewritten to be honest about what is and isn't validated.
- README headline numbers (Hohmann Δv, L1 distance, Halo period, energy drift, Lambert Δv) now show reference value + computed value + match percentage + link to validation report. Every claim is reproducible from a single command.
- Mission template "Lunar Gateway" renamed to "Lunar Gateway (Δv-accounting demo)". Added in-code LIMITATION comment explaining that the Δv values are hardcoded textbook constants, not computed from ephemeris (Earth-Moon transfer is not a heliocentric Lambert problem).

### Added
- `tests/test_validation.cpp` — automated test that cross-checks **every numerical claim in the README and validation docs**. Runs as part of `make test`. Currently 14 assertions covering reports 02, 03, 05, 06, 07, 08. If any claim becomes inaccurate, this test fails and CI breaks.
- `CHANGELOG.md` — this file.
- `LIMITATIONS.md` — 9 categories of explicit non-use cases, including: not for operational mission design, not safety-critical, simplified atmospheric model, single-revolution Lambert only, mission templates are toy models, etc.
- `FEATURE_STATUS.md` — per-module status across 30 features, classified into 4 levels: Implemented / Experimental / Research prototype / Toy model. Each entry has a one-line validation note.
- `STRATEGY.md` — strategic roadmap from orbital engine toward broader Universe-OS direction (vehicle layer, uncertainty layer, infrastructure network, resource flow, agents).
- `docs/validation/` — 8 reproducible validation reports (Kepler solver, Earth-Mars Hohmann, energy conservation, DE440 vs JPL Horizons, Lagrange points, Halo orbit + Jacobi conservation, frame transforms, Lambert solver). Each report has fixed structure: Claim / Reference / Method / Command / Actual output / Errors / Notes / Limitations / Future validation needed.
- `.github/workflows/ci.yml` — CI workflow that runs `make` (must produce zero warnings) → `make test` (must pass 13/13) → smoke test of CLI commands on every push.
- `CHANGELOG.md` (this file).

### Removed
- **State Transition Matrix** (`CR3BPStateSTM`, `propagate_cr3bp_stm`, internal STM helpers): disagreed with finite-difference Jacobian by ~100×. Removed from API rather than ship broken code. Halo orbit corrector uses finite differences (slower but verified). Header file documents the removal.
- `propagate_to_y_crossing` no longer accepts an STM output parameter.
- Stability index computation (was always 0.0 placeholder due to STM bug).

### Fixed
- README numerical claims now match actual command output. Previously "5.656 km/s" Lambert reference was off by 0.014 km/s; corrected to "5.642 km/s" matching current code.
- Mission engine `apply_burn`: now returns `BurnResult` struct with explicit `complete` flag, so partial burns due to fuel exhaustion are properly detected and reported.

### Honest about
- DE440 validation: only **one date, one body** (Earth at J2000) tested
- Halo orbit: only **one Az amplitude** (50,000 km) verified
- Frame transforms: precession/nutation **not** modeled (constant J2000 obliquity only)
- Lambert: **single-revolution only**, not validated against Izzo or Gooding
- Drag at LEO altitude: simple exponential model **underestimates** real thermosphere

---

## [Phase 1-10] - 2026-04-12

Original development phases (in single squashed commit). See `STRATEGY.md` for the phase-by-phase walkthrough. The codebase reached:
- 22 source files, 24 headers, ~6500 lines of C++17
- Zero external dependencies
- Implemented: Kepler solver, N-body, Verlet/RK4/DOPRI5 integrators, JPL DE440 parser, Lambert solver, porkchop, gravity assists, CR3BP, Lagrange points, Halo orbits, mission engine, Monte Carlo, solar transfer network

Subsequent rewrite (above) addressed credibility and validation rather than adding features.
