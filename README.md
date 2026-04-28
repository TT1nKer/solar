# Solar System Simulator

A NASA-grade solar system simulation engine written in C++17. Zero external dependencies. From `git clone && make` to sub-meter ephemeris accuracy, Halo orbit computation, and complete mission simulation.

```
22 source files | 24 headers | ~6500 lines C++ | 0 dependencies | 0 warnings
```

## Quick Start

```bash
git clone <repo> && cd solar
make
./solar bodies                              # list all 17 celestial bodies
./solar ephemeris Earth 2026-04-12          # where is Earth today?
./solar transfer Earth Mars 2026-04-12      # Hohmann transfer to Mars
./solar mission mars                        # run a full Mars mission simulation
```

### With DE440 Precision Ephemeris (optional, sub-meter accuracy)

```bash
# Download JPL DE440 data (~30MB for 1950-2050 coverage)
mkdir -p data
curl -o data/header.440 https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/de440/header.440
curl -o data/ascp01950.440 https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/de440/ascp01950.440
cat data/header.440 data/ascp01950.440 > data/de440.asc

# Use it (all commands automatically gain sub-meter accuracy)
./solar --de data/de440.asc ephemeris Earth 2026-04-12
```

## Architecture

```
solar/
├── include/solar/          # 24 public headers
│   ├── body.h              # Vec3, State, OrbitalElements, Body, PhysicalProperties
│   ├── constants.h         # G, AU, MU_*, C_LIGHT, J2000, obliquity
│   ├── kepler.h            # Kepler equation, element<->state conversions
│   ├── ephemeris.h         # Planet/moon data, position queries, JD conversion
│   ├── jpl_de.h            # JPL DE440 binary reader (Chebyshev interpolation)
│   ├── moons.h             # Moon data (Luna, Phobos, Deimos, Galilean, Titan)
│   ├── integrator.h        # RK4, Verlet, DOPRI5 adaptive
│   ├── force_model.h       # ForceModel ABC (modular physics)
│   ├── gravity.h           # Newtonian point-mass gravity
│   ├── j2.h                # J2 oblateness perturbation
│   ├── spherical_harmonics.h # Full Jn zonal harmonics (J2-J6)
│   ├── gr_correction.h     # Schwarzschild 1PN general relativity
│   ├── srp.h               # Solar radiation pressure + shadow model
│   ├── drag.h              # Atmospheric drag
│   ├── atmosphere.h        # Exponential density (Earth, Mars, Venus, Titan)
│   ├── nbody.h             # N-body simulation engine
│   ├── spacecraft.h        # Spacecraft with fuel, thrust, Tsiolkovsky burns
│   ├── transfer.h          # Hohmann, bi-elliptic, Lambert solver
│   ├── trajectory.h        # Porkchop plots, gravity assists, multi-flyby
│   ├── cr3bp.h             # CR3BP, Lagrange points, Halo orbits
│   ├── frame.h             # Mat3, ICRF/Ecliptic/BodyFixed transforms
│   ├── iau_rotation.h      # IAU body rotation (Earth, Moon, Mars)
│   ├── time_scale.h        # TDB/TT/UTC/TAI conversions, leap seconds
│   └── mission.h           # Mission simulation engine
├── src/                    # 22 implementation files
├── cli/main.cpp            # CLI with 18 commands + 8 physics flags
├── tests/                  # Validation tests
├── data/                   # DE440 ephemeris data (user-downloaded)
└── Makefile                # `make` builds libsolar.a + ./solar
```

## Commands

| Command | Description |
|---------|-------------|
| `bodies` | List all 17 bodies with radius, gravity, atmosphere |
| `info <body>` | Detailed properties + orbital elements + children |
| `moons <planet>` | Moon table with periods and properties |
| `ephemeris <body> <date>` | Position/velocity at date (dual frame for moons) |
| `orbits` | Orbital elements table for all planets |
| `simulate <days> <dt> [interval]` | N-body simulation with configurable physics |
| `energy <days> <dt>` | Energy conservation diagnostic |
| `transfer <from> <to> <date>` | Hohmann transfer calculation |
| `lambert <from> <to> <dep> <arr>` | Lambert problem (realistic trajectory) |
| `porkchop <from> <to> <d1> <d2> <a1> <a2>` | Porkchop plot (TSV, pipeable to gnuplot) |
| `launch-window <from> <to> ...` | Find optimal launch windows per synodic period |
| `flyby <body> <date> <v> <rp>` | Gravity assist calculator |
| `multi-flyby <from> <to> ... --via ...` | Chained multi-flyby trajectory |
| `frame <body> <date> <frame>` | Position in ICRF / Ecliptic / BodyFixed |
| `time <date>` | Time scale conversions (UTC/TAI/TT/TDB) |
| `lagrange <primary> <secondary>` | Lagrange points L1-L5 |
| `halo <primary> <secondary> <L> <Az>` | Halo orbit computation |
| `mission <template>` | Full mission simulation (mars / gateway / grand-tour) |

### Physics Flags

```bash
./solar simulate 365 3600 --j2 --gr --adaptive   # all the physics
./solar energy 365 3600 --harmonics               # J2-J6 spherical harmonics
./solar simulate 30 60 --srp --drag               # with radiation + drag
```

## Development Phases

### Phase 1: Core Mathematics
Kepler equation solver (Newton-Raphson), orbital element conversions, vis-viva equation, N-body gravitational simulation (O(N^2) direct summation), Velocity Verlet symplectic integrator, RK4, Lambert solver (universal variable with Stumpff functions), Hohmann and bi-elliptic transfers.

**Validation**: Kepler solver 16/16 tests pass. Hohmann Earth-Mars: 5.59 km/s, 259 days (textbook match). Energy drift 7.7e-6 over 1 year with Verlet.

### Phase 2: Solar System Model
17 celestial bodies: Sun + 8 planets + Luna, Phobos, Deimos, Io, Europa, Ganymede, Callisto, Titan. Hierarchical parent-child tree (parent_index). Physical properties for each body (radius, surface gravity, escape velocity, atmosphere, pressure). Spacecraft struct with fuel tracking and Tsiolkovsky burns.

**Validation**: All bodies load correctly. Moon orbital period 27.3 days. Titan atmosphere 1.45 atm. Backward compatible with Phase 1.

### Phase 3: Modular Physics Engine
ForceModel abstract base class with `compute()` and `potential_energy()`. NBodySim holds composable `vector<unique_ptr<ForceModel>>`. Concrete models: NewtonianGravity, J2Perturbation, GRCorrection (Schwarzschild 1PN). Auto-defaults to gravity if no forces added.

**Validation**: Energy drift identical before/after refactor. J2 adds measurable perturbation. GR adds ~1e-11 correction. Forces compose by summation.

### Phase 4: JPL DE440 Precision Ephemeris
Reader for JPL Development Ephemeris ASCII format. Chebyshev polynomial interpolation via Clenshaw's algorithm. ICRF-to-ecliptic frame rotation. Barycentric-to-heliocentric conversion. EMB+Moon decomposition for Earth position. Automatic fallback to Keplerian without DE file.

**Validation**: Earth position at J2000 matches JPL Horizons to **< 1 meter**. Velocity to **< 0.01 mm/s**. 8/8 DE tests pass.

### Phase 5: Coordinate Frames + Time System
Mat3 rotation matrix type. Frame enum (EclipticJ2000, ICRF, BodyFixed). `transform_state()` routes through ICRF hub. IAU rotation models for Earth, Moon, Mars (prime meridian angle, pole direction). Time scales: TDB, TT, UTC, TAI with leap second table (1972-2017). Epoch struct with automatic conversion.

**Validation**: TT-UTC = 69.184s (37 leap seconds + 32.184). Sub-solar point longitude matches UTC midnight. ICRF/Ecliptic roundtrip preserves magnitude.

### Phase 6: Higher-Order Force Models
Atmospheric density model (exponential: Earth, Mars, Venus, Titan). Spherical harmonics gravity (J2-J6 for Earth/Mars/Jupiter/Saturn). Solar radiation pressure with cylindrical shadow model. Atmospheric drag with co-rotation correction.

**Validation**: Spherical harmonics J2 regression matches J2Perturbation to 2e-16 (machine precision). SRP at 1 AU = 4.7e-9 km/s^2 (correct). J3-J6 contributes 3.8e-5 of J2 (physical).

### Phase 7: Adaptive Integrator
Dormand-Prince RK5(4) embedded pair with automatic step size control. Error estimation via 4th-order solution. Step controller with safety factor and growth limiter. Generic DOPRI5 for arbitrary ODE systems.

**Validation**: Angular momentum conservation 200x better than fixed Verlet (9.7e-11 vs 1.8e-8). Rejected steps = 2 over 1 year (well-tuned controller). Energy drift comparable to symplectic Verlet (expected for non-symplectic methods).

### Phase 8: Trajectory Optimization
Porkchop plot generator (departure x arrival date sweep). Optimal launch window finder per synodic period. Gravity assist calculator (patched-conic hyperbolic flyby). Multi-flyby trajectory evaluator with v-infinity matching.

**Validation**: Earth-Mars optimal window 2026-10-26, dv=5.64 km/s. 4 windows found at 780-day intervals (Earth-Mars synodic period). Jupiter flyby turn angle physically reasonable.

### Phase 9: Three-Body Dynamics (CR3BP)
Circular Restricted Three-Body Problem equations of motion in rotating frame. Lagrange point computation (Newton-Raphson for L1-L3, analytical for L4-L5). Halo orbit computation via continuation + finite-difference differential correction. Jacobi constant conservation tracking. State Transition Matrix propagation.

**Validation**: Sun-Earth L1 at 1.492e6 km from Earth (published: ~1.5e6 km). Sun-Earth L1 Halo period = **177.89 days** (published: ~177.8 days). Jacobi constant conserved to **2.6e-14** over 10 periods.

### Phase 10: Mission Simulation
Mission engine: event sequencing (launch, coast, burn, flyby, orbit insertion). Delta-v budget tracking via Tsiolkovsky rocket equation. Fuel accounting through mission phases. Three pre-built templates: Mars Direct, Lunar Gateway, Voyager Grand Tour.

**Validation**: Mars Direct 2026: departure 3.62 km/s + MOI 2.05 km/s = 5.67 km/s total, 312 days, SUCCESS. Lunar Gateway: TLI 3.13 + LOI 0.82 = 3.95 km/s, 4 days, SUCCESS. Grand Tour: Jupiter flyby 116 deg turn + Saturn flyby 62 deg, free dv gains of 13.5 + 14.5 km/s.

## Accuracy

| Metric | Without DE440 | With DE440 |
|--------|:---:|:---:|
| Earth position | ~7000 km (0.005%) | **< 1 meter** |
| Earth velocity | ~12 m/s | **< 0.01 mm/s** |
| N-body energy (1yr) | 7.7e-6 relative | 7.7e-6 relative |
| Jacobi constant (CR3BP) | — | 2.6e-14 relative |
| Halo orbit period | — | 0.05% of published |

## Physical Models

| Model | Implementation | Reference |
|-------|---------------|-----------|
| Gravity | O(N^2) direct summation, mu=GM precision | Newton |
| Kepler equation | Newton-Raphson, 1e-12 tolerance | Battin (1999) |
| Spherical harmonics | Recursive Legendre, J2-J6 zonal | Montenbruck & Gill |
| General relativity | Schwarzschild 1PN correction | Soffel et al. (2003) |
| Solar radiation pressure | 1361 W/m^2, cylindrical shadow | Montenbruck & Gill |
| Atmospheric drag | Exponential with co-rotation | Vallado (2007) |
| Lambert solver | Universal variable, Stumpff functions | Battin (1999) |
| CR3BP | Normalized rotating frame, Jacobi integral | Szebehely (1967) |
| Halo orbits | Continuation + differential correction | Howell (1984) |
| Integrators | Verlet (symplectic), RK4, DOPRI5 (adaptive) | Hairer et al. (1993) |
| Ephemeris | JPL DE440 Chebyshev interpolation | Folkner et al. (2014) |
| Time scales | TDB-TT Fairhead & Bretagnon, leap seconds | IAU/IERS |
| Body rotation | IAU 2015 pole/prime meridian models | Archinal et al. (2018) |

## Build

```bash
make            # builds libsolar.a + ./solar
make test       # runs validation suite
make clean      # removes build artifacts
```

Requirements: C++17 compiler (GCC 9+, Clang 10+). No external libraries.

## License

Research and educational use.
