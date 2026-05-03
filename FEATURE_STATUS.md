# Feature Status

Status definitions:
- **Implemented**: Has unit tests + sample run + matches reference values
- **Experimental**: Code runs, basic regression only, no external validation
- **Research prototype**: Algorithm implemented but needs independent validation
- **Toy model**: Demonstrates concept but not physically complete

| Module | Status | Validation |
|---|---|---|
| Kepler equation solver | Implemented | 16/16 unit tests, Newton-Raphson converges to 1e-12 |
| Orbital elements ↔ state | Implemented | Round-trip test passes |
| Hohmann transfer | Implemented | Earth-Mars: 5.59 km/s, 259d (matches textbook) |
| N-body Verlet integrator | Implemented | Energy drift 7.7e-6 over 1 year (9 bodies) |
| RK4 integrator | Implemented | Standard regression |
| DOPRI5 adaptive integrator | Experimental | Angular momentum 200x better than Verlet, but no analytic test |
| DE440 ASCII parser | Experimental | Earth at J2000 matches Horizons to <1m, only 1 test case |
| Frame transforms (ICRF/Ecliptic) | Implemented | Roundtrip preserves magnitude |
| Body-fixed frame (IAU rotation) | Experimental | Sub-solar point sanity check only |
| Time scales (TDB/TT/UTC/TAI) | Implemented | Leap seconds 1972-2017 hardcoded |
| J2 perturbation | Experimental | Order-of-magnitude check |
| Spherical harmonics (J2-J6) | Experimental | Matches J2 model to machine precision; no Jn-only validation |
| GR Schwarzschild correction | Experimental | Order-of-magnitude check, no Mercury precession test |
| Solar radiation pressure | Experimental | 1 AU magnitude correct; shadow logic untested |
| Atmospheric drag | Experimental | Hand calculation match; ISS density underestimated (single-exp) |
| Lambert solver | Experimental | Earth-Mars regression matches |
| Porkchop plot | Experimental | Min dv at 2026 window matches expected ~5.6 km/s |
| Launch window finder | Experimental | Synodic period 780d (correct) |
| Gravity assist | Experimental | Hyperbolic geometry only; no published case validation |
| Multi-flyby | Toy model | v-infinity mismatch shown but not optimized |
| CR3BP equations of motion | Experimental | Jacobi constant conserved to 2.6e-14 |
| Lagrange points L1-L5 | Implemented | L1 Sun-Earth = 1.49e6 km matches published |
| Halo orbit (Lyapunov + continuation) | Research prototype | Period matches published to ~0.1% but no family verified |
| State Transition Matrix | Known broken | Disagrees with finite differences by ~100x |
| Mission simulation engine | Toy model | Demonstrates sequencing + Tsiolkovsky only |
| Mars Direct template | Toy model | dv budget plausible, no operational realism |
| Lunar Gateway template | Toy model | Uses hardcoded TLI Δv (3.13 km/s) and LOI Δv (0.82 km/s); does not actually compute Earth-Moon trajectory |
| Voyager Grand Tour template | Toy model | Demonstrates flyby concept only |
| Monte Carlo runner | Experimental | Produces statistics but limited perturbation model |
| Sensitivity analysis | Experimental | One-at-a-time only, not global |
| Solar transfer network | Experimental | Lambert edges + Dijkstra; single epoch only |

## Headline numbers (with context)

**Earth position at J2000 vs JPL Horizons**: position error <1 m, velocity error <0.01 mm/s.
*Caveat*: This is one date, one body. Not a comprehensive validation.

**Sun-Earth L1 Halo period**: 177.89 days (published ~177.8 days).
*Caveat*: One specific Az amplitude. The full Halo family has not been verified.

**Jacobi constant conservation**: 2.6e-14 over 10 periods.
*Caveat*: Specific to the converged orbit. Other initial conditions may show different drift.

**Energy conservation (1 year, dt=3600s, Verlet)**: 7.7e-6 relative drift.
*Caveat*: Solar system 9-body case only. Different timesteps and force models give different results.

## What's NOT validated

- Multi-rev Lambert transfers
- Long-term ephemeris drift (decades)
- Non-Earth ephemeris precision (only Earth tested at J2000)
- IAU rotation models for body-fixed frames at non-J2000 epochs
- Lunar Halo orbits at L1
- Resonant flyby sequences
- Tesseral spherical harmonics (not implemented)
- Drag at realistic LEO altitudes (single-exp model insufficient)
