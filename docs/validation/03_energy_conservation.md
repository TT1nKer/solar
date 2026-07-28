# Validation: Energy Conservation

## Claim
The N-body simulator's Velocity Verlet integrator has a 5.0e-12 relative endpoint energy error over 1 year for the 9-body solar system.

## Method
Integrate the 9-body solar system (Sun + 8 planets) for 1 year using Velocity Verlet at dt = 3600 s.
Compare initial and final total energy E = KE + PE.

## Theoretical expectation
- Velocity Verlet is **symplectic**: energy errors are bounded oscillations, not secular drift.
- The bound depends on dt and orbital frequencies.
- The reported value is an endpoint comparison, not a bound over the full trajectory.

## Command
```bash
./solar energy 365 3600 --planets-only
```

## Actual output
```
# Energy conservation test: 365 days, dt=3600s, 9 bodies
# Forces: newtonian_gravity
# Initial energy: -1.980931e+29
# Initial |L|:    3.132030e+37
# Final energy:   -1.980931e+29
# Final |L|:      3.132030e+37
# Relative energy drift:  5.002712e-12
# Relative |L| drift:     1.206212e-15
```

## Error
- **Relative endpoint energy error**: 5.00e-12
- **Relative endpoint angular-momentum error**: 1.21e-15
- Earlier reports mixed `mu`-based accelerations with `G*m` diagnostics. The current diagnostic derives effective inertial masses from `mu/G`, matching the integrated equations.

## Notes
- The "drift" here is the difference between initial and final energy. For a true symplectic integrator with bounded oscillation, this should be small but not zero (it depends on phase of the oscillation at start and end).
- A better test would track energy over many evaluations and confirm the error stays **bounded** rather than grows linearly.
- Higher-order symplectic methods (e.g., Yoshida 4th-order) would give smaller errors at the same dt.

## Comparison: Verlet vs DOPRI5 (non-symplectic)

```bash
./solar energy 365 3600 --planets-only --adaptive
```

| Method | Energy drift | Angular momentum drift | Steps |
|---|---|---|---|
| Verlet (fixed dt=3600s) | 5.00e-12 | 1.21e-15 | 8,760 |
| DOPRI5 (max chunk=3600s, tol=1e-10) | 2.84e-15 | 1.51e-15 | 8,760 |

Observations:
- Both methods are near floating-point limits for this endpoint comparison.
- This result does not establish long-term boundedness or distinguish secular from oscillatory error.

## Limitations
- Only one test case (1 year, dt=3600s, 9 bodies)
- Single-point comparison (initial vs final) doesn't distinguish bounded oscillation from secular drift
- Would need to track energy at intermediate points to confirm symplectic behavior
- No comparison with longer integrations (decades, centuries) where the difference between Verlet and DOPRI5 should grow

## Future work
- Add a test that tracks energy at every output interval and reports max-min range
- Compare with literature values (e.g., Wisdom-Holman maps for solar system integration)
- Test at multiple dt values to verify expected scaling
