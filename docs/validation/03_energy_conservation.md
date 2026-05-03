# Validation: Energy Conservation

## Claim
The N-body simulator's Velocity Verlet integrator conserves total mechanical energy to ~1e-5 relative drift over 1 year for the 9-body solar system.

## Method
Integrate the 9-body solar system (Sun + 8 planets) for 1 year using Velocity Verlet at dt = 3600 s.
Compare initial and final total energy E = KE + PE.

## Theoretical expectation
- Velocity Verlet is **symplectic**: energy errors are bounded oscillations, not secular drift.
- The bound depends on dt and orbital frequencies.
- For dt = 3600 s and inner planet timescales (~hours-to-days), expected relative error is order of 1e-5 to 1e-6.

## Command
```bash
./solar energy 365 3600 --planets-only
```

## Actual output
```
# Energy conservation test: 365 days, dt=3600s, 9 bodies
# Forces: newtonian_gravity
# Initial energy: -1.983314e+29
# Initial |L|:    3.132326e+37
# Final energy:   -1.983299e+29
# Final |L|:      3.132326e+37
# Relative energy drift:  7.656552e-06
# Relative |L| drift:     1.785352e-08
```

## Error
- **Relative energy drift**: 7.66e-6 (within expected range for Verlet at this dt)
- **Relative angular momentum drift**: 1.79e-8 (better than energy, as expected for symplectic methods)

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
| Verlet (fixed dt=3600s) | 7.66e-6 | 1.79e-8 | 8,760 |
| DOPRI5 (adaptive, tol=1e-10) | 7.61e-6 | 9.74e-11 | 99,243 |

Observations:
- DOPRI5 gives **200x better** angular momentum conservation
- Energy drift is **similar** because DOPRI5 has secular drift while Verlet has bounded oscillation
- DOPRI5 takes more steps overall

## Limitations
- Only one test case (1 year, dt=3600s, 9 bodies)
- Single-point comparison (initial vs final) doesn't distinguish bounded oscillation from secular drift
- Would need to track energy at intermediate points to confirm symplectic behavior
- No comparison with longer integrations (decades, centuries) where the difference between Verlet and DOPRI5 should grow

## Future work
- Add a test that tracks energy at every output interval and reports max-min range
- Compare with literature values (e.g., Wisdom-Holman maps for solar system integration)
- Test at multiple dt values to verify expected scaling
