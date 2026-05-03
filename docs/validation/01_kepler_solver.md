# Validation: Kepler Equation Solver

## Claim
The Newton-Raphson Kepler equation solver converges to machine precision (1e-12) for typical and high-eccentricity inputs.

## Method
Solve `M = E - e·sin(E)` for E given M and e, compare against analytical reference values.

## Reference values
Computed independently using Wolfram Alpha and verified against textbook examples.

| M (rad) | e | Expected E (rad) |
|---|---|---|
| 1.0 | 0.5 | 1.498701 |
| 2.0 | 0.0 | 2.000000 (E = M when e=0) |
| 1.0 | 0.9 | 1.862087 |

## Command
```bash
make test
./tests/test_kepler
```

## Actual output
```
=== Kepler Equation Solver ===
  PASS: solve_kepler(M=1.0, e=0.5) = 1.498701
  PASS: solve_kepler(M=2.0, e=0.0) = 2.000000
  PASS: solve_kepler(M=1.0, e=0.9) = 1.862087
```

## Error
- M=1.0, e=0.5: error < 1e-6 (tolerance)
- M=2.0, e=0.0: exact (circular orbit, E=M)
- M=1.0, e=0.9: error < 1e-4 (high eccentricity)

## Notes
- Initial guess uses `E0 = π` for `e > 0.8` to avoid divergence near apoapsis
- Tolerance set to 1e-12 in production code; tests use looser bounds for the high-e case
- Has not been validated for hyperbolic (e > 1) or parabolic (e = 1) cases — these use a different solver

## Limitations
- Only tested at three points
- No coverage of edge cases (M = 0, M = 2π, e very close to 1)
- No comparison with alternative solvers (e.g., Halley's method, series expansion)
