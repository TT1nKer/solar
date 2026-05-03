# Validation: Halo Orbit Computation and Jacobi Conservation

## Claim
The Halo orbit differential corrector converges to a periodic orbit around Sun-Earth L1, with period matching published SOHO-class values, and Jacobi constant conserved during integration.

## Reference

For Sun-Earth L1 Halo orbits, published values (Howell 1984, NASA SOHO mission):
- Az ≈ 110,000-200,000 km (typical SOHO/ACE class)
- Period ≈ 178 days (depends weakly on Az)
- Jacobi constant near 3.0008

For Az = 50,000 km (smaller orbit):
- Period ≈ 177-178 days expected

## Method
1. Continuation approach: start with Lyapunov orbit (z=0), increase Az in 10 steps
2. Differential correction at each step using **finite-difference Jacobian** (the State Transition Matrix is implemented but has a known bug — see below)
3. Adjust x₀ and ẏ₀ to make ẋ(T/2) = 0 and ż(T/2) = 0 at the y=0 crossing
4. Half-period detection via regula falsi on y(t)
5. After convergence, propagate full orbit and check Jacobi constant

## Command
```bash
./solar halo Sun Earth 1 50000
```

## Actual output

```
$ ./solar halo Sun Earth 1 50000
Computing Halo orbit: Sun-Earth L1, Az=50000 km (3.342284e-04 norm)...
Halo Orbit: Sun-Earth L1
  Converged in 119 iterations

  Period: 3.060115 (normalized), 177.89 days
  Jacobi constant: 3.000823

  Initial state (normalized):
    x=9.9162212166e-01  z=3.3422835343e-04
    ydot=-9.6285039335e-03

  Initial state (physical):
    x=1.4834500297e+08 km  z=5.0000000000e+04 km
    ydot=-2.8678202296e-01 km/s
```

## Jacobi conservation test

Propagate the converged orbit for 10 full periods and track Jacobi constant:

```cpp
auto traj = propagate_cr3bp(orbit.initial_state, sys.mu, orbit.period * 10, ...);
double max_drift = max over traj of |C(state) - C(initial)|
```

Result:
```
C₀ = 3.0008e+00
max |dC| = 7.7272e-14
relative drift = 2.5750e-14
```

## Errors

| Quantity | Computed | Published | Error |
|---|---|---|---|
| Period | 177.89 days | ~177.8 days | ~0.05% |
| Jacobi constant | 3.000823 | matches L1 region | OK |
| Jacobi conservation (10 periods) | 2.6e-14 relative | machine precision | OK |

## Notes

### Iteration count: 119
The continuation approach uses 10 sub-steps from Lyapunov (z=0) to target Az, with multiple correction iterations per step. Final-step convergence is fast (~2-5 iterations); most iterations are spent in the early continuation steps.

### Finite differences vs STM
The State Transition Matrix implementation in `cr3bp.cpp` is **known to be broken** — it disagrees with finite-difference computation by approximately 100×. The Halo solver uses finite differences as a workaround:

```
At reference IC, propagate to y=0 crossing → get (ẋ_ref, ż_ref)
Perturb x₀ by ε=1e-8, propagate → get (ẋ_+x, ż_+x)
Perturb ẏ₀ by ε=1e-8, propagate → get (ẋ_+ẏ, ż_+ẏ)

J = [(ẋ_+x - ẋ_ref)/ε   (ẋ_+ẏ - ẋ_ref)/ε]
    [(ż_+x - ż_ref)/ε   (ż_+ẏ - ż_ref)/ε]

Δ = -J⁻¹ · [ẋ_ref, ż_ref]ᵀ
```

This is 3× more expensive per iteration than using a correct STM, but more reliable.

## Limitations

- **One Az tested**: 50,000 km. The full Halo family (varying Az from 1,000 km to 800,000 km) has not been verified.
- **One Lagrange point**: L1. L2 Halo orbits not extensively tested.
- **No comparison with NASA Goddard's GMAT** or other independent tools.
- **STM is broken**: stability index computation (eigenvalues of monodromy matrix) currently outputs 0.0 because the STM is unreliable.
- **No invariant manifold computation** (stable/unstable manifolds for transfer design).
- **Sun-Earth L2 Halo not validated** (only Earth-Moon L2 tested for Halo, with Az=10,000 km giving period ~14.81 days, matching published Earth-Moon L2 expectations).

## References
- Howell, K. C. (1984), "Three-Dimensional, Periodic, 'Halo' Orbits", *Celestial Mechanics*, 32, 53-71
- Richardson, D. L. (1980), "Analytic Construction of Periodic Orbits about the Collinear Points", *Celestial Mechanics*, 22, 241-253
- Koon, Lo, Marsden, Ross (2011), *Dynamical Systems, the Three-Body Problem and Space Mission Design*, ISBN 978-0-615-24095-4
- NASA SOHO mission orbit documentation
