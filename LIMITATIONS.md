# Limitations

This project is **educational and experimental**. It is not intended to replace professional astrodynamics tools such as GMAT, STK, SPICE, or Orekit.

## Known limitations

### Numerical accuracy
- Numerical accuracy claims are only valid for the **documented test cases**.
- Comparisons against JPL Horizons have only been performed for a small number of bodies and dates.
- No comprehensive cross-validation against multiple reference epochs or bodies.

### DE440 ephemeris parser
- Supports the **old JPL ASCII format only** (`header.440` + `ascpXXXX.440`).
- Does not support the SPICE SPK/BSP binary format.
- Tested only for the 1950-2050 data file (`ascp01950.440`).
- No checksum verification of input files.

### Force models
- Atmospheric drag uses **simplified single-exponential** atmosphere models. Real thermosphere density at LEO is significantly different.
- Solar radiation pressure uses **cylindrical shadow model** (no penumbra).
- Spherical harmonics include **zonal terms only** (Jn). No tesseral (Cnm/Snm) terms.
- General relativity correction is **1PN Schwarzschild only**, no Lense-Thirring or higher-order terms.

### Integrators
- Verlet is symplectic (good for long-term) but only 2nd order.
- DOPRI5 is non-symplectic and has secular energy drift over very long integrations.
- No symplectic methods of higher order than Verlet (e.g., no Yoshida).
- No Bulirsch-Stoer or Adams-Bashforth-Moulton.

### Trajectory tools
- Lambert solver uses **bisection** rather than the more efficient Newton-Raphson with Battin's variables.
- Porkchop plots assume **single-revolution** transfers only (no Type II/III/IV).
- Multi-flyby trajectories use **simple v-infinity matching**, not Lambert + flyby co-optimization.
- No B-plane targeting or interplanetary navigation analysis.

### CR3BP and Halo orbits
- Halo orbit computation is a **research prototype**. Needs independent validation against published families (e.g., NASA Goddard's GMAT or Howell's results).
- Differential correction uses finite differences rather than the State Transition Matrix (the STM implementation has a known bug).
- No stability analysis (no monodromy matrix eigenvalue computation).
- No invariant manifold computation.

### Mission simulation
- Mission templates are **simplified educational scenarios**. They are intended to demonstrate sequencing and delta-v accounting, **not** to produce operational mission designs.
- Assumes **impulsive burns** (no finite burn duration modeling).
- No launch vehicle constraints or launch site geometry.
- No navigation error or covariance propagation.
- No inclination targeting or plane change analysis.
- No capture strategy beyond simple circular orbit insertion.
- No aerobraking sequences.
- Lunar transfers use heliocentric Lambert (incorrect — should use patched-conic Earth-Moon).

### Uncertainty quantification
- Monte Carlo only perturbs **thrust, Isp, mass, and departure date**.
- No model uncertainty (perturbation force model errors).
- No initial state covariance propagation.
- Sensitivity analysis uses **one-at-a-time** method, not Sobol or other global methods.

### Software engineering
- Single repository, no formal release versioning.
- No continuous integration setup yet (planned).
- Limited test coverage outside the validated cases.
- No fuzzing or property-based testing.
- No formal API stability guarantees.

## Use cases this project is suitable for

- Learning orbital mechanics concepts
- Exploring trade-offs in transfer design
- Sandbox for experimenting with perturbation models
- Educational demonstrations of CR3BP dynamics
- Numerical method comparisons (Verlet vs RK vs DOPRI5)

## Use cases this project is NOT suitable for

- Real mission design or operational trajectory planning
- Spacecraft navigation or operations
- Publication-quality scientific results without independent validation
- Any safety-critical application
- Comparison with operational systems where every meter matters
