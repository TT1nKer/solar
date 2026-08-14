# Dynamic Collapse Track — Milestone 01: Spherical Dust Collapse with a Newton → Post-Newtonian → General-Relativistic Transition

Status: spec (not yet implemented)
Branch: codex/dynamic-collapse-bh

## Objective

Prove Solar's simulation capability on the cloud-to-black-hole collapse path,
starting from the mathematical core: pressureless (dust) spherical collapse.
The same physical trajectory is integrated at three increasing orders of
theory — Newtonian free-fall, first post-Newtonian (1PN) corrected free-fall,
and the exact Oppenheimer-Snyder (OS) general-relativistic solution — with a
smooth, quantitatively controlled transition between regimes, terminating in
genuine black-hole formation.

## Why dust first

- Every regime has an exact solution, so every stage carries a pointwise
  verification anchor (the project's standing methodology).
- The Newtonian and OS surface trajectories share the same cycloid form,
  which makes the Newton → GR bridge explicit and continuous.
- The 1PN equation of motion is the epsilon-expansion of the Schwarzschild
  geodesic equation, so the transition is "one trajectory, increasing order
  of accuracy", not a hand-off between unrelated models.
- Solar already owns the two hard engines this needs: the generic
  Hamiltonian geodesic integrator (time-dependent-metric-capable) and the
  1PN central-body correction ForceModel (gr_correction.h, currently
  spacecraft-oriented).

## Regime 1 — Newtonian spherical dust (exact)

Shell equation, with M(r) the mass enclosed by the shell (conserved per
shell until shell crossing):

    d^2 r / dt^2 = -G M(r) / r^2

Homogeneous ball of mass M and initial radius R0: cycloid solution

    r(t) = (R0/2) (1 + cos eta),   t = (eta + sin eta) sqrt(R0^3 / (8 G M))

Free-fall time:

    t_ff = pi sqrt(R0^3 / (8 G M)) = sqrt(3 pi / (32 G rho0))

Inhomogeneous radial profiles: each mass shell keeps its own cycloid with
its own enclosed mass and its own free-fall time until shells cross
(caustics). Shell crossing is detected, reported, and compared with the
Lagrangian caustic condition (two adjacent shells reaching equal radius).

Verification:
- per-shell cycloid match (relative error vs exact, arbitrary tolerance scan)
- global energy integral per shell: 1/2 rdot^2 - G M(r)/r = E0 = const
- mass conservation: enclosed-mass function invariant in time
- free-fall time vs analytic t_ff (any density normalization)
- shell-crossing time for two known profiles (e.g. inner-density peak)

## Regime 2 — 1PN corrected dust (validated against exact GR)

Radial 1PN equation of motion (test-particle limit about mass M; the
standard 1PN two-body acceleration reduces to this for v || n):

    a_r = -(G M / r^2) [ 1 - 4 G M / (r c^2) - 3 v^2 / c^2 ]

    (flag: coefficient verification is a code task — the claim to prove is
     that integrating this EOM reproduces the exact Schwarzschild radial
     geodesic to O(epsilon^2), epsilon = G M / (r c^2).)

Applicability criterion per shell:

    epsilon(r) = 2 G M(r) / (r c^2),    beta(r) = |v(r)| / c

Blending policy (smooth, per-shell):

- epsilon < eps_lo (~1e-4): pure Newtonian
- eps_lo <= epsilon <= eps_hi (~1e-2): a = a_N + a_1PN * smoothstep
- epsilon > eps_hi: hand the shell to Regime 3 (GR), matched in r, dr/dt
  and enclosed mass

Verification:
- weak-field: 1PN trajectory vs exact Schwarzschild radial geodesic
  (Solar's existing Schwarzschild integrator), error scaling ~ epsilon^2
- 1PN vs Newtonian deviation at fixed epsilon matches the analytic formula
- blended trajectory is C1-continuous across both switch boundaries
- energy-like PN conserved quantity (1PN energy integral) drift < gate

## Regime 3 — Oppenheimer-Snyder collapse (exact GR)

Interior: closed FLRW dust ball

    ds^2 = -d tau^2 + a(tau)^2 [ d chi^2 / (1 - chi^2) + chi^2 d Omega^2 ]

    a(tau) = (a_max/2) (1 + cos eta),   tau = (a_max/2) (eta + sin eta)

Exterior: Schwarzschild with matched mass M. Surface radius R(tau) = a(tau)
chi0 follows the same cycloid functional form as the Newtonian case.

Quantities to derive and verify in code (flagged as derivation tasks, with
the known bounds stated up front):

- event-horizon birth: the horizon first appears at the center BEFORE the
  surface crosses R = 2 G M / c^2 (exact birth time to be derived from the
  null-geodesic matching; MTW treatment)
- outgoing-null freeze-out: photons emitted from the surface with
  R > 2 G M / c^2 escape with redshift diverging as R -> 2 G M / c^2; the
  last escaping photon leaves when R = 2 G M / c^2 exactly
- observed light-curve tail: F(t) ~ exp(-t / tau) with tau ~ 4 G M / c^3
  (coefficient to be confirmed against the null-geodesic computation)

Implementation: OS interior as a Solar Metric implementation (time-dependent
FLRW), so the existing generic Hamiltonian integrator traces null and
timelike geodesics through the collapsing interior — the first
time-dependent metric in Solar, and the direct proof that the integrator
handles dynamical spacetimes.

Verification:
- interior geodesics vs analytic FLRW null/timelike solutions
- matched interior/exterior metric continuity at the surface
- horizon birth time and freeze-out redshift divergence
- light-curve exponential tail constant
- emitted-photon energy vs analytic redshift formula

## Black-hole-formation conditions (the "guaranteed BH" parameter window)

These are the quantitative gates the cloud-to-core initial conditions must
satisfy for the end state to be a black hole; each is a runtime assertion
in the final pipeline, not an assumption:

1. Core mass gate: final collapsing core mass >= ~3 M_sun (above the
   ~2.2 M_sun cold TOV limit), so no neutron-star bounce in the
   realistic-gas extension.
2. Angular-momentum gate: J <= G M^2 / c (Kerr extremal bound); the
   collapse path must keep the core's J well below the bound (disk
   formation carries the excess J outward).
3. Fragmentation gate: cooling time t_cool >= t_ff across the core
   (Gammie-type criterion), otherwise the core fragments into low-mass
   stars instead of one massive object. Satisfied by low metallicity
   (Z <= 1e-4 Z_sun, H2-cooling suppressed) or an optically thick core.
4. Cloud budget: total cloud mass ~1e4 M_sun, initial radius ~1 pc,
   supersonic turbulence injected with the Kritsuk power spectrum
   P(k) ~ k^-2 and Mach ~ 5-10 (log-normal density PDF) — the initial
   conditions of the later full-scale stage.

For THIS milestone (dust) conditions 1-2 reduce to: total mass M, initial
compactness 2 G M / (R0 c^2) << 1, J = 0 (spherical), and the collapse is
guaranteed to terminate at the OS horizon.

## End-to-end transition verification (the "gradual Newton → GR" claim)

- conserved quantities carried across regime switches: enclosed mass per
  shell, total mass, angular momentum (0 in milestone 01), and the shell
  energy (with the PN energy integral in Regime 2)
- continuity diagnostics at both switch boundaries: r, dr/dt, acceleration
  ratio |a_1PN|/|a_N| = O(epsilon) at eps_hi, and the blended trajectory
  vs exact OS at the hand-off differing by O(epsilon_hi^2)
- full-run outputs: surface-radius history, horizon birth time, freeze-out
  light curve, final black-hole mass vs initial cloud mass

## Deliverables

- include/solar/dynamics/spherical_collapse.h (shell models, blending
  policy, horizon tracker, surface emission, observables)
- OS interior as a Metric implementation (time-dependent FLRW)
- tests/relativity/test_dynamic_collapse_*.cpp (one executable per regime
  plus one end-to-end transition test)
- docs/validation/dynamic_collapse_01_dust.md with measured numbers
- RELATIVITY_STATUS.md update (new track, phase state)

## Out of scope for milestone 01 (explicit)

Pressure/gas dynamics, rotation, turbulence injection, 3D geometry,
radiation beyond the surface blackbody, magnetic fields, and any
non-spherical GR. Each is a later milestone with its own anchor.
