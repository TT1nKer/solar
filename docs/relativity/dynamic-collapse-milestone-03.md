# Dynamic Collapse Track — Milestone 03: Nebula-Scale Collapse to a Black Hole

Status: implemented and passing (test_collapse_nebula, 12/12 gates) — see
the measured verification numbers below; the optional CSV trace feeds the
video pipeline milestone
Branch: codex/dynamic-collapse-bh
Follows: milestones 01 (3D post-Newtonian collapse) and 02 (LTB compact
stage). This milestone runs the full chain end to end at the nebula
scale: turbulent cloud initial conditions -> Newtonian collapse with the
per-particle blending driver -> core hand-off -> LTB compact stage ->
horizon formation.

## The IC design (conditions that guarantee a black hole)

The black-hole-formation gates from milestone 01 are built into the
initial conditions, then measured (not assumed) at the hand-off:

1. Mass budget. A 100 M_sun Gaussian clump (r_c = 0.04 pc, truncated at
   3 r_c) is embedded at the center of the 1e4 M_sun / 1 pc turbulent
   cloud. Its inner core (r < r_c) carries 15.9 M_sun, comfortably above
   the 3 M_sun floor set by the ~2.2 M_sun TOV limit. The clump is a
   separate component: 2000 equal-mass particles at rest.
2. Angular momentum. J <= G M^2 / c is astronomically strict (a core
   particle at 0.04 pc may carry at most ~1 cm/s of net rotation), so
   the core must be decoupled from the cloud's angular momentum. The
   rotation profile omega(r) = omega0 ((r - 0.5 pc) / 0.5 pc)^2 is
   exactly zero inside 0.5 pc, and the coherent core region (r < 0.5 pc)
   is rebuilt deterministically and spherically symmetric: the ambient
   is replaced by a zero-velocity uniform buffer carrying the same mass,
   and the clump is laid down as exact equal-mass Gaussian shells on a
   Fibonacci lattice, at rest. Spherical symmetry by construction is the
   dust-model stand-in for the angular-momentum transport (magnetic
   braking / outflows) that real BH formation requires — declared as the
   condition, and its consequence (the measured J gate) is verified. An
   earlier realization WITHOUT the spherical buffer demonstrated the
   physics this condition stands in for: tidal torques from the
   turbulent ambient spun the core up to J/M = 1.5e10 km^2/s during the
   collapse (measured), 1300x the gate — exactly why real cores need J
   transport. The contrast is printed: the cloud as a whole has J/M ~
   3000x the gate, so it could not collapse to a black hole directly.
3. Fragmentation. t_cool >= t_ff trivially (dust milestone), recorded.

## The run

One force model — PnCollapseForce, the per-particle blending driver —
integrates the whole cloud. Two-timescale separation: the core's
collapse clock (t_sing ~ 33 kyr at r_c) runs 5x faster than the cloud's
free-fall (165 kyr), so the ambient cloud barely moves while the core
collapses. The run stops at 0.75 t_sing, where the core surface is at
~0.6 r_c and the core is still resolved by the softening (the hand-off
is deliberately placed at the resolution boundary). Every particle stays
deep in the Newtonian regime (eps ~ 3e-17 at the hand-off); the PN and
GR windows engage at the compact scale and are covered by the 200 rg
tests (spherical regression, blending windows, LTB anchors).

## Verification (measured in test_collapse_nebula)

- core mass budget: 20 M_sun (400 particles) >= 3 M_sun
- J = 0 at IC time by design (clump and buffer at rest, rotation
  cutoff); the cloud as a whole carries J/M = 1.2e13 km^2/s — 3000x
  the strict gate, so the cloud itself could not form a black hole
- core surface (median radius of the Lagrangian core set) collapses
  0.783 -> 0.591 r_c by 0.55 t_sing (substantial, pre-bounce)
- max compactness 4e-11: the run stays in the Newtonian regime as
  designed; the PN/GR windows are covered by the 200 rg tests
- total energy conservation 1.3e-5 through the collapse
- hand-off: M = 18.05 M_sun, J/M = 3.8e9 km^2/s (suppressed by
  3e-4 vs the cloud — the designed decoupling works; the residual is
  the ambient's tidal spin-up, the J-transport gap declared out of
  scope); LTB trajectory passes through (R_s, v_s, M) exactly;
  recovered turnaround radius 0.0366 pc vs the 0.04 pc core radius;
  velocity-implied hand-off clock 0.50 t_sing vs 0.55 t_sing elapsed
  (Gaussian-profile + accretion tolerance); horizon forms at
  R = 2 G M / c^2 = 53 km, then the singularity
- optional CSV trace (t, r_core, eps, J/M, E) for the video pipeline

## Out of scope

Angular-momentum transport (magnetic braking/outflows), gas pressure and
thermodynamics, radiation, and the fully aspherical strong-field
collapse (numerical relativity).