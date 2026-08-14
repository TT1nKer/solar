# Dynamic Collapse Track — Milestone 01 (revised): 3D Post-Newtonian Collapse of a Turbulent Rotating Cloud

Status: Regime A complete — Barnes-Hut tree gravity, the free-fall cycloid anchor, and the turbulent-cloud initial-condition generator (Kritsuk spectrum, log-normal density, rotation) are implemented and passing. Regime B in progress: the field-based 1PN force is in with its analytic anchors (c -> infinity reduction, radial spot check, perihelion precession) and the spherical-limit collapse regression passing; remaining is the per-particle blending driver and the eps_hi hand-off declaration.
Branch: codex/dynamic-collapse-bh
Supersedes: the spherical-dust draft of the same milestone name.

## Position statement

Spherical symmetry and uniform density are textbook limits, not the model.
The real collapse is three-dimensional, turbulent, rotating, and
multi-mass. This milestone makes the 3D dynamics the primary object and
demotes the spherical solutions to the role they deserve: exact special
cases used as regression anchors, nothing more.

## Honest physics map

- Fully aspherical collapse inside the strong-field GR regime requires
  numerical relativity; that is a separate, long-horizon program and is
  NOT claimed here.
- Asphericity, rotation, turbulence, fragmentation, and disk formation all
  dominate the WEAK-FIELD and TRANSITION regimes — which is exactly where
  the post-Newtonian expansion is a systematic, error-controlled theory.
- The black-hole-formation gates (below) force the terminal core toward
  low angular momentum and quasi-sphericity: this is a physical
  consequence (a high-J core fragments or disks out instead of collapsing
  to a black hole), not a modeling shortcut.

## Milestone 01 — 3D post-Newtonian collapse

### Regime A: Newtonian 3D self-gravitating collapse (exact anchors)

- Tree-code / FMM gravity, multi-mass particles, per-particle adaptive
  time steps, GPU-friendly traversal.
- Initial conditions: supersonic turbulence injection with the Kritsuk
  power spectrum P(k) ~ k^-2 (Mach 5-10), log-normal density PDF, a
  clump mass function (log-normal core + power-law tail), and net
  rotation (beta = E_rot / |E_grav| ~ 1-5%).
- Boundary: open (no periodic box — collapse, not driven-turbulence box).

Anchors (spherical special cases only):
- homogeneous ball free-fall: per-shell cycloid and t_ff = sqrt(3 pi / (32 G rho0))
- radially inhomogeneous dust: per-shell cycloid with per-shell enclosed mass
- energy/angular-momentum/mass conservation gates per particle and global
- N-doubling convergence on the density peak evolution

### Regime B: post-Newtonian corrections (the gradual Newton -> GR path)

Full 3D 1PN accelerations in the standard PN gauge, evaluated per
particle through the Newtonian field g and potential Phi of the other
bodies (the same Barnes-Hut octree walk as Regime A):

    a_i = g_i [1 + 4 Phi_i / c^2 - v_i^2 / c^2] - 4 (g_i . v_i) v_i / c^2

with per-particle compactness and velocity

    eps_i = -Phi_i / c^2 = G M_enc / (r_i c^2),   beta_i = v_i / c.

This field form is exactly the classical test-particle acceleration for
a single source (all analytic anchors carry over), and it reproduces the
mean-field spherical limit r'' = -(G M / r^2) [1 - 4 G M / (r c^2) -
5 v^2 / c^2] for an extended cloud. The pairwise test-particle form was
rejected on exactly this point: its static term scales per pair and
misses the enclosed-mass O(G M / (r c^2)) correction, so it cannot
recover the spherical limit (documented in commit dae5598). 2PN (with
the true cross-body terms of the Einstein-Infeld-Hoffmann Lagrangian)
is a follow-up layer on the same field formulation.

Blending policy (per particle, smooth):
- eps < eps_lo (~1e-4): pure Newtonian
- eps_lo <= eps <= eps_hi (~0.05): Newton + 1PN (+2PN near the top),
  smoothstep-weighted
- eps > eps_hi: terminal core hand-off to the compact-stage milestone
  (LTB), matched in position, velocity, and enclosed mass

Verification (status):
- DONE — c -> infinity limit: PN force -> Newtonian, deviation ~1e-15
  (tests/relativity/test_pn_gravity.cpp)
- DONE — radial spot check against the standard-gauge closed form at
  machine precision (same test)
- DONE — perihelion precession: 6 pi G M / (a (1-e^2) c^2) per orbit,
  measured 0.0100596 vs 0.0103569 rad/orbit (2.9%, Verlet discretization)
- DROPPED — 1PN conserved energy integral: the standard-gauge acceleration
  is not the gradient of the harmonic Lagrangian, so no exact first
  integral exists in this form; energy diagnostics stay Newtonian
- DONE — SPHERICAL-LIMIT REGRESSION: a 2048-particle spherical cloud run
  through the full 3D PN pipeline tracks the radial 1PN shell model to
  0.94% (gate 2%) through the weak-field phase (0.35 t_ff, eps 5e-3 ->
  5.4e-3), reproduces the signed standard-gauge coordinate-time lag
  behind the Newtonian cycloid (0.135% of R0 measured vs 0.167%
  analytic), and shows growing surface compactness
  (tests/relativity/test_collapse_pn_spherical.cpp)
- PENDING — at eps_hi, match the OS/LTB surface trajectory to the
  declared PN error order: that is the milestone-02 hand-off test

### Regime C (next milestone, declared here for continuity): LTB compact stage

Exact GR collapse with an ARBITRARY radial density profile
(Lemaître-Tolman-Bondi dust), i.e. the non-uniform general-relativistic
counterpart of Regime A's inhomogeneous shells. Horizon formation,
freeze-out redshift divergence, and the exponential light-curve tail are
derived and verified there. The 3D PN core hands off to LTB once its
compactness crosses eps_hi; the hand-off is justified by the
angular-momentum gate (J << G M^2 / c).

## Black-hole-formation gates (runtime assertions)

1. Core mass >= ~3 M_sun (above the ~2.2 M_sun TOV limit) — no neutron-star
   bounce when gas pressure is added in a later milestone.
2. Angular momentum J <= G M^2 / c at hand-off, and J tracked per shell.
3. Fragmentation suppression: t_cool >= t_ff across the core (low
   metallicity Z <= 1e-4 Z_sun or optically thick core) — enforced in the
   gas milestone; for dust milestone 01 this gate is trivially satisfied
   and recorded.
4. Cloud budget: total mass ~1e4 M_sun, initial radius ~1 pc, turbulent
   spectrum as above.

## Deliverables

- include/solar/dynamics/pn_collapse.h: 3D PN force laws, per-particle
  blending policy, compactness diagnostics, conservation monitors
- tree/FMM gravity core for multi-mass particles (Solar nbody extension)
- tests/relativity/test_dynamic_collapse_*.cpp: one executable per anchor
  (Newtonian special cases, PN error scaling, spherical-limit regression)
- docs/validation/dynamic_collapse_01_pn.md with measured numbers
- RELATIVITY_STATUS.md update (new track, phase state)

## Out of scope for milestone 01 (explicit)

Gas pressure and thermodynamics, magnetic fields, radiation transport,
full numerical relativity, and the LTB compact stage itself (milestone 02).