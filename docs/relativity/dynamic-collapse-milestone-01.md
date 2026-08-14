# Dynamic Collapse Track — Milestone 01 (revised): 3D Post-Newtonian Collapse of a Turbulent Rotating Cloud

Status: Regime A in progress — Barnes-Hut tree gravity, the free-fall cycloid anchor, and the turbulent-cloud initial-condition generator (Kritsuk spectrum, log-normal density, rotation) are implemented and passing
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

Full 3D 1PN (then 2PN) N-body accelerations from the standard
post-Newtonian expansion (Blanchet/Will forms): velocity-dependent and
cross-body terms, per pair, with per-particle compactness

    eps_i = max_j [ G m_j / (r_ij c^2) ],   beta_i = v_i / c

Blending policy (per particle, smooth):
- eps < eps_lo (~1e-4): pure Newtonian
- eps_lo <= eps <= eps_hi (~0.05): Newton + 1PN (+2PN near the top),
  smoothstep-weighted
- eps > eps_hi: terminal core hand-off to the compact-stage milestone
  (LTB), matched in position, velocity, and enclosed mass

Verification:
- c -> infinity limit: PN trajectory -> Newtonian, deviation O(eps)
- two-body weak field: 1PN/2PN vs exact Schwarzschild geodesic (Solar
  already integrates these), error scaling ~ eps^2 / eps^3
- perihelion precession: Mercury 43 arcsec/century (existing Solar
  gr_correction cross-check, now as a two-body PN property)
- 1PN conserved energy integral drift < gate
- SPHERICAL-LIMIT REGRESSION: a spherical dust cloud run through the full
  3D PN pipeline must reproduce the radial 1PN EOM and, at eps_hi, match
  the OS/LTB surface trajectory to the declared PN error order

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
