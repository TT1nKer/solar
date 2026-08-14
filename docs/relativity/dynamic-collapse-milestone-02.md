# Dynamic Collapse Track — Milestone 02: LTB Compact Stage (exact GR dust collapse)

Status: implemented and passing — LTBCollapse with the full OS analytic anchor set (cycloid, t_ff, simultaneous crunch, horizon formation, t_obs log divergence, redshift/luminosity tail) and the N-body -> LTB hand-off test; see the verification section for measured numbers
Branch: codex/dynamic-collapse-bh
Follows: milestone 01 (3D post-Newtonian collapse), which ends at the
eps_hi hand-off: the terminal core is quasi-spherical and low-angular-
momentum (J <= G M^2 / c), handed over in position, velocity, and
enclosed mass.

## Position statement

The strong-field terminal stage is exact GR. The Lemaître-Tolman-Bondi
(LTB) dust solution is the exact general-relativistic counterpart of the
Newtonian inhomogeneous-shell collapse of Regime A: comoving pressureless
shells, each following its own Friedmann-like equation with its own
enclosed mass and energy function. Its uniform-density limit is the
Oppenheimer-Snyder (OS) collapse. Every anchor here is a closed-form
analytic result — this stage is pure math verification, no numerics
guessing.

## Physics

Per-shell equation of motion (bound, released from rest at turnaround):

    (dR/dtau)^2 = 2 G M(r) / R + 2 E(r),   E(r) = - G M(r) / r

with R the areal radius of the shell labeled r, M(r) the enclosed mass,
tau the shell's proper time. Parametric solution (collapse angle
theta in [0, pi]):

    R = (r / 2) (1 + cos theta)
    tau = t_ff(r) (theta + sin theta) / pi,
    t_ff(r) = pi sqrt(r^3 / (8 G M(r)))

Free-fall time per shell; the singularity is R = 0 at tau = t_ff(r).
OS limit (uniform rho0): t_ff(r) = pi sqrt(3 / (32 pi G rho0)) is
independent of r — all shells crunch simultaneously, and the surface
follows the exact Newtonian cycloid in proper time (this is the
statement that the proper-time cycloid of Regime B survives the
transition to GR unchanged).

Apparent horizon: a shell is trapped when R = 2 G M(r). For OS the
outermost shell is trapped first (2M/R grows like r^2 at fixed collapse
angle); the trapped region propagates inward and encloses the center at
the singularity — no naked singularity.

Exterior observables (surface shell, mass M): Schwarzschild-coordinate
time of the surface

    dt/dtau = E_inf / (1 - 2 G M / R),   E_inf = sqrt(1 - 2 G M / R0)

radial null-ray arrival at a distant observer

    t_obs = t_e + (r_obs - R_e) + 2 G M ln((r_obs - 2 G M) / (R_e - 2 G M))

redshift 1 + z = dt_obs / dtau_e, and the luminosity of a constant
rest-frame emitter L_obs = L_em / (1 + z)^2. Near the horizon
t_obs ~ 3 (G M / c^3) ln(1/delta) with delta = R_e / (2 G M / c^2) - 1,
so L_obs decays exponentially in t_obs with e-fold time
3 G M / (2 c^3). (The full Ames-Thorne 3 sqrt(3) G M / c^3 tail
additionally folds in the rest-frame luminosity evolution of the
collapsing emitter; the (1+z)^-2 factor is the redshift-squared piece.)
The surface freezes on the horizon in coordinate/observer time while
crossing it in finite proper time — the "gradual Newton -> GR"
signature at its strongest.

## Verification (measured)

All numbers from tests/relativity/test_ltb_os.cpp and
test_ltb_handoff.cpp (10 M_sun ball, R0 = 200 G M / c^2, t_ff = 0.1548 s):

- DONE — OS surface: cycloid in proper time; independent RK4 integration
  of the shell EOM agrees to 6.9e-12 relative; collapse time t_ff matches
  pi sqrt(R0^3 / (8 G M)) to 1e-12. The proper-time cycloid is exactly
  the Regime A/B Newtonian cycloid — GR changes coordinate time, not the
  dust's proper-time dynamics.
- DONE — simultaneous crunch: singularity-time spread across shells = 0
- DONE — horizon formation: surface trapped at the analytic collapse
  angle (tau_h / t_ff = 0.9996), R = 2 G M / c^2 at that moment, trapped
  region grows from 0 shells (surface horizon) to 100% near the crunch
- DONE — t_obs log divergence: the measured divergence coefficient is
  3 G M / c^3 to 0.3% (time-dilation + null-ray terms), structure ratio
  within 0.3% of ln(1/delta)
- DONE — redshift/luminosity: the (1+z)^-2 tail fits
  exp(-2 c^3 t_obs / (3 G M)) to 0.15%
- DONE — hand-off: an LTB shell reconstructed from the measured N-body
  surface (R_s, v_s, M) recovers the initial cloud radius to 1.7%,
  passes through the hand-off state exactly, and completes the collapse
  with horizon formation at 0.974 t_ff and singularity at 0.975 t_ff.
  Position and velocity imply the same collapse clock to 3% (the outer
  shell carries finite-N force noise, ~5% in radius at N = 2048,
  decreasing with N). Per-shell reconstruction recovers the initial
  radial profile to 2.1% mean error over 15 radial bins.
- NOTE — the full 3D aspherical GR collapse remains out of scope
  (numerical relativity); LTB is exact for the quasi-spherical terminal
  core the black-hole-formation gates select.

## Deliverables

- include/solar/dynamics/ltb_collapse.h: LTBCollapse (per-shell
  trajectories, singularity/horizon times, shell-crossing monitor) and
  the OS surface observables (t_obs, redshift, luminosity tail)
- tests/relativity/test_ltb_os.cpp, tests/relativity/test_ltb_handoff.cpp
- docs/validation/dynamic_collapse_02_ltb.md with measured numbers
- RELATIVITY_STATUS.md update

## Out of scope

Gas pressure (TOV-modified terminal stage), angular momentum (Kerr
exterior), radiation back-reaction, and the full 3D non-spherical GR
collapse (numerical relativity — a separate long-horizon program).