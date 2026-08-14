#pragma once

#include <vector>

namespace solar {
namespace dynamics {

// Lemaître-Tolman-Bondi (LTB) dust collapse: the exact general-
// relativistic terminal stage of the collapse pipeline (milestone 02).
// Comoving pressureless shells, each following its own Friedmann-like
// equation with its own enclosed mass and energy:
//
//   (dR/dtau)^2 = 2 G M(r) / R + 2 E(r)        (E < 0: bound)
//
// R is the shell's areal radius, tau its proper time, r the label
// (areal radius at the start of the collapse / turnaround), M(r) the
// enclosed mass, E(r) the energy per unit mass (km^2/s^2). A shell
// released from rest at R = r has E = - G M(r) / r. The parametric
// solution in the collapse angle theta in [0, pi]:
//
//   R = (R_max / 2) (1 + cos theta),   R_max = G M / (-E)
//   tau = t_ff (theta + sin theta) / pi,
//   t_ff = pi sqrt((G M) / (-2 E)^3) = pi sqrt(r^3 / (8 G M)) from rest
//
// The singularity is R = 0 at tau = t_ff. The uniform-density limit is
// the Oppenheimer-Snyder collapse: t_ff(r) is independent of r, all
// shells crunch simultaneously, and the surface follows the exact
// Newtonian cycloid in proper time.
class LTBCollapse {
public:
    struct Shell {
        double radius_km = 0.0;        // label r (start-of-collapse radius)
        double mass_enclosed_kg = 0.0; // M(r)
        double energy_km2_s2 = 0.0;    // E(r) < 0 for bound collapse
    };

    explicit LTBCollapse(std::vector<Shell> shells);

    // Areal radius of shell i at proper time tau (s from turnaround).
    // Valid for tau in [0, singularity_time(i)].
    double shell_radius(std::size_t i, double tau) const;

    // Proper-time free-fall to the singularity (s).
    double singularity_time(std::size_t i) const;

    // Proper time at which shell i crosses its trapped-surface radius
    // R = 2 G M / c^2; +inf when the shell starts inside the horizon
    // or never forms one. Requires the shell's own mass in the
    // condition, matching the LTB apparent horizon.
    double horizon_time(std::size_t i) const;

    // Sanity of the model: singularity times non-decreasing in r (no
    // shell crossing before the crunch) and every shell starts outside
    // its trapped radius.
    bool well_behaved() const;

private:
    std::vector<Shell> shells_;
};

// --- Oppenheimer-Snyder (uniform density) surface observables --------
// The surface shell of a uniform ball: cycloid in proper time, frozen
// star in observer time.

double os_collapse_time(double radius0_km, double mass_kg);

double os_surface_radius(double radius0_km, double mass_kg, double tau_s);

double os_horizon_time(double radius0_km, double mass_kg);

// Observer time of arrival of a photon emitted from the surface at
// proper time tau_s (must be before horizon crossing):
//   t_obs = t_e + (r_obs - R_e) + (2 G M / c^2)
//           ln((r_obs - 2 G M / c^2) / (R_e - 2 G M / c^2))
// diverges logarithmically as R_e -> 2 G M / c^2.
double os_observed_time(double radius0_km, double mass_kg, double tau_s,
                        double observer_radius_km);

// Surface redshift 1 + z = dt_obs / dtau_e (closed form).
double os_surface_redshift(double radius0_km, double mass_kg, double tau_s,
                           double observer_radius_km);

// Observed luminosity of the surface for unit emitted luminosity:
// L_obs = (1 + z)^-2, decaying exponentially in t_obs with e-fold time
// 3 sqrt(3) G M / c^2 ... measured in units where that statement holds
// at late observer times (Ames-Thorne tail).
double os_luminosity(double radius0_km, double mass_kg, double tau_s,
                     double observer_radius_km);

} // namespace dynamics
} // namespace solar
