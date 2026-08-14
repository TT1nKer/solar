#include "solar/dynamics/ltb_collapse.h"

#include "solar/constants.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace solar {
namespace dynamics {

namespace {

constexpr double pi = 3.14159265358979323846;

// Invert tau = t_ff (theta + sin theta) / pi for theta in [0, pi].
// Newton iteration with the end-point asymptotic theta = pi - cbrt(6 s)
// where s = pi (1 - tau / t_ff), since 1 + cos theta vanishes there.
double collapse_angle(double tau, double t_ff) {
    if (tau <= 0.0) return 0.0;
    const double s = pi * tau / t_ff;
    if (s >= pi) return pi;
    double theta = (s <= 2.0) ? 0.5 * s : std::cbrt(3.0 * s);  // theta+sin(theta) ~ 2 theta for small s
    if (s > 2.5) {
        theta = pi - std::cbrt(6.0 * (pi - s));
    }
    for (int iteration = 0; iteration < 40; ++iteration) {
        const double f = theta + std::sin(theta) - s;
        const double fp = 1.0 + std::cos(theta);
        if (std::fabs(f) < 1.0e-14) break;
        if (fp < 1.0e-8) {  // near the end point, use the asymptotic
            theta = pi - std::cbrt(6.0 * (pi - s));
            continue;
        }
        const double step = f / fp;
        theta -= step;
        if (std::fabs(step) < 1.0e-14) break;
    }
    if (theta > pi) theta = pi;
    if (theta < 0.0) theta = 0.0;
    return theta;
}

} // namespace

LTBCollapse::LTBCollapse(std::vector<Shell> shells)
    : shells_(std::move(shells)) {}

double LTBCollapse::shell_radius(std::size_t i, double tau) const {
    const Shell& shell = shells_[i];
    const double gm = constants::G * shell.mass_enclosed_kg;  // km^3/s^2
    if (shell.energy_km2_s2 >= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double r_max = gm / (-shell.energy_km2_s2);  // km
    const double t_ff = singularity_time(i);
    const double theta = collapse_angle(tau, t_ff);
    return 0.5 * r_max * (1.0 + std::cos(theta));
}

double LTBCollapse::singularity_time(std::size_t i) const {
    const Shell& shell = shells_[i];
    const double gm = constants::G * shell.mass_enclosed_kg;
    const double e2 = -2.0 * shell.energy_km2_s2;  // km^2/s^2
    // Parametric time scale: tau = (G M / (-2 E)^{3/2}) (eta - sin eta).
    return pi * gm / std::pow(e2, 1.5);
}

double LTBCollapse::horizon_time(std::size_t i) const {
    const Shell& shell = shells_[i];
    const double gm = constants::G * shell.mass_enclosed_kg;
    const double r_g2 = 2.0 * gm / (constants::C_LIGHT * constants::C_LIGHT);
    const double r_max = gm / (-shell.energy_km2_s2);
    if (r_max <= r_g2) return std::numeric_limits<double>::infinity();
    // R = (R_max / 2) (1 + cos theta) = r_g2.
    const double cos_theta = 2.0 * r_g2 / r_max - 1.0;
    const double theta = std::acos(std::min(1.0, std::max(-1.0, cos_theta)));
    return singularity_time(i) * (theta + std::sin(theta)) / pi;
}

bool LTBCollapse::well_behaved() const {
    double previous = -1.0;
    for (std::size_t i = 0; i < shells_.size(); ++i) {
        const double t_sing = singularity_time(i);
        if (t_sing < previous * (1.0 - 1.0e-12)) return false;
        previous = t_sing;
        const Shell& shell = shells_[i];
        const double gm = constants::G * shell.mass_enclosed_kg;
        const double r_g2 =
            2.0 * gm / (constants::C_LIGHT * constants::C_LIGHT);
        if (r_g2 >= shell.radius_km) return false;
    }
    return true;
}

namespace {

LTBCollapse os_single_shell(double radius0_km, double mass_kg) {
    return LTBCollapse({LTBCollapse::Shell{
        radius0_km, mass_kg,
        -constants::G * mass_kg / radius0_km}});
}

} // namespace

double os_collapse_time(double radius0_km, double mass_kg) {
    return os_single_shell(radius0_km, mass_kg).singularity_time(0);
}

double os_surface_radius(double radius0_km, double mass_kg, double tau_s) {
    return os_single_shell(radius0_km, mass_kg).shell_radius(0, tau_s);
}

double os_horizon_time(double radius0_km, double mass_kg) {
    return os_single_shell(radius0_km, mass_kg).horizon_time(0);
}

double os_observed_time(double radius0_km, double mass_kg, double tau_s,
                        double observer_radius_km) {
    const double gm = constants::G * mass_kg;               // km^3/s^2
    const double c = constants::C_LIGHT;                    // km/s
    const double r_g2 = 2.0 * gm / (c * c);                 // km
    const double r_e = os_surface_radius(radius0_km, mass_kg, tau_s);
    if (r_e <= r_g2) {
        return std::numeric_limits<double>::infinity();
    }
    const double e_inf = std::sqrt(1.0 - r_g2 / radius0_km);
    const double t_ff = os_collapse_time(radius0_km, mass_kg);
    const double theta_e = collapse_angle(tau_s, t_ff);
    // Integrate dt/dtau = E_inf / (1 - r_g2 / R) numerically down to
    // delta_* = 1e-2; beyond that the integrand follows the horizon
    // asymptotic E_inf / (2 delta) with d(delta)/dtau = |R_dot| / r_g2,
    // which integrates in closed form (the log tail).
    const double delta_star = 1.0e-2;
    const double theta_star = std::acos(std::min(
        1.0, 2.0 * r_g2 * (1.0 + delta_star) / radius0_km - 1.0));
    double t_e = 0.0;
    if (theta_e > theta_star) {
        const double tau_star = t_ff * (theta_star + std::sin(theta_star)) / pi;
        const std::size_t steps = 65536;
        const double h = tau_star / static_cast<double>(steps);
        for (std::size_t k = 0; k <= steps; ++k) {
            const double tau_k = k * h;
            const double r_k =
                os_surface_radius(radius0_km, mass_kg, tau_k);
            const double f = e_inf / (1.0 - r_g2 / r_k);
            if (k == 0 || k == steps) {
                t_e += f;
            } else {
                t_e += (k % 2 == 1) ? 4.0 * f : 2.0 * f;
            }
        }
        t_e *= h / 3.0;
        // Tail: R = r_g2 (1 + delta), |R_dot| = c E_inf at the horizon
        // plus O(delta) corrections.
        const double delta_e = r_e / r_g2 - 1.0;
        const double r_dot_star =
            -(0.5 * radius0_km * std::sin(theta_e)) /
            ((t_ff / pi) * (1.0 + std::cos(theta_e)));
        // Integrand dt/dtau = E_inf / (2 delta) with d(delta)/dtau =
        // |R_dot| / r_g2: integral = E_inf r_g2 ln(...) / (2 |R_dot|).
        t_e += 0.5 * e_inf * (r_g2 / (-r_dot_star)) *
               std::log(delta_star / delta_e);
    } else {
        const std::size_t steps = 65536;
        const double h = tau_s / static_cast<double>(steps);
        for (std::size_t k = 0; k <= steps; ++k) {
            const double tau_k = k * h;
            const double r_k =
                os_surface_radius(radius0_km, mass_kg, tau_k);
            const double f = e_inf / (1.0 - r_g2 / r_k);
            if (k == 0 || k == steps) {
                t_e += f;
            } else {
                t_e += (k % 2 == 1) ? 4.0 * f : 2.0 * f;
            }
        }
        t_e *= h / 3.0;
    }
    // The null-ray terms are lengths; convert to seconds via /c.
    return t_e + (observer_radius_km - r_e) / c +
           (r_g2 / c) *
               std::log((observer_radius_km - r_g2) / (r_e - r_g2));
}

double os_surface_redshift(double radius0_km, double mass_kg, double tau_s,
                           double observer_radius_km) {
    (void)observer_radius_km;  // closed form is observer-independent
    const double gm = constants::G * mass_kg;
    const double c = constants::C_LIGHT;
    const double r_g2 = 2.0 * gm / (c * c);
    const double e_inf = std::sqrt(1.0 - r_g2 / radius0_km);
    const double r_e = os_surface_radius(radius0_km, mass_kg, tau_s);
    if (r_e <= r_g2) return std::numeric_limits<double>::infinity();
    const LTBCollapse model = os_single_shell(radius0_km, mass_kg);
    const double t_ff = model.singularity_time(0);
    const double theta = collapse_angle(tau_s, t_ff);
    const double r_dot =
        -(0.5 * radius0_km * std::sin(theta)) /
        ((t_ff / pi) * (1.0 + std::cos(theta)));  // km/s, negative
    const double dt_dtau = e_inf / (1.0 - r_g2 / r_e);
    const double dtobs_dre = -1.0 / c - r_g2 / (c * (r_e - r_g2));
    return dt_dtau + dtobs_dre * r_dot;  // 1 + z
}

double os_luminosity(double radius0_km, double mass_kg, double tau_s,
                     double observer_radius_km) {
    const double one_plus_z = os_surface_redshift(
        radius0_km, mass_kg, tau_s, observer_radius_km);
    return 1.0 / (one_plus_z * one_plus_z);
}

} // namespace dynamics
} // namespace solar