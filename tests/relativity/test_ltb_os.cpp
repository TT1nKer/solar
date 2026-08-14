#include "solar/constants.h"
#include "solar/dynamics/ltb_collapse.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

using solar::dynamics::LTBCollapse;

namespace {

int failures = 0;

void check(const std::string& name, bool condition) {
    std::cout << (condition ? "  PASS: " : "  FAIL: ") << name << '\n';
    if (!condition) ++failures;
}

constexpr double pi = 3.14159265358979323846;

} // namespace

int main() {
    const double c = solar::constants::C_LIGHT;
    const double mass = 10.0 * 1.98892e30;           // 10 M_sun
    const double gm = solar::constants::G * mass;    // km^3/s^2
    const double rg = gm / (c * c);                  // km
    const double radius0 = 200.0 * rg;               // compact ball
    const double t_ff = pi * std::sqrt(radius0 * radius0 * radius0 /
                                       (8.0 * gm));
    const double r_obs = 1.0e6 * radius0;

    // --- 1. collapse time ----------------------------------------------
    const double model_tff = solar::dynamics::os_collapse_time(radius0, mass);
    std::cout << "  t_ff = " << model_tff << " s (analytic " << t_ff
              << " s)" << '\n';
    check("OS collapse time matches pi sqrt(R0^3 / 8 G M)",
          std::abs(model_tff - t_ff) < 1.0e-12 * t_ff);

    // --- 2. surface trajectory vs independent ODE integration ----------
    // RK4 on the shell equation dR/dtau = -sqrt(2 G M / R - 2 G M / R0),
    // independent of the parametric solution.
    {
        double r_ode = radius0;
        double v_ode = 0.0;
        const double dt = t_ff / 4.0e6;
        double max_relative = 0.0;
        for (double tau = 0.0; tau < 0.99 * t_ff; tau += dt * 100.0) {
            const double r_model =
                solar::dynamics::os_surface_radius(radius0, mass, tau);
            max_relative = std::max(max_relative,
                                    std::abs(r_model - r_ode) / r_ode);
            // advance the ODE by 100 tiny RK4 steps to the next sample
            for (int step = 0; step < 100; ++step) {
                const double a = -gm / (r_ode * r_ode);
                const double r1 = r_ode + 0.5 * dt * v_ode;
                const double v1 = v_ode + 0.5 * dt * a;
                const double a1 = -gm / (r1 * r1);
                const double r2 = r_ode + 0.5 * dt * v1;
                const double v2 = v_ode + 0.5 * dt * a1;
                const double a2 = -gm / (r2 * r2);
                const double r3 = r_ode + dt * v2;
                const double v3 = v_ode + dt * a2;
                const double a3 = -gm / (r3 * r3);
                r_ode += dt / 6.0 * (v_ode + 2.0 * v1 + 2.0 * v2 + v3);
                v_ode += dt / 6.0 * (a + 2.0 * a1 + 2.0 * a2 + a3);
            }
        }
        std::cout << "  surface vs RK4: max relative difference "
                  << max_relative << '\n';
        check("surface trajectory matches the independently integrated EOM",
              max_relative < 1.0e-6);
    }

    // --- 3. OS simultaneity --------------------------------------------
    {
        const std::vector<double> fractions{0.5, 0.75, 1.0};
        std::vector<LTBCollapse::Shell> shells;
        for (const double fraction : fractions) {
            const double r = fraction * radius0;
            const double m = mass * fraction * fraction * fraction;
            shells.push_back(LTBCollapse::Shell{
                r, m, -solar::constants::G * m / r});
        }
        const LTBCollapse model(shells);
        double spread = 0.0;
        for (std::size_t i = 0; i < shells.size(); ++i) {
            spread = std::max(spread,
                              std::abs(model.singularity_time(i) - t_ff));
        }
        std::cout << "  singularity-time spread across shells = "
                  << spread << " s" << '\n';
        check("uniform density: all shells crunch simultaneously",
              spread < 1.0e-12 * t_ff);
        check("LTB model is well behaved", model.well_behaved());
    }

    // --- 4. horizon formation ------------------------------------------
    {
        const double cos_theta_h = 4.0 * gm / (c * c * radius0) - 1.0;
        const double theta_h = std::acos(cos_theta_h);
        const double tau_h_analytic =
            t_ff * (theta_h + std::sin(theta_h)) / pi;
        const double tau_h = solar::dynamics::os_horizon_time(radius0, mass);
        std::cout << "  horizon time = " << tau_h << " s (analytic "
                  << tau_h_analytic << " s, t_ff = " << t_ff << " s)"
                  << '\n';
        check("surface crosses R = 2 G M at the analytic collapse angle",
              std::abs(tau_h - tau_h_analytic) < 1.0e-12 * t_ff);
        const double r_at_horizon =
            solar::dynamics::os_surface_radius(radius0, mass, tau_h);
        check("surface radius at the horizon time is 2 G M / c^2",
              std::abs(r_at_horizon - 2.0 * rg) < 1.0e-9 * rg);
        check("horizon forms before the singularity", tau_h < t_ff);

        // Trapped-region growth in the 3-shell model: none trapped at the
        // surface horizon time, all trapped near the crunch.
        const std::vector<double> fractions{0.5, 0.75, 1.0};
        std::vector<LTBCollapse::Shell> shells;
        for (const double fraction : fractions) {
            const double r = fraction * radius0;
            const double m = mass * fraction * fraction * fraction;
            shells.push_back(LTBCollapse::Shell{
                r, m, -solar::constants::G * m / r});
        }
        const LTBCollapse model(shells);
        auto trapped = [&](double tau) {
            int count = 0;
            for (std::size_t i = 0; i < shells.size(); ++i) {
                const double r_g2_i = 2.0 * solar::constants::G *
                    shells[i].mass_enclosed_kg / (c * c);
                if (model.shell_radius(i, tau) < r_g2_i) ++count;
            }
            return count;
        };
        const int trapped_at_horizon = trapped(tau_h);
        const int trapped_at_crunch = trapped(t_ff * (1.0 - 1.0e-9));
        std::cout << "  trapped shells: " << trapped_at_horizon
                  << " at the surface horizon, " << trapped_at_crunch
                  << " near the crunch" << '\n';
        check("trapped region starts at the surface and reaches the center",
              trapped_at_horizon == 0 && trapped_at_crunch == 3);
    }

    // --- 5. observer-time log divergence -------------------------------
    {
        const double delta1 = 1.0e-2;
        const double delta2 = 1.0e-3;
        const double delta3 = 1.0e-4;
        auto observed = [&](double delta) {
            const double r_e = 2.0 * rg * (1.0 + delta);
            const double cos_theta = 2.0 * r_e / radius0 - 1.0;
            const double theta = std::acos(cos_theta);
            const double tau = t_ff * (theta + std::sin(theta)) / pi;
            return solar::dynamics::os_observed_time(
                radius0, mass, tau, r_obs);
        };
        const double t1 = observed(delta1);
        const double t2 = observed(delta2);
        const double t3 = observed(delta3);
        const double ratio = (t1 - t2) / (t2 - t3);
        const double expected_ratio =
            std::log(delta1 / delta2) / std::log(delta2 / delta3);
        std::cout << "  t_obs ratios: measured " << ratio << " vs "
                  << expected_ratio << " (log divergence)" << '\n';
        check("t_obs diverges logarithmically in R_e - 2 G M",
              std::abs(ratio - expected_ratio) < 0.015 * expected_ratio);
        // t1 < t2 (earlier emission, farther from the horizon), so the
        // slope is negative; its magnitude is the divergence coefficient.
        const double slope = (t1 - t2) / std::log(delta1 / delta2);
        const double expected_slope = -3.0 * gm / (c * c * c);
        std::cout << "  divergence slope = " << slope << " s vs "
                  << expected_slope << " s (-3 G M / c^3)" << '\n';
        check("log-divergence coefficient is 3 G M / c^3",
              std::abs(slope - expected_slope) < 0.03 * std::abs(expected_slope));
    }

    // --- 6. redshift and luminosity tail -------------------------------
    {
        // Fit ln L vs t_obs over late emission times; the (1+z)^-2 model
        // predicts ln L = const - 2 c^3 t_obs / (3 G M) + O(delta).
        std::vector<double> t_obs, ln_l;
        for (const double delta : {5.0e-3, 2.0e-3, 1.0e-3, 5.0e-4,
                                   2.0e-4, 1.0e-4}) {
            const double r_e = 2.0 * rg * (1.0 + delta);
            const double cos_theta = 2.0 * r_e / radius0 - 1.0;
            const double theta = std::acos(cos_theta);
            const double tau = t_ff * (theta + std::sin(theta)) / pi;
            t_obs.push_back(solar::dynamics::os_observed_time(
                radius0, mass, tau, r_obs));
            ln_l.push_back(std::log(solar::dynamics::os_luminosity(
                radius0, mass, tau, r_obs)));
        }
        // Linear least squares.
        const std::size_t n = t_obs.size();
        double mean_t = 0.0, mean_l = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            mean_t += t_obs[i];
            mean_l += ln_l[i];
        }
        mean_t /= static_cast<double>(n);
        mean_l /= static_cast<double>(n);
        double cov = 0.0, var = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            cov += (t_obs[i] - mean_t) * (ln_l[i] - mean_l);
            var += (t_obs[i] - mean_t) * (t_obs[i] - mean_t);
        }
        const double slope = cov / var;
        const double expected_slope = -2.0 * c * c * c / (3.0 * gm);
        std::cout << "  ln L slope = " << slope << " 1/s vs "
                  << expected_slope << " 1/s (-2 c^3 / 3 G M)" << '\n';
        check("luminosity tail decays as exp(-2 c^3 t_obs / (3 G M))",
              std::abs(slope - expected_slope) <
                  0.05 * std::abs(expected_slope));
    }

    std::cout << (failures == 0 ? "PASS: OS collapse anchors"
                                : "FAILURES PRESENT")
              << '\n';
    return failures == 0 ? 0 : 1;
}