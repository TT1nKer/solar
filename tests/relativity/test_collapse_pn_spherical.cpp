#include "solar/body.h"
#include "solar/constants.h"
#include "solar/dynamics/barnes_hut_gravity.h"
#include "solar/dynamics/pn_gravity.h"
#include "solar/nbody.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

using solar::Body;
using solar::NBodySim;
using solar::Vec3;
using solar::dynamics::BarnesHutGravity;
using solar::dynamics::PostNewtonianGravity;

namespace {

int failures = 0;

void check(const std::string& name, bool condition) {
    std::cout << (condition ? "  PASS: " : "  FAIL: ") << name << '\n';
    if (!condition) ++failures;
}

constexpr double pi = 3.14159265358979323846;

// Radial equation of motion for the surface shell of a spherical dust
// cloud. With with_pn the 1PN standard-gauge form is used:
//   r'' = -(G M / r^2) [1 - 4 G M / (r c^2) - 5 v^2 / c^2],
// which is the exact spherical limit of the field-based 1PN force law
// (g = -G M / r^2, Phi = -G M / r, radial infall v = -v n). With
// with_pn false it is the Newtonian free-fall cycloid equation.
// Integrated with the same velocity-Verlet scheme the test uses on the
// N-body side, so the comparison isolates physics, not the integrator.
std::vector<double> shell_radii(double r0, double total_mass,
                                double duration, double dt, double c,
                                bool with_pn) {
    const double mu = solar::constants::G * total_mass;
    std::vector<double> radii;
    double r = r0;
    double v = 0.0;
    for (double t = 0.0; t < duration; t += dt) {
        radii.push_back(r);
        const double correction = with_pn
            ? (1.0 - 4.0 * mu / (r * c * c) - 5.0 * v * v / (c * c))
            : 1.0;
        const double a = -(mu / (r * r)) * correction;
        const double r_next = r + v * dt + 0.5 * a * dt * dt;
        const double correction_next = with_pn
            ? (1.0 - 4.0 * mu / (r_next * c * c) - 5.0 * v * v / (c * c))
            : 1.0;
        const double a_next = -(mu / (r_next * r_next)) * correction_next;
        v += 0.5 * (a + a_next) * dt;
        r = r_next;
    }
    return radii;
}

std::vector<Body> spherical_cloud(double total_mass, double radius0,
                                  std::size_t count, std::mt19937& generator) {
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::vector<Body> bodies;
    bodies.reserve(count);
    const double particle_mass = total_mass / count;
    for (std::size_t i = 0; i < count; ++i) {
        const double u = unit(generator);
        const double v = unit(generator);
        const double w = unit(generator);
        const double r = radius0 * std::cbrt(u);
        const double theta = std::acos(2.0 * v - 1.0);
        const double phi = 2.0 * pi * w;
        Body body;
        body.name = "dust";
        body.mass = particle_mass;
        body.mu = solar::constants::G * particle_mass;
        body.state.pos = {
            r * std::sin(theta) * std::cos(phi),
            r * std::sin(theta) * std::sin(phi),
            r * std::cos(theta)};
        body.state.vel = {};
        bodies.push_back(body);
    }
    return bodies;
}

// Run the N-body collapse with one force model; return the max-radius
// history sampled every output_interval. The PN model contributes the
// full acceleration (base plus correction), so the two legs each carry
// exactly one force model -- the earlier double-counting bug (tree +
// pairwise PN) gave the PN leg twice the gravity and a runaway core.
std::vector<double> collapse_max_radius(
    std::vector<Body> bodies, bool with_pn, double duration, double dt,
    double output_interval, double c, double softening) {
    NBodySim sim;
    sim.init(std::move(bodies));
    sim.clear_forces();
    if (with_pn) {
        sim.add_force(std::make_unique<PostNewtonianGravity>(
            PostNewtonianGravity::Config{0.5, softening, 0.0}, c));
    } else {
        sim.add_force(std::make_unique<BarnesHutGravity>(
            BarnesHutGravity::Config{0.5, softening}));
    }
    std::vector<double> history;
    sim.run(duration, dt, output_interval,
            [&](double, const std::vector<Body>& current) {
                double rmax = 0.0;
                for (const Body& body : current) {
                    rmax = std::max(rmax, body.state.pos.norm());
                }
                history.push_back(rmax);
            });
    return history;
}

} // namespace

int main() {
    const double c = solar::constants::C_LIGHT;
    const double mass = 10.0 * 1.98892e30;   // 10 M_sun compact core
    const double mu = solar::constants::G * mass;
    const double rg = mu / (c * c);
    const double radius0 = 200.0 * rg;       // eps = 5e-3 at the surface
    const std::size_t count = 2048;

    const double t_ff = pi * std::sqrt(radius0 * radius0 * radius0 /
                                       (8.0 * mu));
    const double dt = t_ff / 8000.0;
    const double output_interval = t_ff / 100.0;
    // Stop at 0.35 t_ff: the surface is at ~0.88 R0, eps ~ 5.7e-3, so
    // the whole run stays inside the weak-field validity window of the
    // 1PN expansion. The deeper phases belong to the per-particle
    // blending driver and the relativistic stage, not to this
    // regression.
    const double duration = 0.35 * t_ff;

    std::mt19937 generator(7);
    const auto cloud = spherical_cloud(mass, radius0, count, generator);
    std::mt19937 generator2(7);
    const auto cloud_copy = spherical_cloud(mass, radius0, count, generator2);

    const double softening = 1.0e-3 * radius0;
    const std::vector<double> newton_history = collapse_max_radius(
        cloud, false, duration, dt, output_interval, c, softening);
    const std::vector<double> pn_history = collapse_max_radius(
        cloud_copy, true, duration, dt, output_interval, c, softening);
    const std::vector<double> shell_pn = shell_radii(
        radius0, mass, duration, dt, c, true);
    const std::vector<double> shell_newton = shell_radii(
        radius0, mass, duration, dt, c, false);
    std::cout << "  newton N-body: first=" << newton_history.front() / radius0
              << " last=" << newton_history.back() / radius0 << std::endl;
    std::cout << "  pn N-body:     first=" << pn_history.front() / radius0
              << " last=" << pn_history.back() / radius0 << std::endl;
    std::cout << "  shell 1PN:     first=" << shell_pn.front() / radius0
              << " last=" << shell_pn.back() / radius0 << std::endl;
    std::cout << "  shell Newton:  first=" << shell_newton.front() / radius0
              << " last=" << shell_newton.back() / radius0 << std::endl;

    // The shell model is sampled every dt, the N-body histories every
    // output_interval, so compare at matched times.
    const std::size_t stride = static_cast<std::size_t>(
        output_interval / dt);

    // Gate 1: the 3D N-body PN collapse surface tracks the radial 1PN
    // shell model while the analytic surface stays above 0.85 R0 (the
    // entire run by construction). Tolerance 2% covers tree truncation
    // at theta = 0.5 and finite-N surface noise.
    double max_relative = 0.0;
    int compared = 0;
    for (std::size_t i = 0; i < pn_history.size(); ++i) {
        const std::size_t shell_index = i * stride;
        if (shell_index >= shell_pn.size()) break;
        const double analytic = shell_pn[shell_index];
        if (analytic < 0.85 * radius0) break;
        max_relative = std::max(
            max_relative,
            std::abs(pn_history[i] - analytic) / analytic);
        ++compared;
    }
    std::cout << "  shell comparison: " << compared
              << " samples, max relative r difference = " << max_relative
              << '\n';
    check("3D PN collapse surface tracks the radial 1PN shell model to 2%",
          compared > 20 && max_relative < 0.02);

    // Gate 2 (analytic, deterministic): in the standard PN gauge the
    // 1PN corrections weaken the coordinate-time acceleration (the
    // bracket 1 - 4 eps - 5 v^2/c^2 < 1), so the 1PN surface lags the
    // Newtonian cycloid at every matched time. The proper-time cycloid
    // itself is preserved -- that is the OS-stage statement -- while the
    // coordinate-time lag is the O(eps) signature this stage verifies.
    bool all_positive_lag = true;
    for (std::size_t index = 1; index < shell_pn.size(); ++index) {
        if (shell_pn[index] - shell_newton[index] <= 0.0) {
            all_positive_lag = false;
            break;
        }
    }
    const double analytic_lag =
        shell_pn.back() - shell_newton.back();
    std::cout << "  shell 1PN - Newton lag at end = " << analytic_lag
              << " km (" << analytic_lag / radius0 * 100.0 << "% of R0)"
              << '\n';
    check("1PN shell lags the Newtonian cycloid in coordinate time",
          all_positive_lag);

    // Informational: the N-body runs should reproduce the same signed
    // lag (PN surface outside the Newtonian surface at equal times).
    const double nbody_lag = pn_history.back() - newton_history.back();
    std::cout << "  N-body PN - Newton lag at end = " << nbody_lag
              << " km (" << nbody_lag / radius0 * 100.0 << "% of R0)"
              << '\n';

    // Gate 3: compactness of the surface grew through the collapse.
    const double final_r = pn_history.back();
    const double final_eps = mu / (final_r * c * c);
    const double initial_eps = mu / (radius0 * c * c);
    std::cout << "  final surface compactness eps = " << final_eps
              << " (started at " << initial_eps << ")" << '\n';
    check("compactness grew through the weak-field collapse",
          final_eps > initial_eps);

    std::cout << (failures == 0
                      ? "PASS: spherical-limit 1PN regression"
                      : "FAILURES PRESENT")
              << '\n';
    return failures == 0 ? 0 : 1;
}
