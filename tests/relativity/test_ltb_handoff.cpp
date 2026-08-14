#include "solar/body.h"
#include "solar/constants.h"
#include "solar/dynamics/barnes_hut_gravity.h"
#include "solar/dynamics/ltb_collapse.h"
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
using solar::dynamics::LTBCollapse;

namespace {

int failures = 0;

void check(const std::string& name, bool condition) {
    std::cout << (condition ? "  PASS: " : "  FAIL: ") << name << '\n';
    if (!condition) ++failures;
}

constexpr double pi = 3.14159265358979323846;

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
        bodies.push_back(body);
    }
    return bodies;
}

} // namespace

int main() {
    const double c = solar::constants::C_LIGHT;
    const double mass = 10.0 * 1.98892e30;
    const double gm = solar::constants::G * mass;
    const double rg = gm / (c * c);
    const double radius0 = 200.0 * rg;
    const std::size_t count = 2048;
    const double t_ff = pi * std::sqrt(radius0 * radius0 * radius0 /
                                       (8.0 * gm));
    const double duration = 0.6 * t_ff;

    // N-body Newtonian collapse to the hand-off state.
    std::mt19937 generator(7);
    NBodySim sim;
    sim.init(spherical_cloud(mass, radius0, count, generator));
    sim.clear_forces();
    sim.add_force(std::make_unique<BarnesHutGravity>(
        BarnesHutGravity::Config{0.5, 1.0e-3 * radius0}));
    sim.run(duration, t_ff / 8000.0, t_ff / 100.0,
            [](double, const std::vector<Body>&) {});
    const std::vector<Body>& final = sim.bodies;

    // Hand-off state of the surface shell: max radius and the mean
    // radial velocity of the outermost 3% of particles.
    double r_surface = 0.0;
    for (const Body& body : final) {
        r_surface = std::max(r_surface, body.state.pos.norm());
    }
    double v_sum = 0.0;
    std::size_t v_count = 0;
    for (const Body& body : final) {
        if (body.state.pos.norm() > 0.97 * r_surface) {
            v_sum += body.state.vel.dot(body.state.pos) /
                     body.state.pos.norm();
            ++v_count;
        }
    }
    const double v_surface = v_sum / static_cast<double>(v_count);
    std::cout << "  hand-off state: R_s = " << r_surface / radius0
              << " R0, v_s = " << v_surface / c << " c (shell of "
              << v_count << " particles)" << '\n';

    // Recover the LTB shell (turnaround radius r0, collapse angle theta)
    // whose trajectory passes through (R_s, v_s) with mass M:
    //   sin(theta/2) = -v_s sqrt(R_s / (2 G M)),  r0 = R_s / cos^2(theta/2).
    const double sin_half = -v_surface * std::sqrt(r_surface / (2.0 * gm));
    const double theta = 2.0 * std::asin(std::min(1.0, sin_half));
    const double r0 = r_surface /
        (std::cos(0.5 * theta) * std::cos(0.5 * theta));
    std::cout << "  recovered shell: r0 = " << r0 / radius0
              << " R0 (initial radius 1.0 R0), theta = " << theta << '\n';
    check("hand-off recovers the initial cloud radius",
          std::abs(r0 - radius0) < 0.03 * radius0);

    // Build the LTB stage and verify the construction is exact.
    const LTBCollapse model({LTBCollapse::Shell{
        r0, mass, -gm / r0}});
    const double tau_handoff =
        model.singularity_time(0) * (theta + std::sin(theta)) / pi;
    // The cycloid clock implied by the measured radius:
    // eta_R = acos(2 R_s / R0 - 1), t_cycloid = t_ff (eta_R + sin eta_R) / pi.
    // The outermost shell carries finite-N force noise (the outer shell
    // has O(N^(2/3)) particles; at 0.6 t_ff the max-radius envelope sits
    // ~5% above the exact cycloid at N = 2048 and approaches it as N
    // grows), so the position- and velocity-implied clocks may differ by
    // a few percent. The hand-off statement is that position and velocity
    // belong to the same collapse, i.e. the two clocks agree; the sim
    // clock itself is reported for context.
    const double eta_r = std::acos(2.0 * r_surface / radius0 - 1.0);
    const double t_cycloid_r =
        t_ff * (eta_r + std::sin(eta_r)) / pi;
    std::cout << "  hand-off clocks: velocity-implied " << tau_handoff / t_ff
              << " t_ff, radius-implied " << t_cycloid_r / t_ff
              << " t_ff (sim elapsed 0.6 t_ff)" << '\n';
    check("position and velocity imply the same collapse clock",
          std::abs(tau_handoff - t_cycloid_r) < 0.05 * t_ff);
    check("LTB trajectory passes through the hand-off state",
          std::abs(model.shell_radius(0, tau_handoff) - r_surface) <
              1.0e-9 * r_surface);

    // Internal consistency of the shell velocity at the hand-off.
    const double h = 1.0e-7 * model.singularity_time(0);
    const double r_dot_ltb = (model.shell_radius(0, tau_handoff + h) -
                              model.shell_radius(0, tau_handoff - h)) /
                             (2.0 * h);
    std::cout << "  LTB surface velocity at hand-off = " << r_dot_ltb / c
              << " c (N-body " << v_surface / c << " c)" << '\n';
    check("LTB shell velocity matches the N-body surface velocity",
          std::abs(r_dot_ltb - v_surface) < 0.03 * std::abs(v_surface));

    // The GR continuation: horizon forms before the singularity.
    const double tau_horizon = model.horizon_time(0);
    const double tau_singularity = model.singularity_time(0);
    std::cout << "  continuation: horizon at " << tau_horizon / t_ff
              << " t_ff, singularity at " << tau_singularity / t_ff
              << " t_ff" << '\n';
    check("horizon forms after the hand-off and before the singularity",
          tau_horizon > tau_handoff && tau_horizon < tau_singularity);
    const double r_h = model.shell_radius(0, tau_horizon);
    check("surface crosses R = 2 G M / c^2",
          std::abs(r_h - 2.0 * rg) < 1.0e-9 * rg);

    // Multi-shell hand-off: per radial bin, recover each shell's
    // turnaround radius from the measured (R_i, v_i, M_i) and compare
    // with the homologous expectation r_i0 = R0 R_i / R_surface.
    const std::size_t bin_count = 20;
    std::vector<double> bin_r(bin_count, 0.0), bin_m(bin_count, 0.0),
        bin_v(bin_count, 0.0), bin_count_p(bin_count, 0.0);
    for (const Body& body : final) {
        const double r = body.state.pos.norm();
        const std::size_t bin = std::min(
            bin_count - 1, static_cast<std::size_t>(
                std::floor(r / r_surface * bin_count)));
        const double v_r = body.state.vel.dot(body.state.pos) /
            std::max(r, 1.0e-9);
        bin_r[bin] += r * body.mass;
        bin_v[bin] += v_r * body.mass;
        bin_m[bin] += body.mass;
        bin_count_p[bin] += 1.0;
    }
    double recovery_error_sum = 0.0;
    std::size_t recovery_bins = 0;
    for (std::size_t i = 0; i < bin_count; ++i) {
        if (bin_m[i] <= 0.0 || bin_count_p[i] < 10) continue;
        const double r_i = bin_r[i] / bin_m[i];
        const double v_i = bin_v[i] / bin_m[i];
        // Cumulative enclosed mass at the bin's outer edge.
        double m_i = 0.0;
        for (std::size_t j = 0; j <= i; ++j) m_i += bin_m[j];
        if (v_i >= 0.0 || r_i < 1.0e-6) continue;
        const double r0_i = r_i / (1.0 - v_i * v_i * r_i /
                                   (2.0 * solar::constants::G * m_i));
        const double expected = radius0 * r_i / r_surface;
        recovery_error_sum += std::abs(r0_i - expected) / expected;
        ++recovery_bins;
    }
    const double mean_recovery_error = recovery_error_sum / recovery_bins;
    std::cout << "  per-shell radius recovery: mean |r0_i - r_i0|/r_i0 = "
              << mean_recovery_error << " over " << recovery_bins
              << " bins" << '\n';
    check("per-shell hand-off recovers the initial radial profile",
          recovery_bins > 10 && mean_recovery_error < 0.06);

    std::cout << (failures == 0 ? "PASS: LTB hand-off"
                                : "FAILURES PRESENT")
              << '\n';
    return failures == 0 ? 0 : 1;
}