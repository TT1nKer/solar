#include "solar/body.h"
#include "solar/constants.h"
#include "solar/dynamics/barnes_hut_gravity.h"
#include "solar/dynamics/pn_collapse.h"
#include "solar/dynamics/pn_gravity.h"
#include "solar/nbody.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

using solar::Body;
using solar::NBodySim;
using solar::Vec3;
using solar::dynamics::BarnesHutGravity;
using solar::dynamics::PnCollapseForce;
using solar::dynamics::PostNewtonianGravity;

namespace {

int failures = 0;

void check(const std::string& name, bool condition) {
    std::cout << (condition ? "  PASS: " : "  FAIL: ") << name << '\n';
    if (!condition) ++failures;
}

constexpr double pi = 3.14159265358979323846;

Body make_body(const Vec3& pos, const Vec3& vel, double mass) {
    Body body;
    body.name = "dust";
    body.mass = mass;
    body.mu = solar::constants::G * mass;
    body.state.pos = pos;
    body.state.vel = vel;
    return body;
}

std::vector<Body> random_cloud(std::size_t count, std::mt19937& generator) {
    std::uniform_real_distribution<double> position(-1.0e6, 1.0e6);
    std::uniform_real_distribution<double> velocity(-50.0, 50.0);
    std::uniform_real_distribution<double> mass(1.0e24, 1.0e28);
    std::vector<Body> bodies;
    for (std::size_t i = 0; i < count; ++i) {
        bodies.push_back(make_body(
            {position(generator), position(generator), position(generator)},
            {velocity(generator), velocity(generator), velocity(generator)},
            mass(generator)));
    }
    return bodies;
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
        bodies.push_back(make_body(
            {r * std::sin(theta) * std::cos(phi),
             r * std::sin(theta) * std::sin(phi),
             r * std::cos(theta)},
            {}, particle_mass));
    }
    return bodies;
}

double relative_acc_deviation(const std::vector<Vec3>& a,
                              const std::vector<Vec3>& b) {
    double worst = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        worst = std::max(
            worst, (a[i] - b[i]).norm() / std::max(b[i].norm(), 1.0e-12));
    }
    return worst;
}

std::vector<double> collapse_max_radius(
    std::vector<Body> bodies, int mode, double duration, double dt,
    double output_interval, double softening) {
    NBodySim sim;
    sim.init(std::move(bodies));
    sim.clear_forces();
    switch (mode) {
    case 0:  // pure Newtonian tree gravity
        sim.add_force(std::make_unique<BarnesHutGravity>(
            BarnesHutGravity::Config{0.5, softening}));
        break;
    case 1:  // full 1PN (un-blended)
        sim.add_force(std::make_unique<PostNewtonianGravity>(
            PostNewtonianGravity::Config{0.5, softening, 0.0}, 0.0));
        break;
    default:  // blended Newton -> 1PN windows
        sim.add_force(std::make_unique<PnCollapseForce>(
            PnCollapseForce::Config{0.5, softening, 1.0e-4, 0.05, 0.0},
            0.0));
        break;
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
    const double softening = 1.0e3;  // km, for the random-cloud checks

    // --- 1. Degenerate window (eps_lo = eps_hi = 0) -> full 1PN --------
    {
        std::mt19937 generator(11);
        const auto bodies = random_cloud(64, generator);
        PnCollapseForce blend(
            PnCollapseForce::Config{0.5, softening, 0.0, 0.0, 0.0}, 0.0);
        PostNewtonianGravity pn(
            PostNewtonianGravity::Config{0.5, softening, 0.0}, 0.0);
        std::vector<Vec3> acc_blend(bodies.size(), Vec3{});
        std::vector<Vec3> acc_pn(bodies.size(), Vec3{});
        blend.compute(bodies, 0.0, acc_blend);
        pn.compute(bodies, 0.0, acc_pn);
        const double deviation = relative_acc_deviation(acc_blend, acc_pn);
        std::cout << "  full-1PN reduction max relative deviation = "
                  << deviation << '\n';
        check("eps_lo = eps_hi = 0 reduces to the un-blended 1PN force",
              deviation < 1.0e-12);
    }

    // --- 2. Large eps_lo -> pure Newtonian tree gravity -----------------
    {
        std::mt19937 generator(11);
        const auto bodies = random_cloud(64, generator);
        PnCollapseForce blend(
            PnCollapseForce::Config{0.5, softening, 1.0e9, 1.0e9, 0.0}, 0.0);
        BarnesHutGravity tree(
            BarnesHutGravity::Config{0.5, softening});
        std::vector<Vec3> acc_blend(bodies.size(), Vec3{});
        std::vector<Vec3> acc_tree(bodies.size(), Vec3{});
        blend.compute(bodies, 0.0, acc_blend);
        tree.compute(bodies, 0.0, acc_tree);
        const double deviation = relative_acc_deviation(acc_blend, acc_tree);
        std::cout << "  Newtonian reduction max relative deviation = "
                  << deviation << '\n';
        check("large eps_lo reduces to the Newtonian tree gravity",
              deviation < 1.0e-12);
    }

    // --- 3-4. Compactness diagnostic and hand-off tagging ---------------
    const double mass = 10.0 * 1.98892e30;
    const double mu = solar::constants::G * mass;
    const double rg = mu / (c * c);
    const double radius0 = 200.0 * rg;   // eps0 = 5e-3 at the surface
    const std::size_t count = 1024;
    const double eps0 = mu / (radius0 * c * c);

    std::mt19937 generator(7);
    const auto cloud = spherical_cloud(mass, radius0, count, generator);

    const std::vector<double> eps =
        solar::dynamics::per_particle_compactness(cloud, 1.0e-3 * radius0, c);

    // Inside a uniform sphere Phi(r) = -(G M / R0) (3/2 - r^2 / (2 R0^2)),
    // so eps(r) = eps0 (1.5 - 0.5 (r/R0)^2) ranges from eps0 at the
    // surface to 1.5 eps0 at the center. The single-particle direct sum
    // carries ~3% finite-N surface noise (the 1/r pairs), so the surface
    // anchor uses the shell mean and the center anchor allows the same
    // slack.
    double shell_sum = 0.0;
    std::size_t shell_count = 0;
    std::size_t innermost = 0;
    for (std::size_t i = 0; i < cloud.size(); ++i) {
        if (cloud[i].state.pos.norm() > 0.9 * radius0) {
            shell_sum += eps[i];
            ++shell_count;
        }
        if (cloud[i].state.pos.norm() <
            cloud[innermost].state.pos.norm()) {
            innermost = i;
        }
    }
    const double shell_mean = shell_sum / static_cast<double>(shell_count);
    // Shell particles sit at r in (0.9, 1.0] R0; for a uniform sphere
    // E[r^2/R0^2] = 0.6 (1 - 0.9^5) / (1 - 0.9^3) = 0.9067, so the mean
    // compactness there is eps0 (1.5 - 0.5 * 0.9067) ~ 1.047 eps0.
    const double shell_sq_mean = 0.6 * (1.0 - std::pow(0.9, 5)) /
                                 (1.0 - std::pow(0.9, 3));
    const double expected_shell =
        eps0 * (1.0 - 1.0 / count) * (1.5 - 0.5 * shell_sq_mean);
    std::cout << "  shell-mean eps = " << shell_mean << " expected "
              << expected_shell << " (shell " << shell_count
              << " particles)" << '\n';
    check("outer-shell mean compactness matches the surface value",
          std::abs(shell_mean - expected_shell) < 0.01 * eps0);

    std::cout << "  innermost eps = " << eps[innermost]
              << " (uniform-sphere center limit 1.5 eps0 = "
              << 1.5 * eps0 << ")" << '\n';
    check("center compactness sits at the 1.5 eps0 analytic maximum",
          eps[innermost] > 1.35 * eps0 && eps[innermost] < 1.65 * eps0);

    // With eps_hi = 0.05 nothing is compact enough at t=0 (the whole
    // cloud lives at eps <= 1.5 eps0 = 7.5e-3): the hand-off boundary is
    // a dynamical one, reached only deep in the collapse. With
    // eps_hi = 1.2 eps0 the tag must extract the central region
    // r <= sqrt(0.6) R0 ~ 0.775 R0 (about 46.5% of the particles).
    const std::vector<int> none_yet = solar::dynamics::handoff_candidates(
        cloud, 1.0e-3 * radius0, c, 0.05);
    check("no hand-off candidates at t=0 with eps_hi = 0.05",
          none_yet.empty());

    const std::vector<int> candidates = solar::dynamics::handoff_candidates(
        cloud, 1.0e-3 * radius0, c, 1.2 * eps0);
    // The analytic boundary r = sqrt(0.6) R0 encloses 46.5% of the
    // particles; individual-particle eps carries heavy-tailed 1/r pair
    // noise (~6% in radius at the boundary), so the anchors are the
    // tagged fraction and the mean radii of the two sets.
    double flagged_radius_sum = 0.0;
    for (const int index : candidates) {
        flagged_radius_sum +=
            cloud[static_cast<std::size_t>(index)].state.pos.norm();
    }
    double unflagged_radius_sum = 0.0;
    for (std::size_t i = 0; i < cloud.size(); ++i) {
        const bool flagged = std::find(candidates.begin(), candidates.end(),
                                       static_cast<int>(i)) != candidates.end();
        if (!flagged) unflagged_radius_sum += cloud[i].state.pos.norm();
    }
    const double flagged_mean =
        flagged_radius_sum / static_cast<double>(candidates.size());
    const double unflagged_mean = unflagged_radius_sum /
        static_cast<double>(count - candidates.size());
    const double analytic_fraction = std::pow(0.6, 1.5);
    std::cout << "  hand-off: " << candidates.size() << "/" << count
              << " tagged (analytic fraction " << analytic_fraction << ")"
              << ", mean r: flagged " << flagged_mean / radius0
              << " R0 vs unflagged " << unflagged_mean / radius0
              << " R0" << '\n';
    check("hand-off tagging extracts the central region (r ~ 0.775 R0)",
          candidates.size() > 350 && candidates.size() < 600 &&
              flagged_mean < 0.6 * radius0 &&
              unflagged_mean > 0.75 * radius0);

    // --- 5. Blended collapse dynamics -----------------------------------
    const double t_ff = pi * std::sqrt(radius0 * radius0 * radius0 /
                                       (8.0 * mu));
    const double dt = t_ff / 8000.0;
    const double output_interval = t_ff / 100.0;
    const double duration = 0.35 * t_ff;

    std::mt19937 generator2(7);
    const auto cloud_newton = spherical_cloud(mass, radius0, count, generator2);
    std::mt19937 generator3(7);
    const auto cloud_pn = spherical_cloud(mass, radius0, count, generator3);
    std::mt19937 generator4(7);
    const auto cloud_blend = spherical_cloud(mass, radius0, count, generator4);

    const auto newton_history = collapse_max_radius(
        cloud_newton, 0, duration, dt, output_interval, 1.0e-3 * radius0);
    const auto pn_history = collapse_max_radius(
        cloud_pn, 1, duration, dt, output_interval, 1.0e-3 * radius0);
    const auto blend_history = collapse_max_radius(
        cloud_blend, 2, duration, dt, output_interval, 1.0e-3 * radius0);

    const double final_newton = newton_history.back();
    const double final_pn = pn_history.back();
    const double final_blend = blend_history.back();
    std::cout << "  final surface radii (R0 units): newton "
              << final_newton / radius0 << ", blend "
              << final_blend / radius0 << ", pn " << final_pn / radius0
              << '\n';
    // The PN correction weakens coordinate-time acceleration, so the
    // full-PN surface lags the Newtonian one; the blend interpolates.
    const double slack = 2.0e-4 * radius0;
    check("blended surface lies between Newtonian and full-1PN surfaces",
          final_blend >= final_newton - slack &&
              final_blend <= final_pn + slack);

    // Total compactness sum grows as the cloud contracts
    // (sum eps_i = -2 PE / (M c^2) for equal masses).
    double eps_sum_initial = 0.0;
    for (const double value : eps) eps_sum_initial += value;
    NBodySim blend_sim;
    blend_sim.init(cloud_blend);
    blend_sim.clear_forces();
    blend_sim.add_force(std::make_unique<PnCollapseForce>(
        PnCollapseForce::Config{0.5, 1.0e-3 * radius0, 1.0e-4, 0.05, 0.0},
        0.0));
    blend_sim.run(duration, dt, output_interval,
                  [](double, const std::vector<Body>&) {});
    const std::vector<double> eps_final =
        solar::dynamics::per_particle_compactness(
            blend_sim.bodies, 1.0e-3 * radius0, c);
    double eps_sum_final = 0.0;
    for (const double value : eps_final) eps_sum_final += value;
    std::cout << "  sum eps_i: " << eps_sum_initial << " -> "
              << eps_sum_final << '\n';
    check("total compactness grows through the blended collapse",
          eps_sum_final > eps_sum_initial * 1.02);

    std::cout << (failures == 0
                      ? "PASS: collapse blending driver"
                      : "FAILURES PRESENT")
              << '\n';
    return failures == 0 ? 0 : 1;
}