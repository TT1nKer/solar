#include "solar/dynamics/pn_collapse.h"

#include "solar/constants.h"
#include "solar/dynamics/barnes_hut_octree.h"

#include <cmath>
#include <cstddef>
#include <vector>

namespace solar {
namespace dynamics {

namespace {

double smoothstep(double x) {
    return x * x * (3.0 - 2.0 * x);
}

} // namespace

PnCollapseForce::PnCollapseForce(Config config)
    : config_(config), speed_of_light_km_s_(0.0) {}

PnCollapseForce::PnCollapseForce(Config config, double speed_of_light_km_s)
    : config_(config), speed_of_light_km_s_(speed_of_light_km_s) {}

void PnCollapseForce::compute(
    const std::vector<Body>& bodies,
    double /*time*/,
    std::vector<Vec3>& acc) const {
    if (bodies.empty()) return;
    const double c = speed_of_light_km_s_ > 0.0
        ? speed_of_light_km_s_
        : constants::C_LIGHT;  // km/s
    const double c_sq = c * c;
    const octree_detail::Tree tree(bodies);
    const double softening_sq = config_.softening_km * config_.softening_km;

    for (std::size_t i = 0; i < bodies.size(); ++i) {
        Vec3 field{};
        double potential = 0.0;
        octree_detail::monopole_field(tree, bodies, i, softening_sq,
                                      config_.opening_angle, field,
                                      potential);
        const Vec3& velocity = bodies[i].state.vel;
        const double eps = -potential / c_sq;  // G M_enc / (r c^2)

        // Blending weight across the Newton -> 1PN window. eps >= eps_hi
        // wins over eps <= eps_lo so a degenerate window (eps_hi <=
        // eps_lo) selects the full 1PN force, matching the un-blended
        // PostNewtonianGravity; a very large eps_lo selects pure
        // Newtonian gravity.
        double weight = 0.0;
        if (eps >= config_.eps_hi) {
            weight = 1.0;
        } else if (eps > config_.eps_lo) {
            const double span = config_.eps_hi - config_.eps_lo;
            if (span > 0.0) {
                weight = smoothstep((eps - config_.eps_lo) / span);
            }
        }
        double fade = 1.0;
        if (config_.fade_hi > 0.0) {
            fade = eps < config_.fade_hi
                ? 1.0 - smoothstep(eps / config_.fade_hi)
                : 0.0;
        }
        const double blend = weight * fade;

        const double v_sq = velocity.norm_sq();
        const double g_dot_v = field.dot(velocity);
        // a_1PN - g = g (4 Phi / c^2 - v^2 / c^2) - 4 (g . v) v / c^2.
        const Vec3 correction =
            field * (4.0 * potential / c_sq - v_sq / c_sq) -
            velocity * (4.0 * g_dot_v / c_sq);
        acc[i] += field + correction * blend;
    }
}

double PnCollapseForce::potential_energy(
    const std::vector<Body>& bodies,
    double /*time*/) const {
    double energy = 0.0;
    const double softening_sq = config_.softening_km * config_.softening_km;
    for (std::size_t i = 0; i < bodies.size(); ++i) {
        for (std::size_t j = i + 1; j < bodies.size(); ++j) {
            const double distance_sq =
                (bodies[j].state.pos - bodies[i].state.pos).norm_sq() +
                softening_sq;
            energy -= bodies[i].mu * bodies[j].mass / std::sqrt(distance_sq);
        }
    }
    return energy;
}

std::vector<double> per_particle_compactness(
    const std::vector<Body>& bodies,
    double softening_km,
    double speed_of_light_km_s) {
    const double c = speed_of_light_km_s > 0.0
        ? speed_of_light_km_s
        : constants::C_LIGHT;
    const double c_sq = c * c;
    const double softening_sq = softening_km * softening_km;
    std::vector<double> eps(bodies.size(), 0.0);
    for (std::size_t i = 0; i < bodies.size(); ++i) {
        double potential = 0.0;
        for (std::size_t j = 0; j < bodies.size(); ++j) {
            if (i == j) continue;
            const double distance = std::sqrt(
                (bodies[j].state.pos - bodies[i].state.pos).norm_sq() +
                softening_sq);
            potential -= bodies[j].mu / distance;
        }
        eps[i] = -potential / c_sq;
    }
    return eps;
}

std::vector<int> handoff_candidates(
    const std::vector<Body>& bodies,
    double softening_km,
    double speed_of_light_km_s,
    double eps_hi) {
    const std::vector<double> eps =
        per_particle_compactness(bodies, softening_km, speed_of_light_km_s);
    std::vector<int> candidates;
    for (std::size_t i = 0; i < eps.size(); ++i) {
        if (eps[i] >= eps_hi) candidates.push_back(static_cast<int>(i));
    }
    return candidates;
}

} // namespace dynamics
} // namespace solar
