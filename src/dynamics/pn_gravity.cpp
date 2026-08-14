#include "solar/dynamics/pn_gravity.h"

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

PostNewtonianGravity::PostNewtonianGravity(double speed_of_light_km_s)
    : PostNewtonianGravity(Config{}, speed_of_light_km_s) {}

PostNewtonianGravity::PostNewtonianGravity(Config config,
                                           double speed_of_light_km_s)
    : config_(config), speed_of_light_km_s_(speed_of_light_km_s) {}

void PostNewtonianGravity::compute(
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
    const double theta = config_.opening_angle;

    for (std::size_t i = 0; i < bodies.size(); ++i) {
        const Vec3& position = bodies[i].state.pos;

        // Newtonian field and potential of the OTHER bodies, from the
        // same monopole walk as BarnesHutGravity (same self-exclusion,
        // same opening-angle policy).
        Vec3 field{};
        double potential = 0.0;

        std::vector<int> stack;
        stack.push_back(0);
        while (!stack.empty()) {
            const int node_index = stack.back();
            stack.pop_back();
            const octree_detail::Node& node =
                tree.nodes[static_cast<std::size_t>(node_index)];

            if (node.leaf && node.particle >= 0) {
                const int j = node.particle;
                if (j == static_cast<int>(i)) continue;
                const Vec3 direction =
                    bodies[static_cast<std::size_t>(j)].state.pos - position;
                const double distance_sq = direction.norm_sq() + softening_sq;
                const double inverse = 1.0 / std::sqrt(distance_sq);
                field += direction *
                         (bodies[static_cast<std::size_t>(j)].mu *
                          inverse * inverse * inverse);
                potential -= bodies[static_cast<std::size_t>(j)].mu * inverse;
                continue;
            }
            if (node.total_mass <= 0.0) continue;

            const Vec3 direction = node.com - position;
            const double distance = direction.norm();
            if (distance == 0.0) {
                for (const int child : node.children) {
                    if (child >= 0) stack.push_back(child);
                }
                continue;
            }
            const double side = 2.0 * node.half;
            const bool contains =
                std::fabs(position.x - node.center.x) <= node.half &&
                std::fabs(position.y - node.center.y) <= node.half &&
                std::fabs(position.z - node.center.z) <= node.half;
            if (!contains && side / distance < theta) {
                const double node_mu = constants::G * node.total_mass;
                const double inverse =
                    1.0 / std::sqrt(distance * distance + softening_sq);
                field += direction *
                         (node_mu * inverse * inverse * inverse);
                potential -= node_mu * inverse;
            } else {
                for (const int child : node.children) {
                    if (child >= 0) stack.push_back(child);
                }
            }
        }

        const Vec3& velocity = bodies[i].state.vel;
        const double epsilon = -potential / c_sq;  // G M_enc / (r c^2)
        double fade = 1.0;
        if (config_.fade_hi > 0.0) {
            fade = epsilon < config_.fade_hi
                ? 1.0 - smoothstep(epsilon / config_.fade_hi)
                : 0.0;
        }
        const double v_sq = velocity.norm_sq();
        const double g_dot_v = field.dot(velocity);
        acc[i] += (field * (1.0 + 4.0 * potential / c_sq - v_sq / c_sq) -
                   velocity * (4.0 * g_dot_v / c_sq)) * fade;
    }
}

} // namespace dynamics
} // namespace solar
