#include "solar/dynamics/barnes_hut_gravity.h"

#include "solar/constants.h"
#include "solar/dynamics/barnes_hut_octree.h"

#include <cmath>
#include <cstddef>
#include <vector>

namespace solar {
namespace dynamics {

BarnesHutGravity::BarnesHutGravity(Config config) : config_(config) {}

void BarnesHutGravity::compute(
    const std::vector<Body>& bodies,
    double /*time*/,
    std::vector<Vec3>& acc) const {
    if (bodies.empty()) return;
    const octree_detail::Tree tree(bodies);
    const double softening_sq = config_.softening_km * config_.softening_km;
    const double theta = config_.opening_angle;

    for (std::size_t i = 0; i < bodies.size(); ++i) {
        const Vec3& position = bodies[i].state.pos;

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
                Vec3 direction = bodies[static_cast<std::size_t>(j)].state.pos -
                                 position;
                const double distance_sq =
                    direction.norm_sq() + softening_sq;
                const double inverse = 1.0 / std::sqrt(distance_sq);
                acc[i] += direction *
                          (bodies[static_cast<std::size_t>(j)].mu *
                           inverse * inverse * inverse);
                continue;
            }
            if (node.total_mass <= 0.0) continue;

            Vec3 direction = node.com - position;
            const double distance = direction.norm();
            if (distance == 0.0) {
                for (const int child : node.children) {
                    if (child >= 0) stack.push_back(child);
                }
                continue;
            }
            const double side = 2.0 * node.half;
            // Never approximate a node that contains the target particle:
            // its monopole would include the particle's own mass and the
            // self-term degrades accuracy. Descend instead (standard
            // self-exclusion via geometric containment).
            const bool contains =
                std::fabs(position.x - node.center.x) <= node.half &&
                std::fabs(position.y - node.center.y) <= node.half &&
                std::fabs(position.z - node.center.z) <= node.half;
            if (!contains && side / distance < theta) {
                const double node_mu = constants::G * node.total_mass;
                const double inverse =
                    1.0 / std::sqrt(distance * distance + softening_sq);
                acc[i] += direction *
                          (node_mu * inverse * inverse * inverse);
            } else {
                for (const int child : node.children) {
                    if (child >= 0) stack.push_back(child);
                }
            }
        }
    }
}

double BarnesHutGravity::potential_energy(
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

} // namespace dynamics
} // namespace solar
