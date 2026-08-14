#include "solar/dynamics/barnes_hut_gravity.h"

#include "solar/constants.h"

#include <array>
#include <cmath>
#include <cstddef>

namespace solar {
namespace dynamics {
namespace {

constexpr std::size_t octree_capacity_hint = 64;

struct Node {
    Vec3 center{};      // geometric center of the cell
    double half = 0.0;  // half-side length (km)
    Vec3 com{};         // mass-weighted center of mass
    double total_mass = 0.0;  // kg
    int particle = -1;  // leaf particle index; -1 = internal or empty
    std::array<int, 8> children{};
    bool leaf = true;

    Node() { children.fill(-1); }
};

int octant_index(const Node& node, const Vec3& position) noexcept {
    int index = 0;
    if (position.x >= node.center.x) index |= 1;
    if (position.y >= node.center.y) index |= 2;
    if (position.z >= node.center.z) index |= 4;
    return index;
}

Node child_cell(const Node& parent, int octant) {
    Node child;
    child.half = 0.5 * parent.half;
    child.center = {
        parent.center.x + (octant & 1 ? child.half : -child.half),
        parent.center.y + (octant & 2 ? child.half : -child.half),
        parent.center.z + (octant & 4 ? child.half : -child.half),
    };
    return child;
}

// Insert body_index into the subtree rooted at node_index. Never holds a
// Node reference across nodes.push_back (the vector may reallocate).
void insert_particle(
    std::vector<Node>& nodes,
    int node_index,
    const std::vector<Body>& bodies,
    int body_index) {
    const Vec3& position = bodies[static_cast<std::size_t>(body_index)].state.pos;
    const double mass = bodies[static_cast<std::size_t>(body_index)].mass;

    if (nodes[static_cast<std::size_t>(node_index)].leaf &&
        nodes[static_cast<std::size_t>(node_index)].particle == -1) {
        Node& node = nodes[static_cast<std::size_t>(node_index)];
        node.particle = body_index;
        node.com = position;
        node.total_mass = mass;
        return;
    }

    if (nodes[static_cast<std::size_t>(node_index)].leaf &&
        nodes[static_cast<std::size_t>(node_index)].particle != -1) {
        const int existing =
            nodes[static_cast<std::size_t>(node_index)].particle;
        nodes[static_cast<std::size_t>(node_index)].particle = -1;
        nodes[static_cast<std::size_t>(node_index)].leaf = false;
        for (int octant = 0; octant < 8; ++octant) {
            nodes.push_back(
                child_cell(nodes[static_cast<std::size_t>(node_index)], octant));
            nodes[static_cast<std::size_t>(node_index)]
                .children[static_cast<std::size_t>(octant)] =
                static_cast<int>(nodes.size() - 1);
        }
        const int existing_octant = octant_index(
            nodes[static_cast<std::size_t>(node_index)],
            bodies[static_cast<std::size_t>(existing)].state.pos);
        const int current_octant = octant_index(
            nodes[static_cast<std::size_t>(node_index)], position);
        insert_particle(
            nodes,
            nodes[static_cast<std::size_t>(node_index)]
                .children[static_cast<std::size_t>(existing_octant)],
            bodies, existing);
        insert_particle(
            nodes,
            nodes[static_cast<std::size_t>(node_index)]
                .children[static_cast<std::size_t>(current_octant)],
            bodies, body_index);
        // Rebuild aggregates from the children.
        Node& node = nodes[static_cast<std::size_t>(node_index)];
        node.com = {};
        node.total_mass = 0.0;
        for (const int child : node.children) {
            if (child < 0) continue;
            const Node& c = nodes[static_cast<std::size_t>(child)];
            node.com += c.com * c.total_mass;
            node.total_mass += c.total_mass;
        }
        node.com = node.com * (1.0 / node.total_mass);
        return;
    }

    // Internal node: descend, then refresh aggregates.
    const int octant = octant_index(
        nodes[static_cast<std::size_t>(node_index)], position);
    insert_particle(
        nodes,
        nodes[static_cast<std::size_t>(node_index)]
            .children[static_cast<std::size_t>(octant)],
        bodies, body_index);
    Node& node = nodes[static_cast<std::size_t>(node_index)];
    const double previous_mass = node.total_mass;
    node.com = (node.com * previous_mass + position * mass) /
               (previous_mass + mass);
    node.total_mass = previous_mass + mass;
}

struct Tree {
    std::vector<Node> nodes;

    explicit Tree(const std::vector<Body>& bodies) {
        nodes.reserve(octree_capacity_hint * (bodies.size() + 1));
        if (bodies.empty()) return;
        Vec3 minimum = bodies[0].state.pos;
        Vec3 maximum = bodies[0].state.pos;
        for (const Body& body : bodies) {
            minimum.x = std::min(minimum.x, body.state.pos.x);
            minimum.y = std::min(minimum.y, body.state.pos.y);
            minimum.z = std::min(minimum.z, body.state.pos.z);
            maximum.x = std::max(maximum.x, body.state.pos.x);
            maximum.y = std::max(maximum.y, body.state.pos.y);
            maximum.z = std::max(maximum.z, body.state.pos.z);
        }
        const Vec3 center = (minimum + maximum) * 0.5;
        const double span = std::max(
            {maximum.x - minimum.x, maximum.y - minimum.y,
             maximum.z - minimum.z});
        Node root;
        root.center = center;
        root.half = 0.5 * span * (1.0 + 1.0e-9) + 1.0e-9;
        nodes.push_back(root);
        for (std::size_t index = 0; index < bodies.size(); ++index) {
            insert_particle(nodes, 0, bodies, static_cast<int>(index));
        }
    }
};

} // namespace

BarnesHutGravity::BarnesHutGravity(Config config) : config_(config) {}

void BarnesHutGravity::compute(
    const std::vector<Body>& bodies,
    double /*time*/,
    std::vector<Vec3>& acc) const {
    if (bodies.empty()) return;
    const Tree tree(bodies);
    const double softening_sq = config_.softening_km * config_.softening_km;
    const double theta = config_.opening_angle;

    for (std::size_t i = 0; i < bodies.size(); ++i) {
        const Vec3& position = bodies[i].state.pos;

        std::vector<int> stack;
        stack.push_back(0);
        while (!stack.empty()) {
            const int node_index = stack.back();
            stack.pop_back();
            const Node& node = tree.nodes[static_cast<std::size_t>(node_index)];

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
            const bool contains = std::fabs(position.x - node.center.x) <= node.half &&
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
