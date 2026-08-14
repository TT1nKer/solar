#pragma once

#include "solar/force_model.h"

#include <vector>

namespace solar {
namespace dynamics {

// Barnes-Hut octree gravity with monopole (center-of-mass) approximations
// and Plummer softening. Multi-mass by construction: node monopoles carry
// the mass-weighted center of mass, and leaf pairs use the exact pairwise
// mu-based acceleration.
//
// Accuracy is controlled by the opening angle theta: a node is accepted as
// a monopole when (node side / distance) < theta. theta -> 0 reduces to
// direct summation at O(N^2) cost; typical research values are 0.5-0.8.
class BarnesHutGravity : public ForceModel {
public:
    struct Config {
        double opening_angle = 0.7;   // theta, in (0, 1]
        double softening_km = 0.0;    // Plummer softening length (km)
    };

    BarnesHutGravity() = default;
    explicit BarnesHutGravity(Config config);

    std::string name() const override { return "barnes_hut"; }

    // Adds monopole/leaf accelerations to acc (km/s^2).
    void compute(
        const std::vector<Body>& bodies,
        double time,
        std::vector<Vec3>& acc) const override;

    // Exact direct-sum potential energy WITH the same Plummer softening,
    // so energy-conservation diagnostics are not polluted by the tree
    // approximation. O(N^2); intended for tests and monitoring.
    double potential_energy(
        const std::vector<Body>& bodies,
        double time) const override;

private:
    Config config_;
};

} // namespace dynamics
} // namespace solar
