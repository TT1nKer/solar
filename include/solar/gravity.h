#pragma once

#include "solar/force_model.h"

namespace solar {

class NewtonianGravity : public ForceModel {
public:
    std::string name() const override { return "newtonian_gravity"; }

    // O(N^2) pairwise gravitational acceleration using mu = GM
    void compute(
        const std::vector<Body>& bodies,
        double time,
        std::vector<Vec3>& acc) const override;

    // Gravitational potential energy: sum of -G*mi*mj/rij
    double potential_energy(
        const std::vector<Body>& bodies,
        double time) const override;
};

} // namespace solar
