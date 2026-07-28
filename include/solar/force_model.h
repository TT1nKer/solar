#pragma once

#include "solar/body.h"
#include <vector>
#include <string>

namespace solar {

class ForceModel {
public:
    virtual ~ForceModel() = default;

    // Name for logging/debugging
    virtual std::string name() const = 0;

    // Velocity Verlet is valid only when every active force is position-dependent.
    virtual bool depends_on_velocity() const { return false; }

    // Compute accelerations on all bodies. ADDS to acc (does not replace).
    // acc is pre-sized to bodies.size() and zero-initialized by the caller.
    virtual void compute(
        const std::vector<Body>& bodies,
        double time,
        std::vector<Vec3>& acc) const = 0;

    // Optional: potential energy contribution for energy conservation tracking.
    // Default returns 0.
    virtual double potential_energy(
        const std::vector<Body>& bodies,
        double time) const;
};

} // namespace solar
