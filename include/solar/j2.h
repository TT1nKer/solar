#pragma once

#include "solar/force_model.h"
#include <vector>

namespace solar {

// Specification of a body with J2 oblateness
struct J2Entry {
    int body_index;   // index of the oblate body in the bodies vector
    double J2;        // J2 zonal harmonic coefficient
    double R_eq;      // equatorial radius (km)
};

// J2 oblateness perturbation.
// Adds the acceleration due to the equatorial bulge of oblate bodies.
// Affects all other bodies orbiting within the sphere of influence.
class J2Perturbation : public ForceModel {
public:
    explicit J2Perturbation(std::vector<J2Entry> entries);

    // Convenience: create with standard solar system J2 values
    // (requires bodies vector to look up indices by name)
    static J2Perturbation solar_system_defaults(const std::vector<Body>& bodies);

    std::string name() const override { return "j2_oblateness"; }

    void compute(
        const std::vector<Body>& bodies,
        double time,
        std::vector<Vec3>& acc) const override;

private:
    std::vector<J2Entry> entries_;
};

} // namespace solar
