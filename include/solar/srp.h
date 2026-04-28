#pragma once

#include "solar/force_model.h"
#include <vector>

namespace solar {

struct SRPTarget {
    int body_index;     // spacecraft index in bodies vector
    double Cr;          // reflectivity coefficient (1.0=absorber, 2.0=mirror)
    double area_m2;     // effective area facing Sun (m^2)
};

// Solar Radiation Pressure force model
class SolarRadiationPressure : public ForceModel {
public:
    // sun_index: index of Sun in bodies vector (typically 0)
    SolarRadiationPressure(int sun_index, std::vector<SRPTarget> targets,
                           bool enable_shadow = true);

    std::string name() const override { return "solar_radiation_pressure"; }

    void compute(
        const std::vector<Body>& bodies,
        double time,
        std::vector<Vec3>& acc) const override;

private:
    int sun_index_;
    std::vector<SRPTarget> targets_;
    bool enable_shadow_;

    // Cylindrical shadow: 0.0=full shadow, 1.0=full sunlight
    double shadow_factor(const Vec3& sc_pos, const Vec3& sun_pos,
                         const std::vector<Body>& bodies) const;
};

} // namespace solar
