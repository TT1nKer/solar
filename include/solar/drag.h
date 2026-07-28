#pragma once

#include "solar/force_model.h"
#include <vector>

namespace solar {

struct DragTarget {
    int body_index;     // spacecraft index in bodies vector
    double Cd;          // drag coefficient (typically 2.0-2.5)
    double area_m2;     // cross-sectional area (m^2)
};

// Atmospheric drag force model
// Applies drag to spacecraft in low orbit around bodies with atmosphere
class AtmosphericDrag : public ForceModel {
public:
    explicit AtmosphericDrag(std::vector<DragTarget> targets);

    std::string name() const override { return "atmospheric_drag"; }
    bool depends_on_velocity() const override { return true; }

    void compute(
        const std::vector<Body>& bodies,
        double time,
        std::vector<Vec3>& acc) const override;

private:
    std::vector<DragTarget> targets_;
};

} // namespace solar
