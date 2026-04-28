#pragma once

#include "solar/force_model.h"

namespace solar {

// Post-Newtonian Schwarzschild correction (1PN).
// Adds the leading-order general relativistic precession term
// relative to a central body (typically the Sun).
// This produces Mercury's 43 arcsec/century perihelion advance.
class GRCorrection : public ForceModel {
public:
    explicit GRCorrection(int central_body_index = 0);

    std::string name() const override { return "gr_schwarzschild"; }

    void compute(
        const std::vector<Body>& bodies,
        double time,
        std::vector<Vec3>& acc) const override;

private:
    int central_index_;
};

} // namespace solar
