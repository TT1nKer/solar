#pragma once

#include "solar/relativity/metric.h"

namespace solar::relativity {

/**
 * Minkowski spacetime in Cartesian coordinates (t,x,y,z).
 *
 * Validation:
 *   tests/relativity/test_metrics.cpp
 */
class MinkowskiMetric final : public Metric {
public:
    Chart chart() const noexcept override {
        return Chart::MinkowskiCartesian;
    }
    std::string name() const override { return "minkowski"; }

    Mat4 covariant(const Contravariant4& x) const override;
    Mat4 contravariant(const Contravariant4& x) const override;
    std::array<Mat4, 4>
    contravariant_derivatives(const Contravariant4& x) const override;
    bool valid_point(const Contravariant4& x) const noexcept override;
};

} // namespace solar::relativity
