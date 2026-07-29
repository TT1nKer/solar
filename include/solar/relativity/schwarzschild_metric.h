#pragma once

#include "solar/relativity/metric.h"

namespace solar::relativity {

/**
 * Schwarzschild spacetime in Boyer-Lindquist/spherical coordinates
 * (t,r,theta,phi).
 *
 * Domain:
 *   Exterior r > 2M + margin and away from the polar coordinate axis.
 *
 * Validation:
 *   tests/relativity/test_metrics.cpp
 */
class SchwarzschildBoyerLindquistMetric final : public Metric {
public:
    explicit SchwarzschildBoyerLindquistMetric(
        double mass_M, double horizon_margin_fraction = 1.0e-8);

    Chart chart() const noexcept override {
        return Chart::BoyerLindquist;
    }
    std::string name() const override { return "schwarzschild"; }

    Mat4 covariant(const Contravariant4& x) const override;
    Mat4 contravariant(const Contravariant4& x) const override;
    std::array<Mat4, 4>
    covariant_derivatives(const Contravariant4& x) const;
    std::array<Mat4, 4>
    contravariant_derivatives(const Contravariant4& x) const override;
    bool valid_point(const Contravariant4& x) const noexcept override;

    double mass() const noexcept { return mass_M_; }
    double outer_horizon_radius() const noexcept { return 2.0 * mass_M_; }
    double horizon_margin() const noexcept { return horizon_margin_M_; }

private:
    double mass_M_;
    double horizon_margin_M_;
};

} // namespace solar::relativity
