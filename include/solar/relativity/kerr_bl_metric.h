#pragma once

#include "solar/relativity/metric.h"

namespace solar::relativity {

/**
 * Kerr spacetime in Boyer-Lindquist coordinates (t,r,theta,phi).
 *
 * Parameters:
 *   mass_M is the geometrized mass length. spin_chi=a/M is dimensionless and
 *   must satisfy abs(spin_chi)<1.
 *
 * Domain:
 *   Regular exterior r>r_++margin and away from the polar coordinate axis.
 *
 * Validation:
 *   tests/relativity/test_kerr_bl.cpp
 *   tests/relativity/test_metric_derivatives.cpp
 */
class KerrBoyerLindquistMetric final : public Metric {
public:
    KerrBoyerLindquistMetric(
        double mass_M,
        double spin_chi,
        double horizon_margin_fraction = 1.0e-8);

    Chart chart() const noexcept override {
        return Chart::BoyerLindquist;
    }
    std::string name() const override { return "kerr-bl"; }

    Mat4 covariant(const Contravariant4& x) const override;
    Mat4 contravariant(const Contravariant4& x) const override;
    std::array<Mat4, 4>
    covariant_derivatives(const Contravariant4& x) const;
    std::array<Mat4, 4>
    contravariant_derivatives(const Contravariant4& x) const override;
    bool valid_point(const Contravariant4& x) const noexcept override;

    double mass() const noexcept { return mass_M_; }
    double spin_chi() const noexcept { return spin_chi_; }
    double spin_length() const noexcept { return spin_a_M_; }
    double horizon_margin() const noexcept { return horizon_margin_M_; }
    double outer_horizon_radius() const noexcept { return outer_horizon_M_; }
    double inner_horizon_radius() const noexcept { return inner_horizon_M_; }
    double outer_stationary_limit_radius(double theta) const;

private:
    double mass_M_;
    double spin_chi_;
    double spin_a_M_;
    double horizon_margin_M_;
    double outer_horizon_M_;
    double inner_horizon_M_;
};

} // namespace solar::relativity
