#pragma once

#include "solar/relativity/kerr_constants.h"

namespace solar::relativity::detail {

struct KerrSeparatedPotentialValues {
    double delta;
    double sigma;
    double radial;
    double polar;
    double radial_derivative;
    double polar_derivative;
    double radial_scale;
    double polar_scale;
};

class KerrSeparatedPotentials {
public:
    KerrSeparatedPotentials(
        double mass_M,
        double spin_length_M,
        KerrConstants constants);

    KerrSeparatedPotentialValues evaluate(
        double radius_M,
        double mu) const;

    double mass() const noexcept {
        return mass_M_;
    }

private:
    double mass_M_;
    double spin_length_M_;
    KerrConstants constants_;
    double mass_fourth_power_;
    double mass_squared_;
};

} // namespace solar::relativity::detail
