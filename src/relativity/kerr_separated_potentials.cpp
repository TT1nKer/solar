#include "kerr_separated_potentials.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace solar::relativity::detail {
namespace {

bool all_finite(const KerrConstants& constants) noexcept {
    return std::isfinite(constants.E) &&
           std::isfinite(constants.Lz) &&
           std::isfinite(constants.Q) &&
           std::isfinite(constants.mass_sq);
}

} // namespace

KerrSeparatedPotentials::KerrSeparatedPotentials(
    double mass_M,
    double spin_length_M,
    KerrConstants constants)
    : mass_M_(mass_M),
      spin_length_M_(spin_length_M),
      constants_(constants),
      mass_fourth_power_(mass_M * mass_M * mass_M * mass_M),
      mass_squared_(mass_M * mass_M) {
    if (!std::isfinite(mass_M_) || mass_M_ <= 0.0) {
        throw std::invalid_argument(
            "Kerr separated mass must be positive and finite");
    }
    if (!std::isfinite(spin_length_M_) ||
        !all_finite(constants_)) {
        throw std::invalid_argument(
            "Kerr separated parameters must be finite");
    }
    if (!std::isfinite(mass_fourth_power_) ||
        !std::isfinite(mass_squared_)) {
        throw std::invalid_argument(
            "Kerr separated mass scale is not representable");
    }
}

KerrSeparatedPotentialValues
KerrSeparatedPotentials::evaluate(
    double radius_M,
    double mu) const {
    if (!std::isfinite(radius_M) ||
        !std::isfinite(mu) ||
        std::fabs(mu) > 1.0) {
        throw std::domain_error(
            "Kerr separated potential point is invalid");
    }

    const double radius_sq = radius_M * radius_M;
    const double spin_sq =
        spin_length_M_ * spin_length_M_;
    const double mu_sq = mu * mu;
    const double one_minus_mu_sq = 1.0 - mu_sq;
    const double delta =
        radius_sq - 2.0 * mass_M_ * radius_M + spin_sq;
    const double sigma = radius_sq + spin_sq * mu_sq;
    const double radial_momentum =
        constants_.E * (radius_sq + spin_sq) -
        spin_length_M_ * constants_.Lz;
    const double shifted_angular_momentum =
        constants_.Lz - spin_length_M_ * constants_.E;
    const double radial_bracket =
        constants_.mass_sq * radius_sq +
        shifted_angular_momentum *
            shifted_angular_momentum +
        constants_.Q;
    const double radial =
        radial_momentum * radial_momentum -
        delta * radial_bracket;
    const double polar_coefficient =
        spin_sq *
        (constants_.mass_sq -
         constants_.E * constants_.E);
    const double polar =
        constants_.Q * one_minus_mu_sq -
        mu_sq *
            (polar_coefficient * one_minus_mu_sq +
             constants_.Lz * constants_.Lz);
    const double radial_derivative =
        4.0 * constants_.E * radius_M *
            radial_momentum -
        2.0 * (radius_M - mass_M_) *
            radial_bracket -
        delta * 2.0 * constants_.mass_sq * radius_M;
    const double polar_derivative =
        -2.0 * mu *
        (constants_.Q + polar_coefficient +
         constants_.Lz * constants_.Lz -
         2.0 * polar_coefficient * mu_sq);
    const double radial_scale =
        std::max(mass_fourth_power_, std::fabs(radial));
    const double polar_scale =
        std::max(mass_squared_, std::fabs(polar));

    const KerrSeparatedPotentialValues result{
        delta,
        sigma,
        radial,
        polar,
        radial_derivative,
        polar_derivative,
        radial_scale,
        polar_scale,
    };
    if (!std::isfinite(result.delta) ||
        !std::isfinite(result.sigma) ||
        !std::isfinite(result.radial) ||
        !std::isfinite(result.polar) ||
        !std::isfinite(result.radial_derivative) ||
        !std::isfinite(result.polar_derivative) ||
        !std::isfinite(result.radial_scale) ||
        !std::isfinite(result.polar_scale)) {
        throw std::domain_error(
            "Kerr separated potential is not finite");
    }
    return result;
}

} // namespace solar::relativity::detail
