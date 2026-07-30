#include "solar/relativity/kerr_orbits.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace solar::relativity {
namespace {

void require_valid_sense(OrbitSense sense) {
    if (!is_valid_orbit_sense(sense)) {
        throw std::invalid_argument(
            "Kerr orbit sense is not recognized");
    }
}

double relative_spin_magnitude(
    const KerrBoyerLindquistMetric& metric) {
    return std::fabs(metric.spin_chi());
}

double coordinate_rotation_sign(
    const KerrBoyerLindquistMetric& metric,
    OrbitSense sense) {
    const double spin_sign =
        metric.spin_chi() < 0.0 ? -1.0 : 1.0;
    return sense == OrbitSense::Prograde
               ? spin_sign
               : -spin_sign;
}

CircularOrbitResult orbit_failure(std::string message) {
    return CircularOrbitResult{
        std::nullopt, std::move(message)};
}

} // namespace

double kerr_isco_radius(
    const KerrBoyerLindquistMetric& metric,
    OrbitSense sense) {
    require_valid_sense(sense);
    const double spin = relative_spin_magnitude(metric);
    if (spin == 0.0) {
        return 6.0 * metric.mass();
    }

    const double z1 =
        1.0 +
        std::cbrt(1.0 - spin * spin) *
            (std::cbrt(1.0 + spin) +
             std::cbrt(1.0 - spin));
    const double z2 =
        std::sqrt(3.0 * spin * spin + z1 * z1);
    const double root = std::sqrt(std::max(
        0.0,
        (3.0 - z1) *
            (3.0 + z1 + 2.0 * z2)));
    const double relative_radius =
        3.0 + z2 +
        (sense == OrbitSense::Prograde ? -root : root);
    return metric.mass() * relative_radius;
}

double kerr_equatorial_photon_radius(
    const KerrBoyerLindquistMetric& metric,
    OrbitSense sense) {
    require_valid_sense(sense);
    const double spin = relative_spin_magnitude(metric);
    if (spin == 0.0) {
        return 3.0 * metric.mass();
    }

    const double acos_argument =
        sense == OrbitSense::Prograde ? -spin : spin;
    const double relative_radius =
        2.0 *
        (1.0 +
         std::cos(
             (2.0 / 3.0) *
             std::acos(acos_argument)));
    return metric.mass() * relative_radius;
}

double kerr_marginally_bound_radius(
    const KerrBoyerLindquistMetric& metric,
    OrbitSense sense) {
    require_valid_sense(sense);
    const double spin = relative_spin_magnitude(metric);
    if (spin == 0.0) {
        return 4.0 * metric.mass();
    }

    const double relative_radius =
        sense == OrbitSense::Prograde
            ? 2.0 - spin + 2.0 * std::sqrt(1.0 - spin)
            : 2.0 + spin + 2.0 * std::sqrt(1.0 + spin);
    return metric.mass() * relative_radius;
}

CircularOrbitResult
evaluate_equatorial_circular_timelike_orbit(
    const KerrBoyerLindquistMetric& metric,
    double radius,
    OrbitSense sense) {
    require_valid_sense(sense);
    if (!std::isfinite(radius) || radius <= 0.0) {
        return orbit_failure(
            "circular orbit radius must be positive and finite");
    }
    if (radius <=
        kerr_equatorial_photon_radius(metric, sense)) {
        return orbit_failure(
            "no timelike circular orbit exists at or below "
            "the corresponding photon orbit");
    }

    const double mass = metric.mass();
    const double normalized_radius = radius / mass;
    const double root_radius = std::sqrt(normalized_radius);
    const double radius_three_halves =
        normalized_radius * root_radius;
    const double radius_three_quarters =
        std::sqrt(radius_three_halves);
    const double rotation_sign =
        coordinate_rotation_sign(metric, sense);
    const double spin = metric.spin_chi();
    const double denominator_argument =
        radius_three_halves -
        3.0 * root_radius +
        2.0 * rotation_sign * spin;
    if (!std::isfinite(denominator_argument) ||
        denominator_argument <= 0.0) {
        return orbit_failure(
            "circular timelike normalization is non-positive");
    }

    const double denominator =
        radius_three_quarters *
        std::sqrt(denominator_argument);
    const double angular_denominator =
        mass *
        (radius_three_halves + rotation_sign * spin);
    const double energy =
        (radius_three_halves -
         2.0 * root_radius +
         rotation_sign * spin) /
        denominator;
    const double angular_momentum =
        mass * rotation_sign *
        (normalized_radius * normalized_radius -
         2.0 * rotation_sign * spin * root_radius +
         spin * spin) /
        denominator;
    const double angular_velocity =
        rotation_sign / angular_denominator;
    if (!std::isfinite(energy) ||
        !std::isfinite(angular_momentum) ||
        !std::isfinite(angular_velocity)) {
        return orbit_failure(
            "circular timelike quantities are non-finite");
    }

    return CircularOrbitResult{
        CircularTimelikeOrbit{
            radius,
            angular_velocity,
            energy,
            angular_momentum,
            radius >= kerr_isco_radius(metric, sense)
                ? CircularOrbitStability::Stable
                : CircularOrbitStability::Unstable,
        },
        {},
    };
}

} // namespace solar::relativity
