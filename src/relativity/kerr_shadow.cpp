#include "solar/relativity/kerr_shadow.h"

#include "solar/relativity/kerr_orbits.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace solar::relativity {
namespace {

constexpr long double radicand_tolerance_factor = 128.0L;

void validate_shadow_request(
    double inclination,
    std::size_t samples_per_branch) {
    const double pi = std::acos(-1.0);
    if (!std::isfinite(inclination) ||
        inclination <= 0.0 ||
        inclination >= pi) {
        throw std::invalid_argument(
            "shadow inclination must be finite and inside (0, pi)");
    }
    if (samples_per_branch < 2) {
        throw std::invalid_argument(
            "shadow curve requires at least two samples per branch");
    }
    if (samples_per_branch >
        std::numeric_limits<std::size_t>::max() / 2 + 1) {
        throw std::invalid_argument(
            "shadow sample count is too large");
    }
}

std::vector<ShadowCriticalPoint> schwarzschild_shadow_curve(
    double mass,
    std::size_t samples_per_branch) {
    const std::size_t sample_count =
        2 * samples_per_branch - 2;
    const double radius = std::sqrt(27.0) * mass;
    const double photon_radius = 3.0 * mass;
    const double two_pi = 2.0 * std::acos(-1.0);

    std::vector<ShadowCriticalPoint> curve;
    curve.reserve(sample_count);
    for (std::size_t index = 0;
         index < sample_count;
         ++index) {
        const double angle =
            two_pi * static_cast<double>(index) /
            static_cast<double>(sample_count);
        curve.push_back(ShadowCriticalPoint{
            radius * std::cos(angle),
            radius * std::sin(angle),
            photon_radius,
        });
    }
    return curve;
}

bool is_numerically_zero_beta(
    const ShadowCriticalPoint& point,
    double mass) {
    const long double normalized_beta =
        static_cast<long double>(point.beta / mass);
    const long double tolerance =
        std::sqrt(
            radicand_tolerance_factor *
            std::numeric_limits<long double>::epsilon());
    return std::fabs(normalized_beta) <= tolerance;
}

} // namespace

std::vector<ShadowCriticalPoint> bardeen_shadow_curve(
    const KerrBoyerLindquistMetric& metric,
    double inclination,
    std::size_t samples_per_branch) {
    validate_shadow_request(inclination, samples_per_branch);

    const double spin_chi = metric.spin_chi();
    const double small_spin_limit =
        64.0 * std::sqrt(std::numeric_limits<double>::epsilon());
    if (std::fabs(spin_chi) <= small_spin_limit) {
        return schwarzschild_shadow_curve(
            metric.mass(), samples_per_branch);
    }

    const double first_radius =
        kerr_equatorial_photon_radius(
            metric, OrbitSense::Prograde);
    const double last_radius =
        kerr_equatorial_photon_radius(
            metric, OrbitSense::Retrograde);
    const long double mass = metric.mass();
    const long double first_normalized_radius =
        static_cast<long double>(first_radius) / mass;
    const long double last_normalized_radius =
        static_cast<long double>(last_radius) / mass;
    const long double spin = spin_chi;
    const long double sin_inclination =
        std::sin(static_cast<long double>(inclination));
    const long double cos_inclination =
        std::cos(static_cast<long double>(inclination));
    const long double cotangent =
        cos_inclination / sin_inclination;

    std::vector<ShadowCriticalPoint> upper_branch;
    upper_branch.reserve(samples_per_branch);
    for (std::size_t index = 0;
         index < samples_per_branch;
         ++index) {
        const long double fraction =
            static_cast<long double>(index) /
            static_cast<long double>(samples_per_branch - 1);
        const long double normalized_radius =
            first_normalized_radius +
            fraction *
                (last_normalized_radius -
                 first_normalized_radius);
        const long double photon_radius =
            mass * normalized_radius;
        const long double radius_squared =
            normalized_radius * normalized_radius;
        const long double spin_squared = spin * spin;
        const long double radius_minus_one =
            normalized_radius - 1.0L;

        const long double xi =
            (radius_squared *
                 (normalized_radius - 3.0L) +
             spin_squared *
                 (normalized_radius + 1.0L)) /
            (-spin * radius_minus_one);
        long double eta =
            normalized_radius * radius_squared *
            (4.0L * spin_squared -
             normalized_radius *
                 (normalized_radius - 3.0L) *
                 (normalized_radius - 3.0L)) /
            (spin_squared *
             radius_minus_one * radius_minus_one);
        if (index == 0 ||
            index + 1 == samples_per_branch) {
            // These radii are analytic equatorial photon orbits, so Q=0.
            eta = 0.0L;
        }

        const long double spin_projection =
            spin_squared *
            cos_inclination * cos_inclination;
        const long double angular_projection =
            xi * xi * cotangent * cotangent;
        long double beta_radicand =
            eta + spin_projection - angular_projection;
        const long double term_scale = std::max(
            1.0L,
            std::fabs(eta) +
                std::fabs(spin_projection) +
                std::fabs(angular_projection));
        const long double tolerance =
            radicand_tolerance_factor *
            std::numeric_limits<long double>::epsilon() *
            term_scale;
        if (beta_radicand < -tolerance) {
            continue;
        }
        if (beta_radicand < 0.0L) {
            beta_radicand = 0.0L;
        }

        upper_branch.push_back(ShadowCriticalPoint{
            static_cast<double>(
                -mass * xi / sin_inclination),
            static_cast<double>(
                mass * std::sqrt(beta_radicand)),
            static_cast<double>(photon_radius),
        });
    }

    std::vector<ShadowCriticalPoint> curve = upper_branch;
    curve.reserve(2 * upper_branch.size());
    for (auto point = upper_branch.rbegin();
         point != upper_branch.rend();
         ++point) {
        if (is_numerically_zero_beta(*point, metric.mass())) {
            continue;
        }
        curve.push_back(ShadowCriticalPoint{
            point->alpha,
            -point->beta,
            point->photon_radius,
        });
    }
    return curve;
}

} // namespace solar::relativity
