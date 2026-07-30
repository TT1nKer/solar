#include "solar/relativity/kerr_shadow.h"

#include "kerr_shadow_geometry.h"
#include "solar/relativity/kerr_orbits.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace solar::relativity {
namespace {

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
    const detail::VisiblePhotonInterval visible_interval =
        detail::find_visible_photon_interval(
            spin,
            sin_inclination,
            cos_inclination,
            first_normalized_radius,
            last_normalized_radius);

    std::vector<ShadowCriticalPoint> upper_branch;
    upper_branch.reserve(samples_per_branch);
    for (std::size_t index = 0;
         index < samples_per_branch;
         ++index) {
        const long double fraction =
            static_cast<long double>(index) /
            static_cast<long double>(samples_per_branch - 1);
        const long double normalized_radius =
            visible_interval.first_radius +
            fraction *
                (visible_interval.last_radius -
                 visible_interval.first_radius);
        const long double photon_radius =
            mass * normalized_radius;
        const detail::BardeenQuantities quantities =
            detail::evaluate_bardeen_quantities(
                normalized_radius,
                spin,
                sin_inclination,
                cos_inclination);
        long double beta_radicand =
            quantities.beta_radicand;
        const bool endpoint =
            index == 0 ||
            index + 1 == samples_per_branch;
        if (endpoint) {
            beta_radicand = 0.0L;
        } else if (beta_radicand <
                   -detail::bardeen_radicand_tolerance(
                       quantities)) {
            throw std::domain_error(
                "shadow sample escaped its physical visible interval");
        } else if (beta_radicand < 0.0L) {
            beta_radicand = 0.0L;
        }

        long double normalized_alpha =
            -quantities.xi / sin_inclination;
        if (endpoint &&
            sin_inclination < 0.1L) {
            const long double tip_magnitude_squared =
                quantities.eta +
                spin * spin *
                    cos_inclination * cos_inclination;
            if (!std::isfinite(tip_magnitude_squared) ||
                tip_magnitude_squared < 0.0L) {
                throw std::domain_error(
                    "shadow tip coordinate is not physical");
            }
            const long double tip_magnitude =
                std::sqrt(tip_magnitude_squared) /
                std::fabs(cos_inclination);
            normalized_alpha =
                quantities.xi >= 0.0L
                    ? -tip_magnitude
                    : tip_magnitude;
        }
        const ShadowCriticalPoint point{
            static_cast<double>(
                mass * normalized_alpha),
            static_cast<double>(
                mass * std::sqrt(beta_radicand)),
            static_cast<double>(photon_radius)};
        if (!std::isfinite(point.alpha) ||
            !std::isfinite(point.beta) ||
            !std::isfinite(point.photon_radius)) {
            throw std::overflow_error(
                "shadow curve exceeds the finite output range");
        }
        upper_branch.push_back(point);
    }

    std::vector<ShadowCriticalPoint> curve = upper_branch;
    curve.reserve(2 * upper_branch.size() - 2);
    for (std::size_t index = upper_branch.size() - 1;
         index-- > 1;) {
        const ShadowCriticalPoint& point =
            upper_branch[index];
        curve.push_back(ShadowCriticalPoint{
            point.alpha,
            -point.beta,
            point.photon_radius,
        });
    }
    return curve;
}

} // namespace solar::relativity
