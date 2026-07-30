#include "kerr_shadow_sampling.h"

#include "kerr_shadow_geometry.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace solar::relativity::detail {
namespace {

long double visible_tip_alpha(
    const BardeenQuantities& quantities,
    long double spin,
    long double sin_inclination,
    long double cos_inclination) {
    if (sin_inclination >= 0.1L) {
        return -quantities.xi / sin_inclination;
    }

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
    return quantities.xi >= 0.0L
               ? -tip_magnitude
               : tip_magnitude;
}

long double stable_beta_radicand(
    const BardeenQuantities& quantities,
    long double normalized_alpha,
    long double spin,
    long double cos_inclination) {
    const long double spin_projection =
        spin * spin *
        cos_inclination * cos_inclination;
    const long double angular_projection =
        normalized_alpha * normalized_alpha *
        cos_inclination * cos_inclination;
    const long double beta_radicand =
        quantities.eta +
        spin_projection -
        angular_projection;
    BardeenQuantities stable_quantities = quantities;
    stable_quantities.beta_radicand = beta_radicand;
    stable_quantities.radicand_scale =
        std::max(
            1.0L,
            std::fabs(quantities.eta) +
                std::fabs(spin_projection) +
                std::fabs(angular_projection));
    if (beta_radicand <
        -bardeen_radicand_tolerance(
            stable_quantities)) {
        throw std::domain_error(
            "shadow sample escaped its physical visible interval");
    }
    return std::max(0.0L, beta_radicand);
}

} // namespace

std::vector<NormalizedShadowPoint>
sample_normalized_kerr_shadow_upper_branch(
    long double spin,
    long double sin_inclination,
    long double cos_inclination,
    long double equatorial_first_radius,
    long double equatorial_last_radius,
    std::size_t sample_count) {
    const VisiblePhotonInterval visible_interval =
        find_visible_photon_interval(
            spin,
            sin_inclination,
            cos_inclination,
            equatorial_first_radius,
            equatorial_last_radius);
    const BardeenQuantities first_quantities =
        evaluate_bardeen_quantities(
            visible_interval.first_radius,
            spin,
            sin_inclination,
            cos_inclination);
    const BardeenQuantities last_quantities =
        evaluate_bardeen_quantities(
            visible_interval.last_radius,
            spin,
            sin_inclination,
            cos_inclination);
    const long double first_alpha =
        visible_tip_alpha(
            first_quantities,
            spin,
            sin_inclination,
            cos_inclination);
    const long double last_alpha =
        visible_tip_alpha(
            last_quantities,
            spin,
            sin_inclination,
            cos_inclination);

    std::vector<NormalizedShadowPoint> upper_branch;
    upper_branch.reserve(sample_count);
    for (std::size_t index = 0;
         index < sample_count;
         ++index) {
        const long double fraction =
            static_cast<long double>(index) /
            static_cast<long double>(sample_count - 1);
        const long double normalized_alpha =
            first_alpha +
            fraction * (last_alpha - first_alpha);
        const bool endpoint =
            index == 0 || index + 1 == sample_count;
        const long double normalized_radius =
            index == 0
                ? visible_interval.first_radius
                : index + 1 == sample_count
                      ? visible_interval.last_radius
                      : find_photon_radius_for_screen_alpha(
                            spin,
                            sin_inclination,
                            visible_interval.first_radius,
                            visible_interval.last_radius,
                            normalized_alpha);
        long double beta_radicand = 0.0L;
        if (!endpoint) {
            const BardeenQuantities quantities =
                evaluate_bardeen_quantities(
                    normalized_radius,
                    spin,
                    sin_inclination,
                    cos_inclination);
            beta_radicand =
                stable_beta_radicand(
                    quantities,
                    normalized_alpha,
                    spin,
                    cos_inclination);
        }
        upper_branch.push_back(NormalizedShadowPoint{
            normalized_alpha,
            std::sqrt(beta_radicand),
            normalized_radius,
        });
    }
    return upper_branch;
}

} // namespace solar::relativity::detail
