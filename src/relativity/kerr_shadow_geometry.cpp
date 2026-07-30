#include "kerr_shadow_geometry.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace solar::relativity::detail {
namespace {

constexpr long double radicand_tolerance_factor = 128.0L;
constexpr std::size_t root_iterations = 256;

long double xi_numerator(
    long double radius,
    long double spin_squared) {
    return radius * radius * (radius - 3.0L) +
           spin_squared * (radius + 1.0L);
}

template <typename Function>
long double bisect_sign_change(
    long double first,
    long double last,
    Function function,
    bool first_is_nonnegative) {
    for (std::size_t iteration = 0;
         iteration < root_iterations;
         ++iteration) {
        const long double midpoint =
            first + 0.5L * (last - first);
        if (midpoint == first || midpoint == last) {
            break;
        }
        const bool midpoint_is_nonnegative =
            function(midpoint) >= 0.0L;
        if (midpoint_is_nonnegative == first_is_nonnegative) {
            first = midpoint;
        } else {
            last = midpoint;
        }
    }
    return first + 0.5L * (last - first);
}

} // namespace

BardeenQuantities evaluate_bardeen_quantities(
    long double radius,
    long double spin,
    long double sin_inclination,
    long double cos_inclination) {
    const long double radius_squared = radius * radius;
    const long double spin_squared = spin * spin;
    const long double radius_minus_one = radius - 1.0L;
    const long double xi =
        (radius_squared * (radius - 3.0L) +
         spin_squared * (radius + 1.0L)) /
        (-spin * radius_minus_one);
    const long double eta =
        radius * radius_squared *
        (4.0L * spin_squared -
         radius * (radius - 3.0L) *
             (radius - 3.0L)) /
        (spin_squared *
         radius_minus_one * radius_minus_one);
    const long double cotangent =
        cos_inclination / sin_inclination;
    const long double spin_projection =
        spin_squared *
        cos_inclination * cos_inclination;
    const long double angular_projection =
        xi * xi * cotangent * cotangent;
    return BardeenQuantities{
        xi,
        eta,
        eta + spin_projection - angular_projection,
        std::max(
            1.0L,
            std::fabs(eta) +
                std::fabs(spin_projection) +
                std::fabs(angular_projection)),
    };
}

long double bardeen_radicand_tolerance(
    const BardeenQuantities& quantities) {
    return radicand_tolerance_factor *
           std::numeric_limits<long double>::epsilon() *
           quantities.radicand_scale;
}

VisiblePhotonInterval find_visible_photon_interval(
    long double spin,
    long double sin_inclination,
    long double cos_inclination,
    long double equatorial_first,
    long double equatorial_last) {
    const long double input_angle_tolerance =
        64.0L * std::numeric_limits<double>::epsilon();
    if (std::fabs(cos_inclination) <=
        input_angle_tolerance) {
        return VisiblePhotonInterval{
            equatorial_first, equatorial_last};
    }

    const long double spin_squared = spin * spin;
    const long double central_radius = bisect_sign_change(
        equatorial_first,
        equatorial_last,
        [spin_squared](long double radius) {
            return xi_numerator(radius, spin_squared);
        },
        false);
    const auto radicand =
        [spin, sin_inclination, cos_inclination](
            long double radius) {
            return evaluate_bardeen_quantities(
                       radius,
                       spin,
                       sin_inclination,
                       cos_inclination)
                .beta_radicand;
        };

    const BardeenQuantities first_quantities =
        evaluate_bardeen_quantities(
            equatorial_first,
            spin,
            sin_inclination,
            cos_inclination);
    const BardeenQuantities central_quantities =
        evaluate_bardeen_quantities(
            central_radius,
            spin,
            sin_inclination,
            cos_inclination);
    const BardeenQuantities last_quantities =
        evaluate_bardeen_quantities(
            equatorial_last,
            spin,
            sin_inclination,
            cos_inclination);
    // The equatorial photon radii have eta=0 analytically. Removing the
    // rounded eta residual preserves a sign bracket arbitrarily close to an
    // equatorial viewing angle.
    const long double first_endpoint_radicand =
        first_quantities.beta_radicand -
        first_quantities.eta;
    const long double last_endpoint_radicand =
        last_quantities.beta_radicand -
        last_quantities.eta;
    if (central_quantities.beta_radicand <= 0.0L ||
        first_endpoint_radicand >= 0.0L ||
        last_endpoint_radicand >= 0.0L) {
        throw std::domain_error(
            "shadow visible interval is numerically unresolved");
    }

    const long double first_visible = bisect_sign_change(
        equatorial_first,
        central_radius,
        radicand,
        false);
    const long double last_visible = bisect_sign_change(
        central_radius,
        equatorial_last,
        radicand,
        true);
    if (!std::isfinite(first_visible) ||
        !std::isfinite(last_visible) ||
        first_visible >= last_visible) {
        throw std::domain_error(
            "shadow visible interval is not finite and ordered");
    }
    return VisiblePhotonInterval{
        first_visible, last_visible};
}

} // namespace solar::relativity::detail
