#pragma once

namespace solar::relativity::detail {

struct BardeenQuantities {
    long double xi;
    long double eta;
    long double beta_radicand;
    long double radicand_scale;
};

struct VisiblePhotonInterval {
    long double first_radius;
    long double last_radius;
};

BardeenQuantities evaluate_bardeen_quantities(
    long double radius,
    long double spin,
    long double sin_inclination,
    long double cos_inclination);

long double bardeen_radicand_tolerance(
    const BardeenQuantities& quantities);

VisiblePhotonInterval find_visible_photon_interval(
    long double spin,
    long double sin_inclination,
    long double cos_inclination,
    long double equatorial_first,
    long double equatorial_last);

long double find_photon_radius_for_screen_alpha(
    long double spin,
    long double sin_inclination,
    long double first_radius,
    long double last_radius,
    long double screen_alpha);

} // namespace solar::relativity::detail
