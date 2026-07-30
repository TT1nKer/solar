#pragma once

#include <cstddef>
#include <vector>

namespace solar::relativity::detail {

struct NormalizedShadowPoint {
    long double alpha;
    long double beta;
    long double photon_radius;
};

std::vector<NormalizedShadowPoint>
sample_normalized_kerr_shadow_upper_branch(
    long double spin,
    long double sin_inclination,
    long double cos_inclination,
    long double equatorial_first_radius,
    long double equatorial_last_radius,
    std::size_t sample_count);

} // namespace solar::relativity::detail
