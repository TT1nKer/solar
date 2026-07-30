#include "solar/relativity/thin_disk.h"

#include <cmath>
#include <stdexcept>

namespace solar::relativity {
namespace {

void require_nonnegative_finite(
    double value,
    const char* message) {
    if (!std::isfinite(value) || value < 0.0) {
        throw std::invalid_argument(message);
    }
}

} // namespace

ThinDiskSurfaceEmission::ThinDiskSurfaceEmission(
    double specific_intensity_scale,
    double bolometric_intensity_scale,
    double surface_optical_depth)
    : specific_intensity_scale_(specific_intensity_scale),
      bolometric_intensity_scale_(bolometric_intensity_scale),
      surface_optical_depth_(surface_optical_depth) {
    require_nonnegative_finite(
        specific_intensity_scale,
        "surface specific intensity scale must be finite and non-negative");
    require_nonnegative_finite(
        bolometric_intensity_scale,
        "surface bolometric intensity scale must be finite and non-negative");
    if (std::isnan(surface_optical_depth) ||
        surface_optical_depth < 0.0) {
        throw std::invalid_argument(
            "surface optical depth must be non-negative");
    }
}

ThinDiskSurfaceEmissionSample
ThinDiskSurfaceEmission::evaluate(
    const FluidSample& fluid) const {
    if (!fluid.valid ||
        !std::isfinite(fluid.temperature) ||
        fluid.temperature < 0.0) {
        throw std::domain_error(
            "surface emission requires a valid finite fluid temperature");
    }
    const double temperature_squared =
        fluid.temperature * fluid.temperature;
    const ThinDiskSurfaceEmissionSample sample{
        specific_intensity_scale_ * fluid.temperature,
        bolometric_intensity_scale_ *
            temperature_squared * temperature_squared,
        surface_optical_depth_,
    };
    if (!std::isfinite(sample.emitted_specific_intensity) ||
        !std::isfinite(sample.emitted_bolometric_intensity)) {
        throw std::overflow_error(
            "surface emission exceeds the finite intensity range");
    }
    return sample;
}

} // namespace solar::relativity
