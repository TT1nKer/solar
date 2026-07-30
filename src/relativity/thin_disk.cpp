#include "solar/relativity/thin_disk.h"

#include "thin_disk_geometry.h"
#include "solar/relativity/spacetime_algebra.h"

#include <cmath>
#include <exception>
#include <stdexcept>
#include <string>
#include <utility>

namespace solar::relativity {
namespace {

constexpr std::size_t maximum_crossing_count = 1024;

ThinDiskRecordResult record_failure(
    TransferError error,
    bool closed,
    std::string message) {
    return ThinDiskRecordResult{
        error, false, closed, std::move(message)};
}

bool valid_opacity_mode(DiskOpacityMode mode) noexcept {
    return mode == DiskOpacityMode::Opaque ||
           mode == DiskOpacityMode::SemiTransparent;
}

} // namespace

ThinDiskCrossingRecorder::ThinDiskCrossingRecorder(
    ThinDiskRecorderConfig config,
    AnalyticCircularDiskFluid disk,
    ThinDiskSurfaceEmission emission)
    : config_(config),
      disk_(std::move(disk)),
      emission_(std::move(emission)) {
    if (!valid_opacity_mode(config_.opacity_mode)) {
        throw std::invalid_argument(
            "disk opacity mode is not recognized");
    }
    if (config_.max_crossings == 0 ||
        config_.max_crossings > maximum_crossing_count) {
        throw std::invalid_argument(
            "disk crossing bound must be between one and 1024");
    }
    if (config_.opacity_mode ==
            DiskOpacityMode::SemiTransparent &&
        std::isinf(emission_.surface_optical_depth_)) {
        throw std::invalid_argument(
            "semi-transparent surface depth must be finite");
    }
    crossings_.reserve(config_.max_crossings);
}

ThinDiskRecordResult ThinDiskCrossingRecorder::record(
    const Metric& metric,
    const PhaseSpaceState& photon,
    double observer_frequency) {
    if (closed_) {
        return ThinDiskRecordResult{
            TransferError::None, false, true, {}};
    }

    FluidSample fluid;
    try {
        fluid = disk_.sample(metric, photon.x);
    } catch (const std::exception& error) {
        return record_failure(
            TransferError::InvalidFluidSample,
            closed_,
            std::string("disk fluid failed: ") + error.what());
    }
    if (!fluid.valid) {
        return ThinDiskRecordResult{
            TransferError::None, false, false, {}};
    }
    if (crossings_.size() >= config_.max_crossings) {
        return record_failure(
            TransferError::CrossingLimitReached,
            closed_,
            "thin-disk crossing bound reached");
    }

    const RedshiftResult redshift = evaluate_redshift(
        metric,
        photon,
        fluid.four_velocity,
        observer_frequency);
    if (!redshift) {
        return record_failure(
            redshift.error, closed_, redshift.message);
    }

    detail::ThinDiskSurfaceGeometry geometry;
    try {
        geometry =
            detail::evaluate_thin_disk_surface_geometry(
                metric, photon.x);
        if (!detail::valid_thin_disk_surface_geometry(
                metric,
                photon.x,
                geometry.normal,
                fluid.four_velocity)) {
            return record_failure(
                TransferError::InvalidFluidSample,
                closed_,
                "thin-disk normal is not unit and orthogonal");
        }
    } catch (const std::exception& error) {
        return record_failure(
            TransferError::InvalidMetricPoint,
            closed_,
            std::string("thin-disk geometry failed: ") +
                error.what());
    }

    ThinDiskSurfaceEmissionSample emitted;
    try {
        emitted = emission_.evaluate(fluid);
    } catch (const std::overflow_error& error) {
        return record_failure(
            TransferError::NonFiniteResult,
            closed_,
            error.what());
    } catch (const std::exception& error) {
        return record_failure(
            TransferError::EmissionModelFailure,
            closed_,
            std::string("surface emission failed: ") +
                error.what());
    }

    const double g = redshift.sample.redshift_g;
    const double g_squared = g * g;
    const double observed_specific =
        g_squared * g * emitted.emitted_specific_intensity;
    const double observed_bolometric =
        g_squared * g_squared *
        emitted.emitted_bolometric_intensity;
    const double observed_temperature =
        g * fluid.temperature;
    if (!std::isfinite(observed_temperature) ||
        !std::isfinite(observed_specific) ||
        !std::isfinite(observed_bolometric)) {
        return record_failure(
            TransferError::NonFiniteResult,
            closed_,
            "thin-disk redshifted source is non-finite");
    }

    const ThinDiskCrossing crossing{
        photon.affine,
        photon.x,
        geometry.radius,
        fluid.four_velocity,
        redshift.sample.emitter_frequency,
        g,
        geometry.normal,
        covector_vector_pairing(
            photon.p, geometry.normal) > 0.0,
        crossings_.size(),
        observed_temperature,
        observed_specific,
        observed_bolometric,
    };

    ThinDiskObservedState next = observed_;
    bool next_closed = false;
    if (config_.opacity_mode == DiskOpacityMode::Opaque) {
        next.specific_intensity +=
            observed_.transmission * observed_specific;
        next.bolometric_intensity +=
            observed_.transmission * observed_bolometric;
        next.transmission = 0.0;
        next_closed = true;
    } else {
        const double sheet_absorption =
            -std::expm1(-emitted.surface_optical_depth);
        const double sheet_transmission =
            std::exp(-emitted.surface_optical_depth);
        next.specific_intensity +=
            observed_.transmission *
            observed_specific *
            sheet_absorption;
        next.bolometric_intensity +=
            observed_.transmission *
            observed_bolometric *
            sheet_absorption;
        next.transmission *= sheet_transmission;
    }
    if (!std::isfinite(next.specific_intensity) ||
        !std::isfinite(next.bolometric_intensity) ||
        !std::isfinite(next.transmission)) {
        return record_failure(
            TransferError::NonFiniteResult,
            closed_,
            "thin-disk composition is non-finite");
    }

    crossings_.push_back(crossing);
    observed_ = next;
    closed_ = next_closed;

    return ThinDiskRecordResult{
        TransferError::None, true, closed_, {}};
}

} // namespace solar::relativity
