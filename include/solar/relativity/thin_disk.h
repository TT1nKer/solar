#pragma once

#include "solar/relativity/fluid_model.h"
#include "solar/relativity/radiative_transfer.h"

#include <cstddef>
#include <string>
#include <vector>

namespace solar::relativity {

enum class DiskOpacityMode {
    Opaque,
    SemiTransparent,
};

struct ThinDiskSurfaceEmissionSample {
    double emitted_specific_intensity = 0.0;
    double emitted_bolometric_intensity = 0.0;
    double surface_optical_depth = 0.0;
};

class ThinDiskSurfaceEmission {
public:
    ThinDiskSurfaceEmission(
        double specific_intensity_scale,
        double bolometric_intensity_scale,
        double surface_optical_depth);

    ThinDiskSurfaceEmissionSample evaluate(
        const FluidSample& fluid) const;

private:
    friend class ThinDiskCrossingRecorder;

    double specific_intensity_scale_;
    double bolometric_intensity_scale_;
    double surface_optical_depth_;
};

struct ThinDiskCrossing {
    double affine = 0.0;
    Contravariant4 position;
    double disk_radius = 0.0;
    Contravariant4 emitter_four_velocity;
    double emitter_frequency = 0.0;
    double redshift_g = 0.0;
    Contravariant4 surface_normal;
    bool front_facing = false;
    std::size_t image_order = 0;
    double observed_temperature = 0.0;
    double observed_specific_intensity = 0.0;
    double observed_bolometric_intensity = 0.0;
};

struct ThinDiskObservedState {
    double specific_intensity = 0.0;
    double bolometric_intensity = 0.0;
    double transmission = 1.0;
};

struct ThinDiskRecorderConfig {
    DiskOpacityMode opacity_mode = DiskOpacityMode::Opaque;
    std::size_t max_crossings = 8;
};

struct ThinDiskRecordResult {
    TransferError error = TransferError::None;
    bool recorded = false;
    bool closed = false;
    std::string message;

    explicit operator bool() const noexcept {
        return error == TransferError::None;
    }
};

class ThinDiskCrossingRecorder {
public:
    ThinDiskCrossingRecorder(
        ThinDiskRecorderConfig config,
        AnalyticCircularDiskFluid disk,
        ThinDiskSurfaceEmission emission);

    ThinDiskRecordResult record(
        const Metric& metric,
        const PhaseSpaceState& photon,
        double observer_frequency);

    const std::vector<ThinDiskCrossing>& crossings() const noexcept {
        return crossings_;
    }
    const ThinDiskObservedState& observed() const noexcept {
        return observed_;
    }
    bool closed() const noexcept {
        return closed_;
    }

private:
    ThinDiskRecorderConfig config_;
    AnalyticCircularDiskFluid disk_;
    ThinDiskSurfaceEmission emission_;
    std::vector<ThinDiskCrossing> crossings_;
    ThinDiskObservedState observed_;
    bool closed_ = false;
};

} // namespace solar::relativity
