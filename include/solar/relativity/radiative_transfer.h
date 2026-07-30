#pragma once

#include "solar/relativity/types.h"

#include <limits>
#include <string>

namespace solar::relativity {

class EmissionModel;
class FluidModel;
class Metric;

struct TransferCoefficients {
    double invariant_emissivity = 0.0;
    double invariant_absorption = 0.0;
};

enum class TransferError {
    None,
    NonFiniteInput,
    InvalidObserverFrequency,
    InvalidMetricPoint,
    InvalidFluidSample,
    FourVelocityNotUnitTimelike,
    NonFutureDirectedPhoton,
    EmissionModelFailure,
    InvalidCoefficients,
    InvalidStep,
    NonFiniteResult,
    CrossingLimitReached,
};

struct BackwardTransferState {
    double invariant_intensity = 0.0;
    double transmission = 1.0;
    double optical_depth = 0.0;
};

struct RedshiftSample {
    double normalized_emitter_frequency =
        std::numeric_limits<double>::quiet_NaN();
    double emitter_frequency =
        std::numeric_limits<double>::quiet_NaN();
    double redshift_g =
        std::numeric_limits<double>::quiet_NaN();
};

struct RedshiftResult {
    TransferError error = TransferError::None;
    RedshiftSample sample;
    std::string message;

    explicit operator bool() const noexcept {
        return error == TransferError::None;
    }
};

struct TransferAdvanceResult {
    TransferError error = TransferError::None;
    BackwardTransferState state;
    std::string message;

    explicit operator bool() const noexcept {
        return error == TransferError::None;
    }
};

struct TransferEvaluationResult {
    TransferError error = TransferError::None;
    TransferCoefficients coefficients;
    RedshiftSample redshift;
    std::string message;

    explicit operator bool() const noexcept {
        return error == TransferError::None;
    }
};

TransferAdvanceResult advance_backward_transfer(
    const BackwardTransferState& current,
    const TransferCoefficients& coefficients,
    double positive_ds);

RedshiftResult evaluate_redshift(
    const Metric& metric,
    const PhaseSpaceState& photon,
    const Contravariant4& emitter_four_velocity,
    double observer_frequency);

TransferEvaluationResult evaluate_transfer_coefficients(
    const Metric& metric,
    const PhaseSpaceState& photon,
    double observer_frequency,
    const FluidModel& fluid_model,
    const EmissionModel& emission_model);

double specific_intensity_at_observer(
    double invariant_intensity,
    double observer_frequency);

} // namespace solar::relativity
