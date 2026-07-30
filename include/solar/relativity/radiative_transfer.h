#pragma once

#include <string>

namespace solar::relativity {

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

struct TransferAdvanceResult {
    TransferError error = TransferError::None;
    BackwardTransferState state;
    std::string message;

    explicit operator bool() const noexcept {
        return error == TransferError::None;
    }
};

TransferAdvanceResult advance_backward_transfer(
    const BackwardTransferState& current,
    const TransferCoefficients& coefficients,
    double positive_ds);

double specific_intensity_at_observer(
    double invariant_intensity,
    double observer_frequency);

} // namespace solar::relativity
