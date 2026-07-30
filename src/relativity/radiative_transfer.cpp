#include "solar/relativity/radiative_transfer.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace solar::relativity {
namespace {

TransferAdvanceResult transfer_failure(
    TransferError error,
    const BackwardTransferState& current,
    std::string message) {
    return TransferAdvanceResult{
        error, current, std::move(message)};
}

bool valid_transfer_state(
    const BackwardTransferState& state) {
    if (!std::isfinite(state.invariant_intensity) ||
        state.invariant_intensity < 0.0 ||
        !std::isfinite(state.transmission) ||
        state.transmission < 0.0 ||
        state.transmission > 1.0 ||
        std::isnan(state.optical_depth) ||
        state.optical_depth < 0.0) {
        return false;
    }
    return std::isfinite(state.optical_depth) ||
           (std::isinf(state.optical_depth) &&
            state.transmission == 0.0);
}

bool valid_coefficients(
    const TransferCoefficients& coefficients) {
    return
        std::isfinite(coefficients.invariant_emissivity) &&
        coefficients.invariant_emissivity >= 0.0 &&
        std::isfinite(coefficients.invariant_absorption) &&
        coefficients.invariant_absorption >= 0.0;
}

double attenuation_integral(
    double step_optical_depth,
    double positive_ds,
    double invariant_absorption) {
    if (step_optical_depth == 0.0) {
        return positive_ds;
    }
    if (std::isinf(step_optical_depth)) {
        return 1.0 / invariant_absorption;
    }
    if (std::fabs(step_optical_depth) < 1.0e-6) {
        const double depth_squared =
            step_optical_depth * step_optical_depth;
        return positive_ds *
               (1.0 -
                0.5 * step_optical_depth +
                depth_squared / 6.0);
    }
    return positive_ds *
           (-std::expm1(-step_optical_depth) /
            step_optical_depth);
}

} // namespace

TransferAdvanceResult advance_backward_transfer(
    const BackwardTransferState& current,
    const TransferCoefficients& coefficients,
    double positive_ds) {
    if (!valid_transfer_state(current)) {
        return transfer_failure(
            TransferError::NonFiniteInput,
            current,
            "backward transfer state is invalid");
    }
    if (!valid_coefficients(coefficients)) {
        return transfer_failure(
            TransferError::InvalidCoefficients,
            current,
            "transfer coefficients must be finite and non-negative");
    }
    if (!std::isfinite(positive_ds) || positive_ds < 0.0) {
        return transfer_failure(
            TransferError::InvalidStep,
            current,
            "backward transfer distance must be finite and non-negative");
    }
    if (positive_ds == 0.0 ||
        (coefficients.invariant_emissivity == 0.0 &&
         coefficients.invariant_absorption == 0.0)) {
        return TransferAdvanceResult{
            TransferError::None, current, {}};
    }

    const double step_optical_depth =
        coefficients.invariant_absorption * positive_ds;
    const double attenuation =
        std::isinf(step_optical_depth)
            ? 0.0
            : std::exp(-step_optical_depth);
    const double path_integral = attenuation_integral(
        step_optical_depth,
        positive_ds,
        coefficients.invariant_absorption);
    const double contribution =
        current.transmission *
        coefficients.invariant_emissivity *
        path_integral;
    const double intensity =
        current.invariant_intensity + contribution;
    const double transmission =
        current.transmission * attenuation;
    const double optical_depth =
        current.optical_depth + step_optical_depth;

    if (!std::isfinite(path_integral) ||
        !std::isfinite(contribution) ||
        !std::isfinite(intensity) ||
        !std::isfinite(transmission) ||
        std::isnan(optical_depth)) {
        return transfer_failure(
            TransferError::NonFiniteResult,
            current,
            "backward transfer step produced a non-finite result");
    }

    return TransferAdvanceResult{
        TransferError::None,
        BackwardTransferState{
            intensity,
            transmission,
            optical_depth,
        },
        {},
    };
}

double specific_intensity_at_observer(
    double invariant_intensity,
    double observer_frequency) {
    if (!std::isfinite(invariant_intensity) ||
        invariant_intensity < 0.0) {
        throw std::invalid_argument(
            "invariant intensity must be finite and non-negative");
    }
    if (!std::isfinite(observer_frequency) ||
        observer_frequency <= 0.0) {
        throw std::invalid_argument(
            "observer frequency must be finite and positive");
    }
    const double result =
        invariant_intensity *
        observer_frequency *
        observer_frequency *
        observer_frequency;
    if (!std::isfinite(result)) {
        throw std::overflow_error(
            "observer specific intensity exceeds the finite range");
    }
    return result;
}

} // namespace solar::relativity
