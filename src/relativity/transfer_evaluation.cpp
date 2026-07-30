#include "solar/relativity/radiative_transfer.h"

#include "solar/relativity/emission_model.h"
#include "solar/relativity/fluid_model.h"
#include "solar/relativity/metric.h"
#include "solar/relativity/spacetime_algebra.h"

#include <cmath>
#include <exception>
#include <string>

namespace solar::relativity {
namespace {

constexpr double four_velocity_tolerance = 1.0e-10;

RedshiftResult redshift_failure(
    TransferError error,
    std::string message) {
    return RedshiftResult{
        error, {}, std::move(message)};
}

TransferEvaluationResult evaluation_failure(
    TransferError error,
    std::string message) {
    return TransferEvaluationResult{
        error, {}, {}, std::move(message)};
}

TransferError validate_photon_context(
    const Metric& metric,
    const PhaseSpaceState& photon,
    double observer_frequency) {
    if (!std::isfinite(photon.affine) ||
        !photon.x.v.all_finite() ||
        !photon.p.v.all_finite()) {
        return TransferError::NonFiniteInput;
    }
    if (!std::isfinite(observer_frequency) ||
        observer_frequency <= 0.0) {
        return TransferError::InvalidObserverFrequency;
    }
    if (!metric.valid_point(photon.x)) {
        return TransferError::InvalidMetricPoint;
    }
    return TransferError::None;
}

std::string context_failure_message(TransferError error) {
    if (error == TransferError::InvalidObserverFrequency) {
        return "observer frequency must be finite and positive";
    }
    if (error == TransferError::InvalidMetricPoint) {
        return "photon position is outside the metric domain";
    }
    return "photon phase-space state must be finite";
}

bool valid_fluid_scalars(const FluidSample& fluid) {
    return std::isfinite(fluid.density) &&
           fluid.density >= 0.0 &&
           std::isfinite(fluid.temperature) &&
           fluid.temperature >= 0.0;
}

bool valid_coefficients(
    const TransferCoefficients& coefficients) {
    return
        std::isfinite(coefficients.invariant_emissivity) &&
        coefficients.invariant_emissivity >= 0.0 &&
        std::isfinite(coefficients.invariant_absorption) &&
        coefficients.invariant_absorption >= 0.0;
}

} // namespace

FluidSample VacuumFluid::sample(
    const Metric&,
    const Contravariant4&) const {
    return {};
}

RedshiftResult evaluate_redshift(
    const Metric& metric,
    const PhaseSpaceState& photon,
    const Contravariant4& emitter_four_velocity,
    double observer_frequency) {
    const TransferError context_error =
        validate_photon_context(
            metric, photon, observer_frequency);
    if (context_error != TransferError::None) {
        return redshift_failure(
            context_error,
            context_failure_message(context_error));
    }
    if (!emitter_four_velocity.v.all_finite()) {
        return redshift_failure(
            TransferError::InvalidFluidSample,
            "emitter four-velocity must be finite");
    }

    Mat4 covariant;
    try {
        covariant = metric.covariant(photon.x);
    } catch (const std::exception& error) {
        return redshift_failure(
            TransferError::InvalidMetricPoint,
            std::string("metric evaluation failed: ") +
                error.what());
    }
    const double velocity_norm =
        metric_inner_product(
            covariant,
            emitter_four_velocity,
            emitter_four_velocity);
    if (!std::isfinite(velocity_norm) ||
        std::fabs(velocity_norm + 1.0) >
            four_velocity_tolerance) {
        return redshift_failure(
            TransferError::FourVelocityNotUnitTimelike,
            "emitter four-velocity is not unit timelike");
    }

    const double normalized_emitter_frequency =
        -covector_vector_pairing(
            photon.p, emitter_four_velocity);
    if (!std::isfinite(normalized_emitter_frequency) ||
        normalized_emitter_frequency <= 0.0) {
        return redshift_failure(
            TransferError::NonFutureDirectedPhoton,
            "photon frequency is not positive for the emitter");
    }
    const double emitter_frequency =
        observer_frequency * normalized_emitter_frequency;
    const double redshift_g =
        1.0 / normalized_emitter_frequency;
    if (!std::isfinite(emitter_frequency) ||
        emitter_frequency <= 0.0 ||
        !std::isfinite(redshift_g) ||
        redshift_g <= 0.0) {
        return redshift_failure(
            TransferError::NonFiniteResult,
            "redshift exceeds the finite output range");
    }

    return RedshiftResult{
        TransferError::None,
        RedshiftSample{
            normalized_emitter_frequency,
            emitter_frequency,
            redshift_g,
        },
        {},
    };
}

TransferEvaluationResult evaluate_transfer_coefficients(
    const Metric& metric,
    const PhaseSpaceState& photon,
    double observer_frequency,
    const FluidModel& fluid_model,
    const EmissionModel& emission_model) {
    const TransferError context_error =
        validate_photon_context(
            metric, photon, observer_frequency);
    if (context_error != TransferError::None) {
        return evaluation_failure(
            context_error,
            context_failure_message(context_error));
    }

    FluidSample fluid;
    try {
        fluid = fluid_model.sample(metric, photon.x);
    } catch (const std::exception& error) {
        return evaluation_failure(
            TransferError::InvalidFluidSample,
            std::string("fluid model failed: ") +
                error.what());
    } catch (...) {
        return evaluation_failure(
            TransferError::InvalidFluidSample,
            "fluid model failed with a non-standard exception");
    }
    if (!fluid.valid) {
        return TransferEvaluationResult{};
    }
    if (!valid_fluid_scalars(fluid) ||
        !fluid.four_velocity.v.all_finite()) {
        return evaluation_failure(
            TransferError::InvalidFluidSample,
            "valid fluid sample contains invalid fields");
    }

    const RedshiftResult redshift = evaluate_redshift(
        metric,
        photon,
        fluid.four_velocity,
        observer_frequency);
    if (!redshift) {
        return evaluation_failure(
            redshift.error, redshift.message);
    }

    TransferCoefficients coefficients;
    try {
        coefficients = emission_model.coefficients(
            fluid, photon.p, observer_frequency);
    } catch (const std::exception& error) {
        return evaluation_failure(
            TransferError::EmissionModelFailure,
            std::string("emission model failed: ") +
                error.what());
    } catch (...) {
        return evaluation_failure(
            TransferError::EmissionModelFailure,
            "emission model failed with a non-standard exception");
    }
    if (!valid_coefficients(coefficients)) {
        return evaluation_failure(
            TransferError::InvalidCoefficients,
            "emission model returned invalid transfer coefficients");
    }

    return TransferEvaluationResult{
        TransferError::None,
        coefficients,
        redshift.sample,
        {},
    };
}

} // namespace solar::relativity
