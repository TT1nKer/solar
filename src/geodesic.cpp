#include "solar/relativity/geodesic.h"
#include "solar/integrator.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace solar::relativity {

static std::vector<double> encode_state(const GeodesicState& state) {
    return {
        state.position.t, state.position.r, state.position.theta, state.position.phi,
        state.tangent.t, state.tangent.r, state.tangent.theta, state.tangent.phi,
    };
}

static GeodesicState decode_state(const std::vector<double>& values) {
    if (values.size() != 8)
        throw std::invalid_argument("Geodesic state must contain 8 values");
    return {
        {values[0], values[1], values[2], values[3]},
        {values[4], values[5], values[6], values[7]},
    };
}

static double relative_drift(double value, double reference) {
    double scale = std::max(std::fabs(reference), 1e-30);
    return std::fabs(value - reference) / scale;
}

static GeodesicSample make_sample(
    const SchwarzschildSpacetime& spacetime,
    double affine_parameter,
    const GeodesicState& state) {
    return {
        affine_parameter,
        state,
        metric_norm(spacetime, state),
        conserved_quantities(spacetime, state),
    };
}

static void validate_options(const GeodesicOptions& options) {
    if (!std::isfinite(options.affine_duration) || options.affine_duration <= 0.0 ||
        !std::isfinite(options.initial_step) || options.initial_step <= 0.0 ||
        !std::isfinite(options.absolute_tolerance) || options.absolute_tolerance <= 0.0 ||
        !std::isfinite(options.relative_tolerance) || options.relative_tolerance < 0.0 ||
        !std::isfinite(options.horizon_margin) || options.horizon_margin <= 0.0 ||
        options.max_steps <= 0) {
        throw std::invalid_argument("Invalid geodesic propagation options");
    }
}

GeodesicResult propagate_geodesic(
    const SchwarzschildSpacetime& spacetime,
    const GeodesicState& initial_state,
    GeodesicKind kind,
    const GeodesicOptions& options) {
    validate_options(options);

    GeodesicResult propagation;
    GeodesicSample initial_sample = make_sample(spacetime, 0.0, initial_state);
    if (kind == GeodesicKind::Timelike && initial_sample.metric_norm >= 0.0)
        throw std::invalid_argument("Timelike geodesic requires a negative metric norm");
    if (kind == GeodesicKind::Null && std::fabs(initial_sample.metric_norm) > 1e-9)
        throw std::invalid_argument("Null geodesic initial metric norm must be near zero");

    propagation.samples.push_back(initial_sample);
    std::vector<double> state = encode_state(initial_state);
    double affine_parameter = 0.0;
    double step = std::min(options.initial_step, options.affine_duration);
    double horizon_guard = spacetime.horizon_radius() * (1.0 + options.horizon_margin);

    auto derivative = [&spacetime](double, const std::vector<double>& values) {
        return schwarzschild_geodesic_derivative(spacetime, values);
    };

    for (int attempt = 0; attempt < options.max_steps; ++attempt) {
        if (affine_parameter >= options.affine_duration) {
            propagation.status = PropagationStatus::Completed;
            return propagation;
        }

        GeodesicState current = decode_state(state);
        if (current.position.r <= horizon_guard) {
            propagation.status = PropagationStatus::HorizonReached;
            return propagation;
        }

        double remaining = options.affine_duration - affine_parameter;
        double proposed_step = std::min(step, remaining);

        // Keep explicit Schwarzschild-coordinate stages outside the horizon chart.
        if (current.tangent.r < 0.0) {
            double distance_to_guard = current.position.r - horizon_guard;
            double radial_limit = 0.25 * distance_to_guard / -current.tangent.r;
            proposed_step = std::min(proposed_step, radial_limit);
        }

        double minimum_step = 32.0 * std::numeric_limits<double>::epsilon() *
            std::max({1.0, spacetime.mass_length(), std::fabs(affine_parameter)});
        if (!std::isfinite(proposed_step) || proposed_step <= minimum_step) {
            propagation.status = current.tangent.r < 0.0
                ? PropagationStatus::HorizonReached
                : PropagationStatus::StepSizeUnderflow;
            return propagation;
        }

        GenericAdaptiveResult integration;
        try {
            integration = dopri5_generic_step(
                state, affine_parameter, proposed_step, derivative,
                options.absolute_tolerance, options.relative_tolerance);
        } catch (const std::domain_error&) {
            step = proposed_step * 0.2;
            propagation.rejected_steps++;
            continue;
        }

        step = integration.dt_next;
        if (!integration.accepted) {
            propagation.rejected_steps++;
            continue;
        }

        state = std::move(integration.y);
        affine_parameter += integration.dt_used;
        propagation.accepted_steps++;

        GeodesicState accepted_state = decode_state(state);
        if (accepted_state.position.r <= horizon_guard) {
            propagation.status = PropagationStatus::HorizonReached;
            return propagation;
        }

        GeodesicSample sample = make_sample(spacetime, affine_parameter, accepted_state);
        propagation.max_metric_norm_drift = std::max(
            propagation.max_metric_norm_drift,
            std::fabs(sample.metric_norm - initial_sample.metric_norm));
        propagation.max_relative_energy_drift = std::max(
            propagation.max_relative_energy_drift,
            relative_drift(sample.conserved.energy, initial_sample.conserved.energy));
        propagation.max_relative_angular_momentum_drift = std::max(
            propagation.max_relative_angular_momentum_drift,
            relative_drift(
                sample.conserved.angular_momentum_z,
                initial_sample.conserved.angular_momentum_z));
        propagation.samples.push_back(sample);
    }

    propagation.status = PropagationStatus::StepLimitReached;
    return propagation;
}

} // namespace solar::relativity
