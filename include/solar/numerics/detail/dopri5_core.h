#ifndef SOLAR_NUMERICS_DETAIL_DOPRI5_CORE_H
#define SOLAR_NUMERICS_DETAIL_DOPRI5_CORE_H

#ifndef SOLAR_NUMERICS_DOPRI5_H
#error "include solar/numerics/dopri5.h instead"
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <utility>

namespace solar::numerics::detail {

enum class Dopri5CoreStatus {
    Completed,
    NonFiniteState,
    NonFiniteDerivative,
};

template <typename State>
struct Dopri5CoreResult {
    Dopri5CoreStatus status;
    State state;
    double step_used;
    double next_step;
    double error;
    bool accepted;
    std::array<State, 7> stages;
};

template <typename State>
bool finite_state(const State& state) {
    for (const double component : state) {
        if (!std::isfinite(component)) {
            return false;
        }
    }
    return true;
}

template <typename State>
State linear_combination(
    const State& state,
    double step,
    std::initializer_list<std::pair<double, const State*>> terms) {
    State combined = state;
    for (std::size_t i = 0; i < state.size(); ++i) {
        double increment = 0.0;
        for (const auto& term : terms) {
            increment += term.first * (*term.second)[i];
        }
        combined[i] += step * increment;
    }
    return combined;
}

template <typename State>
Dopri5CoreResult<State> failed_core_step(
    Dopri5CoreStatus status,
    const State& state,
    double step) {
    return Dopri5CoreResult<State>{
        status,
        state,
        step,
        step,
        std::numeric_limits<double>::infinity(),
        false,
        {},
    };
}

template <typename State, typename Rhs>
Dopri5CoreResult<State> dopri5_step_impl(
    const State& state,
    double independent_variable,
    double step,
    const Rhs& rhs,
    const State& absolute_tolerance,
    double relative_tolerance,
    double safety,
    double min_factor,
    double max_factor,
    ErrorNorm error_norm) {
    if (!finite_state(state)) {
        return failed_core_step(
            Dopri5CoreStatus::NonFiniteState, state, step);
    }
    if (!std::isfinite(independent_variable) ||
        !std::isfinite(step) || step == 0.0) {
        return failed_core_step(
            Dopri5CoreStatus::NonFiniteDerivative, state, step);
    }

    const State k1 = rhs(independent_variable, state);
    if (k1.size() != state.size() || !finite_state(k1)) {
        return failed_core_step(
            Dopri5CoreStatus::NonFiniteDerivative, state, step);
    }

    const State y2 = linear_combination(
        state, step, {{1.0 / 5.0, &k1}});
    if (!finite_state(y2)) {
        return failed_core_step(
            Dopri5CoreStatus::NonFiniteDerivative, state, step);
    }
    const double t2 = independent_variable + step * (1.0 / 5.0);
    if (!std::isfinite(t2)) {
        return failed_core_step(
            Dopri5CoreStatus::NonFiniteDerivative, state, step);
    }
    const State k2 = rhs(t2, y2);
    if (k2.size() != state.size() || !finite_state(k2)) {
        return failed_core_step(
            Dopri5CoreStatus::NonFiniteDerivative, state, step);
    }

    const State y3 = linear_combination(
        state, step,
        {{3.0 / 40.0, &k1}, {9.0 / 40.0, &k2}});
    if (!finite_state(y3)) {
        return failed_core_step(
            Dopri5CoreStatus::NonFiniteDerivative, state, step);
    }
    const double t3 = independent_variable + step * (3.0 / 10.0);
    if (!std::isfinite(t3)) {
        return failed_core_step(
            Dopri5CoreStatus::NonFiniteDerivative, state, step);
    }
    const State k3 = rhs(t3, y3);
    if (k3.size() != state.size() || !finite_state(k3)) {
        return failed_core_step(
            Dopri5CoreStatus::NonFiniteDerivative, state, step);
    }

    const State y4 = linear_combination(
        state, step,
        {{44.0 / 45.0, &k1},
         {-56.0 / 15.0, &k2},
         {32.0 / 9.0, &k3}});
    if (!finite_state(y4)) {
        return failed_core_step(
            Dopri5CoreStatus::NonFiniteDerivative, state, step);
    }
    const double t4 = independent_variable + step * (4.0 / 5.0);
    if (!std::isfinite(t4)) {
        return failed_core_step(
            Dopri5CoreStatus::NonFiniteDerivative, state, step);
    }
    const State k4 = rhs(t4, y4);
    if (k4.size() != state.size() || !finite_state(k4)) {
        return failed_core_step(
            Dopri5CoreStatus::NonFiniteDerivative, state, step);
    }

    const State y5 = linear_combination(
        state, step,
        {{19372.0 / 6561.0, &k1},
         {-25360.0 / 2187.0, &k2},
         {64448.0 / 6561.0, &k3},
         {-212.0 / 729.0, &k4}});
    if (!finite_state(y5)) {
        return failed_core_step(
            Dopri5CoreStatus::NonFiniteDerivative, state, step);
    }
    const double t5 = independent_variable + step * (8.0 / 9.0);
    if (!std::isfinite(t5)) {
        return failed_core_step(
            Dopri5CoreStatus::NonFiniteDerivative, state, step);
    }
    const State k5 = rhs(t5, y5);
    if (k5.size() != state.size() || !finite_state(k5)) {
        return failed_core_step(
            Dopri5CoreStatus::NonFiniteDerivative, state, step);
    }

    const State y6 = linear_combination(
        state, step,
        {{9017.0 / 3168.0, &k1},
         {-355.0 / 33.0, &k2},
         {46732.0 / 5247.0, &k3},
         {49.0 / 176.0, &k4},
         {-5103.0 / 18656.0, &k5}});
    if (!finite_state(y6)) {
        return failed_core_step(
            Dopri5CoreStatus::NonFiniteDerivative, state, step);
    }
    const double t6 = independent_variable + step;
    if (!std::isfinite(t6)) {
        return failed_core_step(
            Dopri5CoreStatus::NonFiniteDerivative, state, step);
    }
    const State k6 = rhs(t6, y6);
    if (k6.size() != state.size() || !finite_state(k6)) {
        return failed_core_step(
            Dopri5CoreStatus::NonFiniteDerivative, state, step);
    }

    const State trial = linear_combination(
        state, step,
        {{35.0 / 384.0, &k1},
         {500.0 / 1113.0, &k3},
         {125.0 / 192.0, &k4},
         {-2187.0 / 6784.0, &k5},
         {11.0 / 84.0, &k6}});
    if (!finite_state(trial)) {
        return failed_core_step(
            Dopri5CoreStatus::NonFiniteDerivative, state, step);
    }

    const State k7 = rhs(t6, trial);
    if (k7.size() != state.size() || !finite_state(k7)) {
        return failed_core_step(
            Dopri5CoreStatus::NonFiniteDerivative, state, step);
    }

    double accumulated_error = 0.0;
    double rms_scale = 0.0;
    double rms_scaled_sum_squares = 0.0;
    for (std::size_t i = 0; i < state.size(); ++i) {
        const double error_estimate =
            step * ((71.0 / 57600.0) * k1[i] -
                    (71.0 / 16695.0) * k3[i] +
                    (71.0 / 1920.0) * k4[i] -
                    (17253.0 / 339200.0) * k5[i] +
                    (22.0 / 525.0) * k6[i] -
                    (1.0 / 40.0) * k7[i]);
        const double scale =
            absolute_tolerance[i] +
            relative_tolerance *
                std::max(std::fabs(state[i]), std::fabs(trial[i]));
        const double normalized = std::fabs(error_estimate) / scale;
        if (!std::isfinite(normalized)) {
            return failed_core_step(
                Dopri5CoreStatus::NonFiniteDerivative, state, step);
        }
        if (error_norm == ErrorNorm::RootMeanSquare) {
            if (normalized != 0.0) {
                if (rms_scale < normalized) {
                    const double ratio = rms_scale / normalized;
                    rms_scaled_sum_squares =
                        1.0 +
                        rms_scaled_sum_squares * ratio * ratio;
                    rms_scale = normalized;
                } else {
                    const double ratio = normalized / rms_scale;
                    rms_scaled_sum_squares += ratio * ratio;
                }
            }
        } else {
            accumulated_error = std::max(accumulated_error, normalized);
        }
    }

    const double error =
        error_norm == ErrorNorm::RootMeanSquare
            ? rms_scale *
                  std::sqrt(
                      rms_scaled_sum_squares /
                      static_cast<double>(state.size()))
            : accumulated_error;
    const bool accepted = std::isfinite(error) && error <= 1.0;
    double factor = max_factor;
    if (error > 0.0) {
        factor = safety * std::pow(1.0 / error, 0.2);
        factor = std::clamp(factor, min_factor, max_factor);
    }
    if (!accepted) {
        factor = std::min(factor, 1.0);
    }

    return Dopri5CoreResult<State>{
        Dopri5CoreStatus::Completed,
        trial,
        step,
        std::copysign(std::fabs(step) * factor, step),
        error,
        accepted,
        {{k1, k2, k3, k4, k5, k6, k7}},
    };
}

} // namespace solar::numerics::detail

#endif // SOLAR_NUMERICS_DETAIL_DOPRI5_CORE_H
