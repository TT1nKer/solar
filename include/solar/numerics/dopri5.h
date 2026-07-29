#ifndef SOLAR_NUMERICS_DOPRI5_H
#define SOLAR_NUMERICS_DOPRI5_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <utility>

namespace solar::numerics {

template <std::size_t N>
using StateN = std::array<double, N>;

enum class ErrorNorm {
    RootMeanSquare,
    Maximum,
};

template <std::size_t N>
struct Dopri5Config {
    StateN<N> absolute_tolerance;
    double relative_tolerance;
    double safety;
    double min_factor;
    double max_factor;
    ErrorNorm error_norm;
};

template <std::size_t N>
struct Dopri5StepResult;

template <std::size_t N, typename Rhs>
Dopri5StepResult<N> dopri5_step(
    const StateN<N>& state,
    double independent_variable,
    double step,
    const Rhs& rhs,
    const Dopri5Config<N>& config);

template <std::size_t N>
class Dopri5DenseOutput {
public:
    double start() const noexcept {
        return start_;
    }

    double end() const noexcept {
        return start_ + step_;
    }

    StateN<N> evaluate(double independent_variable) const {
        const double interval_end = end();
        const double lower = std::min(start_, interval_end);
        const double upper = std::max(start_, interval_end);
        if (!std::isfinite(independent_variable) ||
            independent_variable < lower ||
            independent_variable > upper) {
            throw std::out_of_range(
                "DOPRI5 dense output cannot extrapolate");
        }

        const double theta =
            (independent_variable - start_) / step_;
        StateN<N> interpolated = initial_state_;
        for (std::size_t i = 0; i < N; ++i) {
            interpolated[i] +=
                step_ * theta *
                (coefficients_[0][i] +
                 theta * (coefficients_[1][i] +
                          theta * (coefficients_[2][i] +
                                   theta * coefficients_[3][i])));
        }
        return interpolated;
    }

private:
    Dopri5DenseOutput(
        double start,
        double step,
        const StateN<N>& initial_state,
        const std::array<StateN<N>, 4>& coefficients)
        : start_(start),
          step_(step),
          initial_state_(initial_state),
          coefficients_(coefficients) {}

    template <std::size_t M, typename Rhs>
    friend Dopri5StepResult<M> dopri5_step(
        const StateN<M>& state,
        double independent_variable,
        double step,
        const Rhs& rhs,
        const Dopri5Config<M>& config);

    double start_;
    double step_;
    StateN<N> initial_state_;
    std::array<StateN<N>, 4> coefficients_;
};

template <std::size_t N>
struct Dopri5StepResult {
    enum class Status {
        Completed,
        NonFiniteState,
        NonFiniteDerivative,
    };

    Status status;
    StateN<N> state;
    double step_used;
    double next_step;
    double error;
    bool accepted;
    std::optional<Dopri5DenseOutput<N>> dense_output;
};

} // namespace solar::numerics

#include "solar/numerics/detail/dopri5_core.h"
#include "solar/numerics/detail/dopri5_dense.h"

namespace solar::numerics {

template <std::size_t N, typename Rhs>
Dopri5StepResult<N> dopri5_step(
    const StateN<N>& state,
    double independent_variable,
    double step,
    const Rhs& rhs,
    const Dopri5Config<N>& config) {
    if constexpr (N == 0) {
        throw std::invalid_argument("DOPRI5 state must not be empty");
    }
    if (!std::isfinite(independent_variable)) {
        throw std::invalid_argument(
            "DOPRI5 independent variable must be finite");
    }
    if (!std::isfinite(step) || step == 0.0) {
        throw std::invalid_argument(
            "DOPRI5 step must be finite and non-zero");
    }
    if (!std::isfinite(config.relative_tolerance) ||
        config.relative_tolerance <= 0.0) {
        throw std::invalid_argument(
            "DOPRI5 relative tolerance must be positive and finite");
    }
    for (const double tolerance : config.absolute_tolerance) {
        if (!std::isfinite(tolerance) || tolerance <= 0.0) {
            throw std::invalid_argument(
                "DOPRI5 absolute tolerances must be positive and finite");
        }
    }
    if (!std::isfinite(config.safety) ||
        config.safety <= 0.0 || config.safety >= 1.0) {
        throw std::invalid_argument(
            "DOPRI5 safety must be finite and between zero and one");
    }
    if (!std::isfinite(config.min_factor) ||
        !std::isfinite(config.max_factor) ||
        config.min_factor <= 0.0 ||
        config.min_factor > 1.0 ||
        config.max_factor < 1.0 ||
        config.min_factor > config.max_factor) {
        throw std::invalid_argument(
            "DOPRI5 factors must satisfy 0 < min <= 1 <= max");
    }
    if (config.error_norm != ErrorNorm::RootMeanSquare &&
        config.error_norm != ErrorNorm::Maximum) {
        throw std::invalid_argument(
            "DOPRI5 error norm is not recognized");
    }

    const auto core = detail::dopri5_step_impl(
        state,
        independent_variable,
        step,
        rhs,
        config.absolute_tolerance,
        config.relative_tolerance,
        config.safety,
        config.min_factor,
        config.max_factor,
        config.error_norm);

    using Status = typename Dopri5StepResult<N>::Status;
    Status status = Status::Completed;
    if (core.status == detail::Dopri5CoreStatus::NonFiniteState) {
        status = Status::NonFiniteState;
    } else if (
        core.status == detail::Dopri5CoreStatus::NonFiniteDerivative) {
        status = Status::NonFiniteDerivative;
    }

    std::optional<Dopri5DenseOutput<N>> dense_output;
    if (status == Status::Completed) {
        const auto coefficients =
            detail::dopri5_dense_coefficients<N>(core.stages);
        dense_output = Dopri5DenseOutput<N>(
            independent_variable, step, state, coefficients);
    }

    return Dopri5StepResult<N>{
        status,
        core.state,
        core.step_used,
        core.next_step,
        core.error,
        core.accepted,
        std::move(dense_output),
    };
}

} // namespace solar::numerics

#endif // SOLAR_NUMERICS_DOPRI5_H
