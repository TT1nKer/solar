#include "solar/numerics/dopri5.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using solar::numerics::Dopri5Config;
using solar::numerics::Dopri5StepResult;
using solar::numerics::ErrorNorm;
using solar::numerics::StateN;
using StepStatus = Dopri5StepResult<1>::Status;

namespace {

int passed = 0;
int failed = 0;

void check(const char* name, bool condition) {
    std::cout << (condition ? "  PASS: " : "  FAIL: ") << name << '\n';
    condition ? ++passed : ++failed;
}

void check_near(const char* name, double actual, double expected,
                double tolerance) {
    check(name, std::isfinite(actual) &&
                    std::fabs(actual - expected) <= tolerance);
}

template <typename Action>
void check_invalid_argument(const char* name, Action action) {
    try {
        action();
        check(name, false);
    } catch (const std::invalid_argument&) {
        check(name, true);
    }
}

template <typename Action>
void check_out_of_range(const char* name, Action action) {
    try {
        action();
        check(name, false);
    } catch (const std::out_of_range&) {
        check(name, true);
    }
}

Dopri5Config<1> reference_config() {
    return Dopri5Config<1>{
        StateN<1>{{1.0e-12}},
        1.0e-11,
        0.9,
        0.2,
        5.0,
        ErrorNorm::RootMeanSquare,
    };
}

} // namespace

int main() {
    const auto exponential_rhs =
        [](double, const StateN<1>& state) { return state; };

    const StateN<1> initial{{1.0}};
    const Dopri5Config<1> config = reference_config();
    const auto step = solar::numerics::dopri5_step(
        initial, 0.0, 0.01, exponential_rhs, config);

    check("fixed DOPRI step completed",
          step.status == StepStatus::Completed);
    check("fixed DOPRI step accepted", step.accepted);
    check_near("fifth-order exponential step",
               step.state[0], std::exp(0.01), 3.0e-14);
    check("accepted next step preserves direction", step.next_step > 0.0);
    check_near("reported step is the requested step",
               step.step_used, 0.01, 0.0);

    Dopri5Config<1> strict_config = config;
    strict_config.absolute_tolerance[0] = 1.0e-14;
    strict_config.relative_tolerance = 1.0e-14;
    const auto rejected = solar::numerics::dopri5_step(
        initial, 0.0, 1.0, exponential_rhs, strict_config);
    check("large strict step completed",
          rejected.status == StepStatus::Completed);
    check("large strict step rejected", !rejected.accepted);
    check("rejected next step shrinks",
          std::fabs(rejected.next_step) <
              std::fabs(rejected.step_used));

    const auto backward = solar::numerics::dopri5_step(
        initial, 0.0, -0.01, exponential_rhs, config);
    check("backward step completed",
          backward.status == StepStatus::Completed);
    check("backward next step preserves direction",
          backward.next_step < 0.0);
    check_near("backward exponential state",
               backward.state[0], std::exp(-0.01), 3.0e-14);

    check_invalid_argument("zero step rejected", [&] {
        (void)solar::numerics::dopri5_step(
            initial, 0.0, 0.0, exponential_rhs, config);
    });
    check_invalid_argument("non-finite step rejected", [&] {
        (void)solar::numerics::dopri5_step(
            initial, 0.0, std::numeric_limits<double>::infinity(),
            exponential_rhs, config);
    });
    check_invalid_argument("non-finite independent variable rejected", [&] {
        (void)solar::numerics::dopri5_step(
            initial, std::numeric_limits<double>::quiet_NaN(),
            0.1, exponential_rhs, config);
    });

    Dopri5Config<1> invalid_config = config;
    invalid_config.absolute_tolerance[0] = 0.0;
    check_invalid_argument("non-positive absolute tolerance rejected", [&] {
        (void)solar::numerics::dopri5_step(
            initial, 0.0, 0.1, exponential_rhs, invalid_config);
    });

    invalid_config = config;
    invalid_config.relative_tolerance = 0.0;
    check_invalid_argument("non-positive relative tolerance rejected", [&] {
        (void)solar::numerics::dopri5_step(
            initial, 0.0, 0.1, exponential_rhs, invalid_config);
    });

    invalid_config = config;
    invalid_config.safety = 1.0;
    check_invalid_argument("unsafe controller safety rejected", [&] {
        (void)solar::numerics::dopri5_step(
            initial, 0.0, 0.1, exponential_rhs, invalid_config);
    });

    invalid_config = config;
    invalid_config.min_factor = 1.01;
    check_invalid_argument("minimum factor above one rejected", [&] {
        (void)solar::numerics::dopri5_step(
            initial, 0.0, 0.1, exponential_rhs, invalid_config);
    });

    invalid_config = config;
    invalid_config.max_factor = 0.99;
    check_invalid_argument("maximum factor below one rejected", [&] {
        (void)solar::numerics::dopri5_step(
            initial, 0.0, 0.1, exponential_rhs, invalid_config);
    });

    StateN<1> non_finite_state{
        {std::numeric_limits<double>::quiet_NaN()}};
    const auto invalid_state = solar::numerics::dopri5_step(
        non_finite_state, 0.0, 0.1, exponential_rhs, config);
    check("non-finite state has explicit status",
          invalid_state.status == StepStatus::NonFiniteState);
    check("non-finite state is not accepted", !invalid_state.accepted);

    const auto non_finite_rhs =
        [](double, const StateN<1>&) {
            return StateN<1>{{
                std::numeric_limits<double>::infinity()}};
        };
    const auto invalid_derivative = solar::numerics::dopri5_step(
        initial, 0.0, 0.1, non_finite_rhs, config);
    check("non-finite derivative has explicit status",
          invalid_derivative.status ==
              StepStatus::NonFiniteDerivative);
    check("non-finite derivative is not accepted",
          !invalid_derivative.accepted);

    bool received_non_finite_time = false;
    const auto time_guarded_rhs =
        [&received_non_finite_time](double time, const StateN<1>&) {
            received_non_finite_time =
                received_non_finite_time || !std::isfinite(time);
            return StateN<1>{{0.0}};
        };
    const auto overflowed_stage_time =
        solar::numerics::dopri5_step(
            initial,
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            time_guarded_rhs,
            config);
    check("overflowed stage time has explicit status",
          overflowed_stage_time.status ==
              StepStatus::NonFiniteDerivative);
    check("RHS never receives non-finite stage time",
          !received_non_finite_time);

    Dopri5Config<1> dense_config = config;
    dense_config.absolute_tolerance[0] = 1.0e-9;
    dense_config.relative_tolerance = 1.0e-6;
    const auto dense_step = solar::numerics::dopri5_step(
        initial, 0.0, 0.2, exponential_rhs, dense_config);
    check("dense-output step accepted", dense_step.accepted);
    check("dense output is available",
          dense_step.dense_output.has_value());
    const auto dense = dense_step.dense_output.value();
    check_near("dense start", dense.evaluate(0.0)[0],
               1.0, 0.0);
    check_near("dense midpoint", dense.evaluate(0.1)[0],
               std::exp(0.1), 2.0e-7);
    check_near("dense end", dense.evaluate(0.2)[0],
               dense_step.state[0], 2.0e-15);
    check_out_of_range("dense rejects value before interval", [&] {
        (void)dense.evaluate(-1.0e-6);
    });
    check_out_of_range("dense rejects value after interval", [&] {
        (void)dense.evaluate(0.200001);
    });

    const auto backward_dense_step =
        solar::numerics::dopri5_step(
            initial, 0.0, -0.2, exponential_rhs, dense_config);
    check("backward dense-output step accepted",
          backward_dense_step.accepted);
    check("backward dense output is available",
          backward_dense_step.dense_output.has_value());
    const auto backward_dense =
        backward_dense_step.dense_output.value();
    check_near("backward dense start",
               backward_dense.evaluate(0.0)[0], 1.0, 0.0);
    check_near("backward dense midpoint",
               backward_dense.evaluate(-0.1)[0],
               std::exp(-0.1), 2.0e-7);
    check_near("backward dense end",
               backward_dense.evaluate(-0.2)[0],
               backward_dense_step.state[0], 2.0e-15);
    check_out_of_range(
        "backward dense rejects value above interval", [&] {
            (void)backward_dense.evaluate(1.0e-6);
        });
    check_out_of_range(
        "backward dense rejects value below interval", [&] {
            (void)backward_dense.evaluate(-0.200001);
        });

    std::cout << "\n=== Results: " << passed << " passed, "
              << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
