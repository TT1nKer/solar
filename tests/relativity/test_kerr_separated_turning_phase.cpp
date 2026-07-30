#include "../../src/relativity/kerr_separated_turning.h"

#include "solar/numerics/dopri5.h"

#include <cmath>
#include <iostream>
#include <string>

using namespace solar::relativity;
namespace numerics = solar::numerics;

namespace {

int passed = 0;
int failed = 0;

void check(const std::string& name, bool condition) {
    if (condition) {
        std::cout << "  PASS: " << name << "\n";
        ++passed;
    } else {
        std::cerr << "  FAIL: " << name << "\n";
        ++failed;
    }
}

numerics::Dopri5Config<7> test_config() {
    return numerics::Dopri5Config<7>{
        numerics::StateN<7>{{
            1.0e-14,
            1.0e-14,
            1.0e-14,
            1.0e-14,
            1.0e-14,
            1.0e-14,
            1.0e-14,
        }},
        1.0e-13,
        0.9,
        0.2,
        5.0,
        numerics::ErrorNorm::RootMeanSquare,
    };
}

} // namespace

int main() {
    const double carter = 4.0 / 3.0 + 0.01;
    const double polar_frequency_sq = carter + 4.0;
    const detail::KerrSeparatedPotentials potentials(
        1.0,
        0.0,
        KerrConstants{1.0, 2.0, carter, 0.0});
    const auto initial_values =
        potentials.evaluate(10.0, 0.5);
    const detail::KerrTurningPhaseState initial{{
        0.0,
        10.0,
        0.0,
        0.5,
        std::sqrt(initial_values.polar),
        0.0,
        0.0,
    }};
    const auto polar_rhs =
        [polar_frequency_sq](
            double,
            const detail::KerrTurningPhaseState& state) {
            detail::KerrTurningPhaseState derivative{};
            derivative[detail::kPhaseMu] =
                state[detail::kPhasePolarVelocity];
            derivative[detail::kPhasePolarVelocity] =
                -polar_frequency_sq *
                state[detail::kPhaseMu];
            return derivative;
        };
    constexpr double step_size = 0.001;
    auto state = initial;
    double mino = 0.0;
    const auto first_step = numerics::dopri5_step(
        initial,
        0.0,
        step_size,
        polar_rhs,
        test_config());
    bool all_steps_accepted =
        first_step.accepted &&
        first_step.dense_output.has_value();

    const auto before_turn =
        detail::locate_kerr_turning_phase_crossing(
            detail::TurningCoordinate::Polar,
            potentials,
            *first_step.dense_output,
            step_size,
            1.0e-12,
            1.0e-12);
    check(
        "interval ending before the root has no crossing",
        !before_turn.has_value());

    std::optional<numerics::Dopri5DenseOutput<7>>
        crossing_dense;
    state = first_step.state;
    mino = step_size;
    for (std::size_t step_index = 1;
         step_index < 100 && !crossing_dense.has_value();
         ++step_index) {
        const auto phase_step = numerics::dopri5_step(
            state,
            mino,
            step_size,
            polar_rhs,
            test_config());
        all_steps_accepted =
            all_steps_accepted &&
            phase_step.accepted &&
            phase_step.dense_output.has_value();
        if (std::signbit(
                state[detail::kPhasePolarVelocity]) !=
            std::signbit(
                phase_step
                    .state[detail::kPhasePolarVelocity])) {
            crossing_dense = phase_step.dense_output;
        }
        state = phase_step.state;
        mino += step_size;
    }
    check(
        "synthetic polar phase steps are accepted",
        all_steps_accepted);
    const auto crossing =
        crossing_dense.has_value()
            ? detail::locate_kerr_turning_phase_crossing(
                  detail::TurningCoordinate::Polar,
                  potentials,
                  *crossing_dense,
                  crossing_dense->end(),
                  1.0e-12,
                  1.0e-12)
            : std::optional<detail::TurningPhaseCrossing>{};
    const double polar_maximum =
        std::sqrt(carter / polar_frequency_sq);
    const double expected_mino =
        std::acos(0.5 / polar_maximum) /
        std::sqrt(polar_frequency_sq);
    check(
        "polar velocity sign change locates a simple root",
        crossing.has_value() &&
            crossing->status ==
                detail::TurningStatus::Simple);
    check(
        "polar phase root time is accurate",
        crossing.has_value() &&
            std::fabs(
                crossing->mino_parameter -
                expected_mino) <
                1.0e-10);

    auto root_initial = initial;
    root_initial[detail::kPhaseMu] = polar_maximum;
    root_initial[detail::kPhasePolarVelocity] = 0.0;
    const auto root_step = numerics::dopri5_step(
        root_initial,
        0.0,
        step_size,
        polar_rhs,
        test_config());
    const auto duplicate =
        detail::locate_kerr_turning_phase_crossing(
            detail::TurningCoordinate::Polar,
            potentials,
            *root_step.dense_output,
            step_size,
            1.0e-12,
            1.0e-12);
    check(
        "root at step start is not counted twice",
        !duplicate.has_value());

    const auto critical_classification =
        crossing_dense.has_value()
            ? detail::locate_kerr_turning_phase_crossing(
                  detail::TurningCoordinate::Polar,
                  potentials,
                  *crossing_dense,
                  crossing_dense->end(),
                  1.0e-12,
                  10.0)
            : std::optional<detail::TurningPhaseCrossing>{};
    check(
        "crossing derivative gate identifies near-critical roots",
        critical_classification.has_value() &&
            critical_classification->status ==
                detail::TurningStatus::NearCritical);

    std::cout << "\n=== Results: " << passed
              << " passed, " << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
