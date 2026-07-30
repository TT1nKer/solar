#include "solar/relativity/radiative_transfer.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

using namespace solar::relativity;

namespace {

int passed = 0;
int failed = 0;
double maximum_solution_error = 0.0;
double maximum_subdivision_error = 0.0;

void check(
    const char* name,
    bool condition,
    double error = 0.0) {
    std::cout << (condition ? "  PASS: " : "  FAIL: ")
              << name << " (error=" << error << ")\n";
    condition ? ++passed : ++failed;
}

double normalized_error(double actual, double expected) {
    return std::fabs(actual - expected) /
           std::max({1.0, std::fabs(actual), std::fabs(expected)});
}

void check_near(
    const char* name,
    double actual,
    double expected,
    double tolerance) {
    const double error = normalized_error(actual, expected);
    maximum_solution_error =
        std::max(maximum_solution_error, error);
    check(
        name,
        std::isfinite(actual) && error <= tolerance,
        error);
}

BackwardTransferState advance_repeatedly(
    const TransferCoefficients& coefficients,
    double total_ds,
    int subdivisions) {
    BackwardTransferState state{};
    const double ds =
        total_ds / static_cast<double>(subdivisions);
    for (int index = 0; index < subdivisions; ++index) {
        const TransferAdvanceResult result =
            advance_backward_transfer(
                state, coefficients, ds);
        if (!result) {
            return BackwardTransferState{
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(),
            };
        }
        state = result.state;
    }
    return state;
}

void check_subdivision(
    int subdivisions,
    const BackwardTransferState& expected) {
    const BackwardTransferState actual =
        advance_repeatedly(
            TransferCoefficients{3.0, 0.4},
            5.0,
            subdivisions);
    const double error = std::max({
        normalized_error(
            actual.invariant_intensity,
            expected.invariant_intensity),
        normalized_error(
            actual.transmission,
            expected.transmission),
        normalized_error(
            actual.optical_depth,
            expected.optical_depth),
    });
    maximum_subdivision_error =
        std::max(maximum_subdivision_error, error);
    check(
        "constant solution is subdivision invariant",
        error <= 5.0e-14,
        error);
}

void check_error(
    const char* name,
    const BackwardTransferState& state,
    const TransferCoefficients& coefficients,
    double ds,
    TransferError expected) {
    const TransferAdvanceResult result =
        advance_backward_transfer(
            state, coefficients, ds);
    check(
        name,
        !result &&
            result.error == expected &&
            !result.message.empty());
}

} // namespace

int main() {
    const BackwardTransferState initial{};

    const TransferAdvanceResult constant_emission =
        advance_backward_transfer(
            initial,
            TransferCoefficients{2.5, 0.0},
            4.0);
    check("constant emission succeeds", bool(constant_emission));
    check_near(
        "constant emission analytic intensity",
        constant_emission.state.invariant_intensity,
        10.0,
        2.0e-14);
    check_near(
        "zero absorption preserves transmission",
        constant_emission.state.transmission,
        1.0,
        0.0);

    const double expected_transmission = std::exp(-2.0);
    const double expected_intensity =
        (3.0 / 0.4) * (1.0 - expected_transmission);
    const TransferAdvanceResult constant_source =
        advance_backward_transfer(
            initial,
            TransferCoefficients{3.0, 0.4},
            5.0);
    check("constant source succeeds", bool(constant_source));
    check_near(
        "constant source analytic intensity",
        constant_source.state.invariant_intensity,
        expected_intensity,
        5.0e-14);
    check_near(
        "constant absorption analytic transmission",
        constant_source.state.transmission,
        expected_transmission,
        5.0e-15);
    check_near(
        "constant absorption optical depth",
        constant_source.state.optical_depth,
        2.0,
        0.0);

    for (const int subdivisions : {1, 2, 7, 100}) {
        check_subdivision(
            subdivisions, constant_source.state);
    }

    const TransferAdvanceResult optically_thin =
        advance_backward_transfer(
            initial,
            TransferCoefficients{2.0, 1.0e-12},
            1.0);
    const double thin_expected =
        2.0 * (-std::expm1(-1.0e-12)) / 1.0e-12;
    check_near(
        "optically thin branch retains precision",
        optically_thin.state.invariant_intensity,
        thin_expected,
        2.0e-15);

    const TransferAdvanceResult optically_thick =
        advance_backward_transfer(
            initial,
            TransferCoefficients{4.0, 2.0},
            500.0);
    check_near(
        "optically thick source reaches source function",
        optically_thick.state.invariant_intensity,
        2.0,
        2.0e-15);
    check_near(
        "optically thick transmission saturates",
        optically_thick.state.transmission,
        0.0,
        0.0);

    const BackwardTransferState foreground{
        1.25, 0.4, -std::log(0.4)};
    const TransferAdvanceResult background =
        advance_backward_transfer(
            foreground,
            TransferCoefficients{2.0, 0.5},
            3.0);
    const double background_factor =
        1.0 - std::exp(-1.5);
    check_near(
        "foreground-to-background contribution uses prior transmission",
        background.state.invariant_intensity,
        1.25 + 0.4 * (2.0 / 0.5) * background_factor,
        5.0e-14);
    check_near(
        "foreground transmission attenuates farther material",
        background.state.transmission,
        0.4 * std::exp(-1.5),
        5.0e-15);

    const TransferAdvanceResult zero_length =
        advance_backward_transfer(
            foreground,
            TransferCoefficients{7.0, 9.0},
            0.0);
    check(
        "zero-length step is an exact identity",
        zero_length &&
            zero_length.state.invariant_intensity ==
                foreground.invariant_intensity &&
            zero_length.state.transmission ==
                foreground.transmission &&
            zero_length.state.optical_depth ==
                foreground.optical_depth);

    const TransferAdvanceResult vacuum =
        advance_backward_transfer(
            foreground,
            TransferCoefficients{},
            8.0);
    check(
        "vacuum step is an exact identity",
        vacuum &&
            vacuum.state.invariant_intensity ==
                foreground.invariant_intensity &&
            vacuum.state.transmission ==
                foreground.transmission &&
            vacuum.state.optical_depth ==
                foreground.optical_depth);

    check_near(
        "observer intensity restores frequency cube",
        specific_intensity_at_observer(2.0, 3.0),
        54.0,
        0.0);

    const double nan =
        std::numeric_limits<double>::quiet_NaN();
    check_error(
        "negative backward distance rejected",
        initial,
        TransferCoefficients{},
        -1.0,
        TransferError::InvalidStep);
    check_error(
        "non-finite backward distance rejected",
        initial,
        TransferCoefficients{},
        nan,
        TransferError::InvalidStep);
    check_error(
        "negative emissivity rejected",
        initial,
        TransferCoefficients{-1.0, 0.0},
        1.0,
        TransferError::InvalidCoefficients);
    check_error(
        "negative absorption rejected",
        initial,
        TransferCoefficients{0.0, -1.0},
        1.0,
        TransferError::InvalidCoefficients);
    check_error(
        "non-finite coefficient rejected",
        initial,
        TransferCoefficients{nan, 0.0},
        1.0,
        TransferError::InvalidCoefficients);
    check_error(
        "negative accumulated intensity rejected",
        BackwardTransferState{-1.0, 1.0, 0.0},
        TransferCoefficients{},
        1.0,
        TransferError::NonFiniteInput);
    check_error(
        "transmission above one rejected",
        BackwardTransferState{0.0, 1.1, 0.0},
        TransferCoefficients{},
        1.0,
        TransferError::NonFiniteInput);
    check_error(
        "inconsistent saturated state rejected",
        BackwardTransferState{
            0.0,
            0.2,
            std::numeric_limits<double>::infinity()},
        TransferCoefficients{},
        1.0,
        TransferError::NonFiniteInput);

    try {
        (void)specific_intensity_at_observer(1.0, -2.0);
        check("negative observer frequency rejected", false);
    } catch (const std::invalid_argument&) {
        check("negative observer frequency rejected", true);
    }

    std::cout.precision(17);
    std::cout
        << "  max_constant_solution_error="
        << maximum_solution_error
        << " max_subdivision_error="
        << maximum_subdivision_error
        << '\n';
    std::cout << "\n=== Results: " << passed
              << " passed, " << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
