#include "solar/relativity/emission_model.h"
#include "solar/relativity/fluid_model.h"
#include "solar/relativity/minkowski_metric.h"
#include "solar/relativity/radiative_transfer.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

using namespace solar::relativity;

namespace {

int passed = 0;
int failed = 0;
double maximum_redshift_error = 0.0;

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
    maximum_redshift_error =
        std::max(maximum_redshift_error, error);
    check(
        name,
        std::isfinite(actual) && error <= tolerance,
        error);
}

class FixedFluid final : public FluidModel {
public:
    explicit FixedFluid(FluidSample sample)
        : sample_(sample) {}

    FluidSample sample(
        const Metric&,
        const Contravariant4&) const override {
        return sample_;
    }

private:
    FluidSample sample_;
};

PhaseSpaceState photon_state() {
    return PhaseSpaceState{
        0.0,
        Contravariant4{Vec4{{0.0, 0.0, 0.0, 0.0}}},
        Covariant4{Vec4{{-1.0, 1.0, 0.0, 0.0}}},
    };
}

template <typename Callable>
void check_invalid_argument(
    const char* name,
    Callable&& callable) {
    try {
        callable();
        check(name, false);
    } catch (const std::invalid_argument&) {
        check(name, true);
    }
}

} // namespace

int main() {
    const MinkowskiMetric metric;
    const PhaseSpaceState photon = photon_state();
    constexpr double speed = 0.3;
    const double gamma =
        1.0 / std::sqrt(1.0 - speed * speed);
    const FluidSample moving_sample{
        true,
        2.0,
        5.0,
        Contravariant4{
            Vec4{{gamma, gamma * speed, 0.0, 0.0}}},
    };
    const FixedFluid moving_fluid(moving_sample);
    constexpr double observer_frequency = 230.0e9;
    const double expected_normalized_frequency =
        gamma * (1.0 - speed);
    const double expected_emitter_frequency =
        observer_frequency * expected_normalized_frequency;

    const RedshiftResult redshift = evaluate_redshift(
        metric,
        photon,
        moving_sample.four_velocity,
        observer_frequency);
    check("boosted-emitter redshift succeeds", bool(redshift));
    check_near(
        "boosted-emitter normalized frequency",
        redshift.sample.normalized_emitter_frequency,
        expected_normalized_frequency,
        2.0e-14);
    check_near(
        "boosted-emitter physical frequency",
        redshift.sample.emitter_frequency,
        expected_emitter_frequency,
        2.0e-14);
    check_near(
        "boosted-emitter redshift factor",
        redshift.sample.redshift_g,
        1.0 / expected_normalized_frequency,
        2.0e-14);

    const GreyEmission grey(3.0, 0.4);
    const TransferEvaluationResult grey_result =
        evaluate_transfer_coefficients(
            metric,
            photon,
            observer_frequency,
            moving_fluid,
            grey);
    check("grey emission succeeds", bool(grey_result));
    check_near(
        "grey invariant emissivity",
        grey_result.coefficients.invariant_emissivity,
        6.0 /
            (expected_emitter_frequency *
             expected_emitter_frequency),
        2.0e-14);
    check_near(
        "grey invariant absorption",
        grey_result.coefficients.invariant_absorption,
        0.8 * expected_emitter_frequency,
        2.0e-14);

    const VacuumFluid vacuum_fluid;
    const TransferEvaluationResult vacuum_material =
        evaluate_transfer_coefficients(
            metric,
            photon,
            observer_frequency,
            vacuum_fluid,
            grey);
    check(
        "vacuum fluid returns exact zero coefficients",
        vacuum_material &&
            vacuum_material.coefficients.invariant_emissivity == 0.0 &&
            vacuum_material.coefficients.invariant_absorption == 0.0);

    const VacuumEmission vacuum_emission;
    const TransferEvaluationResult vacuum_radiation =
        evaluate_transfer_coefficients(
            metric,
            photon,
            observer_frequency,
            moving_fluid,
            vacuum_emission);
    check(
        "vacuum emission returns exact zero coefficients",
        vacuum_radiation &&
            vacuum_radiation.coefficients.invariant_emissivity == 0.0 &&
            vacuum_radiation.coefficients.invariant_absorption == 0.0);
    check_near(
        "vacuum emission still records physical redshift",
        vacuum_radiation.redshift.redshift_g,
        1.0 / expected_normalized_frequency,
        2.0e-14);

    const DebugPaintEmission debug(1.25, 0.75);
    const TransferEvaluationResult debug_result =
        evaluate_transfer_coefficients(
            metric,
            photon,
            observer_frequency,
            moving_fluid,
            debug);
    check(
        "debug emission preserves configured invariant coefficients",
        debug_result &&
            debug_result.coefficients.invariant_emissivity == 1.25 &&
            debug_result.coefficients.invariant_absorption == 0.75);

    const FluidSample returned_vacuum =
        vacuum_fluid.sample(metric, photon.x);
    check(
        "vacuum fluid sample is explicitly invalid material",
        !returned_vacuum.valid);

    const double nan =
        std::numeric_limits<double>::quiet_NaN();
    check_invalid_argument(
        "negative grey emissivity rejected",
        [] { (void)GreyEmission(-1.0, 0.0); });
    check_invalid_argument(
        "negative grey absorption rejected",
        [] { (void)GreyEmission(0.0, -1.0); });
    check_invalid_argument(
        "non-finite grey coefficient rejected",
        [nan] { (void)GreyEmission(nan, 0.0); });
    check_invalid_argument(
        "negative debug coefficient rejected",
        [] { (void)DebugPaintEmission(0.0, -1.0); });

    std::cout.precision(17);
    std::cout << "  max_redshift_error="
              << maximum_redshift_error << '\n';
    std::cout << "\n=== Results: " << passed
              << " passed, " << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
