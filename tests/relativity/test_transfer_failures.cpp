#include "solar/relativity/emission_model.h"
#include "solar/relativity/fluid_model.h"
#include "solar/relativity/minkowski_metric.h"
#include "solar/relativity/radiative_transfer.h"
#include "solar/relativity/schwarzschild_metric.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace solar::relativity;

namespace {

int passed = 0;
int failed = 0;

void check(const char* name, bool condition) {
    std::cout << (condition ? "  PASS: " : "  FAIL: ")
              << name << '\n';
    condition ? ++passed : ++failed;
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

class ThrowingFluid final : public FluidModel {
public:
    FluidSample sample(
        const Metric&,
        const Contravariant4&) const override {
        throw std::runtime_error("fluid failure");
    }
};

class ThrowingEmission final : public EmissionModel {
public:
    TransferCoefficients coefficients(
        const FluidSample&,
        const Covariant4&,
        double) const override {
        throw std::runtime_error("emission failure");
    }
};

class FixedEmission final : public EmissionModel {
public:
    explicit FixedEmission(TransferCoefficients coefficients)
        : coefficients_(coefficients) {}

    TransferCoefficients coefficients(
        const FluidSample&,
        const Covariant4&,
        double) const override {
        return coefficients_;
    }

private:
    TransferCoefficients coefficients_;
};

PhaseSpaceState future_photon() {
    return PhaseSpaceState{
        0.0,
        Contravariant4{Vec4{{0.0, 0.0, 0.0, 0.0}}},
        Covariant4{Vec4{{-1.0, 1.0, 0.0, 0.0}}},
    };
}

FluidSample unit_fluid() {
    return FluidSample{
        true,
        1.0,
        1.0,
        Contravariant4{Vec4{{1.0, 0.0, 0.0, 0.0}}},
    };
}

void check_evaluation_error(
    const char* name,
    const Metric& metric,
    const PhaseSpaceState& photon,
    double observer_frequency,
    const FluidModel& fluid,
    const EmissionModel& emission,
    TransferError expected) {
    const TransferEvaluationResult result =
        evaluate_transfer_coefficients(
            metric,
            photon,
            observer_frequency,
            fluid,
            emission);
    check(
        name,
        !result &&
            result.error == expected &&
            !result.message.empty());
}

void check_redshift_error(
    const char* name,
    const Metric& metric,
    const PhaseSpaceState& photon,
    const Contravariant4& velocity,
    double observer_frequency,
    TransferError expected) {
    const RedshiftResult result = evaluate_redshift(
        metric,
        photon,
        velocity,
        observer_frequency);
    check(
        name,
        !result &&
            result.error == expected &&
            !result.message.empty());
}

} // namespace

int main() {
    const MinkowskiMetric minkowski;
    const PhaseSpaceState photon = future_photon();
    const VacuumEmission vacuum;
    const FixedFluid valid_fluid(unit_fluid());
    const double nan =
        std::numeric_limits<double>::quiet_NaN();

    FluidSample ignored_vacuum{};
    ignored_vacuum.valid = false;
    ignored_vacuum.density = nan;
    ignored_vacuum.temperature = -1.0;
    ignored_vacuum.four_velocity.v[0] = nan;
    const FixedFluid ignored_vacuum_fluid(ignored_vacuum);
    const TransferEvaluationResult vacuum_result =
        evaluate_transfer_coefficients(
            minkowski,
            photon,
            1.0,
            ignored_vacuum_fluid,
            FixedEmission(
                TransferCoefficients{nan, nan}));
    check(
        "invalid fields in explicit vacuum are not read",
        vacuum_result &&
            vacuum_result.coefficients.invariant_emissivity == 0.0 &&
            vacuum_result.coefficients.invariant_absorption == 0.0);

    FluidSample negative_density = unit_fluid();
    negative_density.density = -1.0;
    check_evaluation_error(
        "negative density rejected",
        minkowski,
        photon,
        1.0,
        FixedFluid(negative_density),
        vacuum,
        TransferError::InvalidFluidSample);

    FluidSample nonfinite_temperature = unit_fluid();
    nonfinite_temperature.temperature = nan;
    check_evaluation_error(
        "non-finite temperature rejected",
        minkowski,
        photon,
        1.0,
        FixedFluid(nonfinite_temperature),
        vacuum,
        TransferError::InvalidFluidSample);

    FluidSample nonunit_velocity = unit_fluid();
    nonunit_velocity.four_velocity.v[0] = 2.0;
    check_evaluation_error(
        "non-unit timelike velocity rejected",
        minkowski,
        photon,
        1.0,
        FixedFluid(nonunit_velocity),
        vacuum,
        TransferError::FourVelocityNotUnitTimelike);

    FluidSample spacelike_velocity = unit_fluid();
    spacelike_velocity.four_velocity.v =
        Vec4{{0.0, 1.0, 0.0, 0.0}};
    check_evaluation_error(
        "spacelike velocity rejected",
        minkowski,
        photon,
        1.0,
        FixedFluid(spacelike_velocity),
        vacuum,
        TransferError::FourVelocityNotUnitTimelike);

    FluidSample nonfinite_velocity = unit_fluid();
    nonfinite_velocity.four_velocity.v[0] = nan;
    check_evaluation_error(
        "non-finite velocity rejected",
        minkowski,
        photon,
        1.0,
        FixedFluid(nonfinite_velocity),
        vacuum,
        TransferError::InvalidFluidSample);

    PhaseSpaceState past_photon = photon;
    past_photon.p.v[0] = 1.0;
    check_evaluation_error(
        "past-directed photon rejected",
        minkowski,
        past_photon,
        1.0,
        valid_fluid,
        vacuum,
        TransferError::NonFutureDirectedPhoton);

    check_evaluation_error(
        "zero observer frequency rejected",
        minkowski,
        photon,
        0.0,
        valid_fluid,
        vacuum,
        TransferError::InvalidObserverFrequency);
    check_evaluation_error(
        "non-finite observer frequency rejected",
        minkowski,
        photon,
        nan,
        valid_fluid,
        vacuum,
        TransferError::InvalidObserverFrequency);

    PhaseSpaceState nonfinite_photon = photon;
    nonfinite_photon.p.v[2] = nan;
    check_evaluation_error(
        "non-finite photon rejected",
        minkowski,
        nonfinite_photon,
        1.0,
        valid_fluid,
        vacuum,
        TransferError::NonFiniteInput);

    const SchwarzschildBoyerLindquistMetric schwarzschild(1.0);
    PhaseSpaceState invalid_point = photon;
    invalid_point.x.v =
        Vec4{{0.0, 1.5, 1.5707963267948966, 0.0}};
    check_evaluation_error(
        "invalid metric point rejected",
        schwarzschild,
        invalid_point,
        1.0,
        valid_fluid,
        vacuum,
        TransferError::InvalidMetricPoint);

    const ThrowingFluid throwing_fluid;
    check_evaluation_error(
        "fluid exception converted to explicit failure",
        minkowski,
        photon,
        1.0,
        throwing_fluid,
        vacuum,
        TransferError::InvalidFluidSample);

    const ThrowingEmission throwing_emission;
    check_evaluation_error(
        "emission exception converted to explicit failure",
        minkowski,
        photon,
        1.0,
        valid_fluid,
        throwing_emission,
        TransferError::EmissionModelFailure);

    check_evaluation_error(
        "negative emissivity output rejected",
        minkowski,
        photon,
        1.0,
        valid_fluid,
        FixedEmission(TransferCoefficients{-1.0, 0.0}),
        TransferError::InvalidCoefficients);
    check_evaluation_error(
        "non-finite absorption output rejected",
        minkowski,
        photon,
        1.0,
        valid_fluid,
        FixedEmission(TransferCoefficients{0.0, nan}),
        TransferError::InvalidCoefficients);

    check_redshift_error(
        "redshift rejects non-unit emitter",
        minkowski,
        photon,
        nonunit_velocity.four_velocity,
        1.0,
        TransferError::FourVelocityNotUnitTimelike);
    check_redshift_error(
        "redshift rejects past-directed photon",
        minkowski,
        past_photon,
        unit_fluid().four_velocity,
        1.0,
        TransferError::NonFutureDirectedPhoton);
    check_redshift_error(
        "redshift rejects invalid observer frequency",
        minkowski,
        photon,
        unit_fluid().four_velocity,
        -1.0,
        TransferError::InvalidObserverFrequency);

    PhaseSpaceState high_frequency_photon = photon;
    high_frequency_photon.p.v[0] = -2.0;
    high_frequency_photon.p.v[1] = 2.0;
    check_redshift_error(
        "redshift reports finite-range overflow explicitly",
        minkowski,
        high_frequency_photon,
        unit_fluid().four_velocity,
        std::numeric_limits<double>::max(),
        TransferError::NonFiniteResult);

    std::cout << "\n=== Results: " << passed
              << " passed, " << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
