#include "solar/relativity/fluid_model.h"
#include "solar/relativity/kerr_bl_metric.h"
#include "solar/relativity/kerr_chart_transform.h"
#include "solar/relativity/kerr_orbits.h"
#include "solar/relativity/kerr_schild_metric.h"
#include "solar/relativity/local_initialization.h"
#include "solar/relativity/minkowski_metric.h"
#include "solar/relativity/spacetime_algebra.h"
#include "solar/relativity/thin_disk.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace solar::relativity;

namespace {

int passed = 0;
int failed = 0;
double max_surface_normal_error = 0.0;
double max_surface_composition_error = 0.0;

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
    check(
        name,
        std::isfinite(actual) && error <= tolerance,
        error);
}

void check_near_and_record(
    const char* name,
    double actual,
    double expected,
    double tolerance,
    double& maximum_error) {
    const double error = normalized_error(actual, expected);
    maximum_error = std::max(maximum_error, error);
    check(
        name,
        std::isfinite(actual) && error <= tolerance,
        error);
}

AnalyticCircularDiskConfig disk_config() {
    return AnalyticCircularDiskConfig{
        1.0,
        0.5,
        OrbitSense::Prograde,
        6.0,
        20.0,
        1.0,
        16.407349347422414,
        0.0,
        1.0e-8,
    };
}

PhaseSpaceState photon_with_emitter_frequency(
    const KerrBoyerLindquistMetric& metric,
    const Contravariant4& point,
    double emitter_frequency,
    double theta_direction = -1.0) {
    const ObserverResult emitter =
        make_equatorial_circular_observer(
            metric, point, OrbitSense::Prograde);
    if (!emitter) {
        return {};
    }
    InitialStateResult photon =
        initialize_local_photon(
            metric,
            *emitter.frame,
            Vec3{{0.0, theta_direction, 0.0}});
    if (!photon) {
        return {};
    }
    photon.state->p.v =
        photon.state->p.v * emitter_frequency;
    photon.state->affine = -3.0;
    return *photon.state;
}

template <typename Action>
void check_invalid_argument(
    const char* name,
    Action action) {
    bool rejected = false;
    try {
        action();
    } catch (const std::invalid_argument&) {
        rejected = true;
    } catch (...) {
    }
    check(name, rejected);
}

void check_two_sheet_composition(
    const KerrBoyerLindquistMetric& metric,
    const Contravariant4& point) {
    const double optical_depth = std::log(2.0);
    ThinDiskCrossingRecorder recorder(
        ThinDiskRecorderConfig{
            DiskOpacityMode::SemiTransparent, 8},
        AnalyticCircularDiskFluid(disk_config()),
        ThinDiskSurfaceEmission(
            0.75, 10.0 / 4096.0, optical_depth));

    PhaseSpaceState first =
        photon_with_emitter_frequency(metric, point, 2.0);
    PhaseSpaceState second =
        photon_with_emitter_frequency(metric, point, 4.0);
    first.affine = -1.0;
    second.affine = -2.0;
    const ThinDiskRecordResult first_result =
        recorder.record(metric, first, 1.0);
    const ThinDiskRecordResult second_result =
        recorder.record(metric, second, 1.0);

    check(
        "two semi-transparent sheets remain open",
        first_result && first_result.recorded &&
            !first_result.closed &&
            second_result && second_result.recorded &&
            !second_result.closed &&
            !recorder.closed());
    check(
        "two semi-transparent crossings are ordered",
        recorder.crossings().size() == 2 &&
            recorder.crossings()[0].image_order == 0 &&
            recorder.crossings()[1].image_order == 1 &&
            recorder.crossings()[0].affine == -1.0 &&
            recorder.crossings()[1].affine == -2.0);
    check_near_and_record(
        "two-sheet specific intensity uses foreground transmission",
        recorder.observed().specific_intensity,
        0.3984375,
        3.0e-13,
        max_surface_composition_error);
    check_near_and_record(
        "two-sheet bolometric intensity uses foreground transmission",
        recorder.observed().bolometric_intensity,
        0.322265625,
        3.0e-13,
        max_surface_composition_error);
    check_near_and_record(
        "two-sheet transmission is multiplicative",
        recorder.observed().transmission,
        0.25,
        3.0e-13,
        max_surface_composition_error);
}

void check_crossing_bound_and_vacuum(
    const KerrBoyerLindquistMetric& metric,
    const Contravariant4& point) {
    ThinDiskCrossingRecorder bounded(
        ThinDiskRecorderConfig{
            DiskOpacityMode::SemiTransparent, 8},
        AnalyticCircularDiskFluid(disk_config()),
        ThinDiskSurfaceEmission(1.0, 1.0, 0.1));
    const PhaseSpaceState photon =
        photon_with_emitter_frequency(metric, point, 2.0);
    bool first_eight_recorded = true;
    for (std::size_t index = 0; index < 8; ++index) {
        PhaseSpaceState crossing = photon;
        crossing.affine = -static_cast<double>(index + 1);
        const ThinDiskRecordResult result =
            bounded.record(metric, crossing, 1.0);
        first_eight_recorded =
            first_eight_recorded &&
            result && result.recorded && !result.closed;
    }
    const ThinDiskRecordResult ninth =
        bounded.record(metric, photon, 1.0);
    check(
        "exactly eight semi-transparent crossings fit",
        first_eight_recorded &&
            bounded.crossings().size() == 8);
    check(
        "ninth crossing returns explicit limit error",
        !ninth &&
            ninth.error ==
                TransferError::CrossingLimitReached &&
            bounded.crossings().size() == 8);

    ThinDiskCrossingRecorder vacuum_then_disk(
        ThinDiskRecorderConfig{
            DiskOpacityMode::Opaque, 8},
        AnalyticCircularDiskFluid(disk_config()),
        ThinDiskSurfaceEmission(1.0, 1.0, 1.0));
    const Contravariant4 outside{
        Vec4{{0.0, 25.0, 0.5 * std::acos(-1.0), 0.2}}};
    const PhaseSpaceState outside_photon =
        photon_with_emitter_frequency(
            metric, outside, 2.0);
    const ThinDiskRecordResult ignored =
        vacuum_then_disk.record(
            metric, outside_photon, 1.0);
    const ThinDiskRecordResult accepted =
        vacuum_then_disk.record(metric, photon, 1.0);
    check(
        "out-of-radius sample does not consume image order",
        ignored && !ignored.recorded &&
            accepted && accepted.recorded &&
            vacuum_then_disk.crossings().size() == 1 &&
            vacuum_then_disk.crossings()[0].image_order == 0);
}

void check_chart_agreement(
    const KerrBoyerLindquistMetric& bl_metric,
    const Contravariant4& bl_point) {
    const KerrSchildCartesianMetric ks_metric(1.0, 0.5);
    const KerrChartTransform transform(1.0, 0.5);
    const PhaseSpaceState bl_photon =
        photon_with_emitter_frequency(
            bl_metric, bl_point, 2.0);
    const PhaseSpaceState ks_photon =
        transform.state_to_kerr_schild(bl_photon);
    const ThinDiskSurfaceEmission emission(
        0.75, 10.0 / 4096.0, 0.7);
    ThinDiskCrossingRecorder bl_recorder(
        ThinDiskRecorderConfig{
            DiskOpacityMode::Opaque, 8},
        AnalyticCircularDiskFluid(disk_config()),
        emission);
    ThinDiskCrossingRecorder ks_recorder(
        ThinDiskRecorderConfig{
            DiskOpacityMode::Opaque, 8},
        AnalyticCircularDiskFluid(disk_config()),
        emission);

    const ThinDiskRecordResult bl_result =
        bl_recorder.record(bl_metric, bl_photon, 1.0);
    const ThinDiskRecordResult ks_result =
        ks_recorder.record(ks_metric, ks_photon, 1.0);
    const bool both_recorded =
        bl_result && bl_result.recorded &&
        ks_result && ks_result.recorded &&
        bl_recorder.crossings().size() == 1 &&
        ks_recorder.crossings().size() == 1;
    check("BL and KS both record the same disk event", both_recorded);
    if (!both_recorded) {
        return;
    }

    const ThinDiskCrossing& bl =
        bl_recorder.crossings().front();
    const ThinDiskCrossing& ks =
        ks_recorder.crossings().front();
    check_near_and_record(
        "BL and KS disk radii agree",
        ks.disk_radius,
        bl.disk_radius,
        1.0e-10,
        max_surface_composition_error);
    check_near_and_record(
        "BL and KS redshifts agree",
        ks.redshift_g,
        bl.redshift_g,
        1.0e-10,
        max_surface_composition_error);

    const Vec4 expected_ks_normal = multiply(
        transform.boyer_lindquist_to_kerr_schild_jacobian(
            bl_point),
        bl.surface_normal.v);
    double normal_component_error = 0.0;
    for (std::size_t component = 0; component < 4; ++component) {
        normal_component_error = std::max(
            normal_component_error,
            std::fabs(
                ks.surface_normal.v[component] -
                expected_ks_normal[component]));
    }
    max_surface_normal_error = std::max(
        max_surface_normal_error, normal_component_error);
    check(
        "KS normal is the full contravariant BL transform",
        normal_component_error <= 1.0e-10,
        normal_component_error);
    check_near_and_record(
        "KS normal is unit spacelike",
        metric_inner_product(
            ks_metric.covariant(ks.position),
            ks.surface_normal,
            ks.surface_normal),
        1.0,
        1.0e-10,
        max_surface_normal_error);
    check_near_and_record(
        "KS normal remains orthogonal to emitter",
        metric_inner_product(
            ks_metric.covariant(ks.position),
            ks.surface_normal,
            ks.emitter_four_velocity),
        0.0,
        1.0e-10,
        max_surface_normal_error);
    check(
        "BL and KS face classification agrees",
        bl.front_facing == ks.front_facing);
}

void check_invalid_configuration(
    const KerrBoyerLindquistMetric& metric,
    const Contravariant4& point) {
    check_invalid_argument(
        "surface rejects negative specific scale",
        [] {
            ThinDiskSurfaceEmission emission(-1.0, 0.0, 0.0);
        });
    check_invalid_argument(
        "surface rejects non-finite bolometric scale",
        [] {
            ThinDiskSurfaceEmission emission(
                0.0,
                std::numeric_limits<double>::infinity(),
                0.0);
        });
    check_invalid_argument(
        "surface rejects negative optical depth",
        [] {
            ThinDiskSurfaceEmission emission(0.0, 0.0, -1.0);
        });
    check_invalid_argument(
        "recorder rejects zero crossing bound",
        [] {
            ThinDiskCrossingRecorder recorder(
                ThinDiskRecorderConfig{
                    DiskOpacityMode::Opaque, 0},
                AnalyticCircularDiskFluid(disk_config()),
                ThinDiskSurfaceEmission(1.0, 1.0, 1.0));
        });
    check_invalid_argument(
        "recorder rejects unknown opacity mode",
        [] {
            ThinDiskCrossingRecorder recorder(
                ThinDiskRecorderConfig{
                    static_cast<DiskOpacityMode>(99), 8},
                AnalyticCircularDiskFluid(disk_config()),
                ThinDiskSurfaceEmission(1.0, 1.0, 1.0));
        });
    check_invalid_argument(
        "semi-transparent mode rejects infinite sheet depth",
        [] {
            ThinDiskCrossingRecorder recorder(
                ThinDiskRecorderConfig{
                    DiskOpacityMode::SemiTransparent, 8},
                AnalyticCircularDiskFluid(disk_config()),
                ThinDiskSurfaceEmission(
                    1.0,
                    1.0,
                    std::numeric_limits<double>::infinity()));
        });

    ThinDiskCrossingRecorder overflow(
        ThinDiskRecorderConfig{
            DiskOpacityMode::Opaque, 8},
        AnalyticCircularDiskFluid(disk_config()),
        ThinDiskSurfaceEmission(
            std::numeric_limits<double>::max(),
            0.0,
            1.0));
    const ThinDiskRecordResult overflow_result =
        overflow.record(
            metric,
            photon_with_emitter_frequency(
                metric, point, 2.0),
            1.0);
    check(
        "surface overflow returns explicit non-finite error",
        !overflow_result &&
            overflow_result.error ==
                TransferError::NonFiniteResult &&
            overflow.crossings().empty());
}

void check_runtime_failure_meaning(
    const KerrBoyerLindquistMetric& metric,
    const Contravariant4& point) {
    ThinDiskCrossingRecorder recorder(
        ThinDiskRecorderConfig{
            DiskOpacityMode::SemiTransparent, 8},
        AnalyticCircularDiskFluid(disk_config()),
        ThinDiskSurfaceEmission(1.0, 1.0, 0.2));
    const PhaseSpaceState valid_photon =
        photon_with_emitter_frequency(
            metric, point, 2.0);
    const ThinDiskRecordResult accepted =
        recorder.record(metric, valid_photon, 1.0);
    const ThinDiskObservedState before_failure =
        recorder.observed();

    PhaseSpaceState non_finite = valid_photon;
    non_finite.x.v[0] =
        std::numeric_limits<double>::quiet_NaN();
    const ThinDiskRecordResult non_finite_result =
        recorder.record(metric, non_finite, 1.0);
    check(
        "non-finite surface photon has context error",
        !non_finite_result &&
            non_finite_result.error ==
                TransferError::NonFiniteInput);

    PhaseSpaceState invalid_point = valid_photon;
    invalid_point.x.v[1] = 1.0;
    const ThinDiskRecordResult invalid_point_result =
        recorder.record(metric, invalid_point, 1.0);
    check(
        "invalid surface metric point has context error",
        !invalid_point_result &&
            invalid_point_result.error ==
                TransferError::InvalidMetricPoint);

    const ThinDiskRecordResult invalid_frequency =
        recorder.record(metric, valid_photon, 0.0);
    check(
        "invalid surface observer frequency is explicit",
        !invalid_frequency &&
            invalid_frequency.error ==
                TransferError::InvalidObserverFrequency);
    check(
        "failed crossings preserve accepted state",
        accepted && accepted.recorded &&
            recorder.crossings().size() == 1 &&
            recorder.observed().specific_intensity ==
                before_failure.specific_intensity &&
            recorder.observed().bolometric_intensity ==
                before_failure.bolometric_intensity &&
            recorder.observed().transmission ==
                before_failure.transmission);

    const KerrBoyerLindquistMetric mismatched_metric(
        1.0, 0.4);
    const ThinDiskRecordResult mismatched =
        recorder.record(
            mismatched_metric,
            photon_with_emitter_frequency(
                mismatched_metric, point, 2.0),
            1.0);
    check(
        "mismatched Kerr parameters remain fluid failure",
        !mismatched &&
            mismatched.error ==
                TransferError::InvalidFluidSample);

    const MinkowskiMetric minkowski;
    PhaseSpaceState flat_photon = valid_photon;
    flat_photon.x =
        Contravariant4{Vec4{{0.0, 8.0, 0.0, 0.0}}};
    const ThinDiskRecordResult unsupported =
        recorder.record(minkowski, flat_photon, 1.0);
    check(
        "unsupported metric remains fluid failure",
        !unsupported &&
            unsupported.error ==
                TransferError::InvalidFluidSample);
}

} // namespace

int main() {
    constexpr double pi = 3.14159265358979323846;
    const KerrBoyerLindquistMetric metric(1.0, 0.5);
    const Contravariant4 point{
        Vec4{{0.0, 8.0, 0.5 * pi, 0.2}}};
    const PhaseSpaceState photon =
        photon_with_emitter_frequency(
            metric, point, 2.0);
    const AnalyticCircularDiskFluid disk(disk_config());
    const ThinDiskSurfaceEmission emission(
        0.75,
        10.0 / 4096.0,
        0.7);
    ThinDiskCrossingRecorder recorder(
        ThinDiskRecorderConfig{
            DiskOpacityMode::Opaque, 8},
        disk,
        emission);

    const ThinDiskRecordResult result =
        recorder.record(metric, photon, 1.0);
    check("opaque disk crossing succeeds", bool(result));
    check(
        "opaque disk crossing is recorded and closes",
        result.recorded && result.closed &&
            recorder.closed());
    check(
        "opaque disk stores one crossing",
        recorder.crossings().size() == 1);

    const ThinDiskCrossing& crossing =
        recorder.crossings().front();
    check_near_and_record(
        "surface emitter frequency",
        crossing.emitter_frequency,
        2.0,
        2.0e-13,
        max_surface_composition_error);
    check_near_and_record(
        "surface redshift",
        crossing.redshift_g,
        0.5,
        2.0e-13,
        max_surface_composition_error);
    check_near_and_record(
        "surface observed temperature uses g",
        crossing.observed_temperature,
        4.0,
        2.0e-13,
        max_surface_composition_error);
    check_near_and_record(
        "surface specific intensity uses g cubed",
        crossing.observed_specific_intensity,
        0.75,
        2.0e-13,
        max_surface_composition_error);
    check_near_and_record(
        "surface bolometric intensity uses g fourth",
        crossing.observed_bolometric_intensity,
        0.625,
        2.0e-13,
        max_surface_composition_error);
    check(
        "first crossing has image order zero",
        crossing.image_order == 0);
    check(
        "north-going photon is front-facing",
        crossing.front_facing);
    check_near(
        "crossing retains affine",
        crossing.affine,
        -3.0,
        0.0);
    check_near(
        "crossing retains disk radius",
        crossing.disk_radius,
        8.0,
        2.0e-14);

    const Mat4 covariant = metric.covariant(point);
    check_near_and_record(
        "surface normal is unit spacelike",
        metric_inner_product(
            covariant,
            crossing.surface_normal,
            crossing.surface_normal),
        1.0,
        1.0e-10,
        max_surface_normal_error);
    check_near_and_record(
        "surface normal is orthogonal to emitter",
        metric_inner_product(
            covariant,
            crossing.surface_normal,
            crossing.emitter_four_velocity),
        0.0,
        1.0e-10,
        max_surface_normal_error);
    check_near_and_record(
        "opaque observed specific intensity is full source",
        recorder.observed().specific_intensity,
        0.75,
        2.0e-13,
        max_surface_composition_error);
    check_near_and_record(
        "opaque observed bolometric intensity is full source",
        recorder.observed().bolometric_intensity,
        0.625,
        2.0e-13,
        max_surface_composition_error);
    check_near_and_record(
        "opaque crossing blocks farther surfaces",
        recorder.observed().transmission,
        0.0,
        0.0,
        max_surface_composition_error);

    const ObserverResult expected_emitter =
        make_equatorial_circular_observer(
            metric, point, OrbitSense::Prograde);
    for (std::size_t component = 0; component < 4; ++component) {
        check_near(
            "crossing retains position component",
            crossing.position.v[component],
            point.v[component],
            0.0);
        check_near(
            "crossing retains emitter velocity component",
            crossing.emitter_four_velocity.v[component],
            expected_emitter.frame->tetrad.basis[0].v[component],
            2.0e-13);
    }

    ThinDiskCrossingRecorder back_face_recorder(
        ThinDiskRecorderConfig{
            DiskOpacityMode::Opaque, 8},
        AnalyticCircularDiskFluid(disk_config()),
        emission);
    const ThinDiskRecordResult back_face_result =
        back_face_recorder.record(
            metric,
            photon_with_emitter_frequency(
                metric, point, 2.0, 1.0),
            1.0);
    check(
        "south-going photon is back-facing",
        back_face_result &&
            back_face_result.recorded &&
            !back_face_recorder.crossings()
                 .front()
                 .front_facing);

    check_two_sheet_composition(metric, point);
    check_crossing_bound_and_vacuum(metric, point);
    check_chart_agreement(metric, point);
    check_invalid_configuration(metric, point);
    check_runtime_failure_meaning(metric, point);

    std::cout.precision(17);
    std::cout << "  max_surface_normal_error="
              << max_surface_normal_error
              << " max_surface_composition_error="
              << max_surface_composition_error << '\n';
    std::cout << "\n=== Results: " << passed
              << " passed, " << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
