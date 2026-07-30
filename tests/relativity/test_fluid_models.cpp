#include "solar/relativity/fluid_model.h"
#include "solar/relativity/kerr_bl_metric.h"
#include "solar/relativity/kerr_chart_transform.h"
#include "solar/relativity/kerr_orbits.h"
#include "solar/relativity/kerr_schild_metric.h"
#include "solar/relativity/minkowski_metric.h"
#include "solar/relativity/spacetime_algebra.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace solar::relativity;

namespace {

int passed = 0;
int failed = 0;
double maximum_velocity_norm_error = 0.0;
double maximum_chart_velocity_error = 0.0;

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

Contravariant4 apply_jacobian(
    const Mat4& jacobian,
    const Contravariant4& vector) {
    Contravariant4 result{};
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            result.v[row] +=
                jacobian[row][column] * vector.v[column];
        }
    }
    return result;
}

double maximum_component_error(
    const Contravariant4& actual,
    const Contravariant4& expected) {
    double maximum = 0.0;
    for (std::size_t component = 0;
         component < 4;
         ++component) {
        maximum = std::max(
            maximum,
            normalized_error(
                actual.v[component],
                expected.v[component]));
    }
    return maximum;
}

AnalyticCircularDiskConfig disk_config(
    double mass,
    double spin,
    OrbitSense sense) {
    return AnalyticCircularDiskConfig{
        mass,
        spin,
        sense,
        std::nullopt,
        30.0 * mass,
        4.0,
        10.0,
        1.5,
        1.0e-8,
    };
}

AnalyticOpticallyThinTorusConfig torus_config() {
    return AnalyticOpticallyThinTorusConfig{
        1.0,
        0.6,
        OrbitSense::Prograde,
        8.0,
        2.0,
        0.2,
        5.0,
        7.0,
        0.4,
        1.0e-4,
    };
}

Contravariant4 bl_point(double radius, double theta) {
    return Contravariant4{
        Vec4{{0.25, radius, theta, -0.4}}};
}

void check_unit_velocity(
    const char* name,
    const Metric& metric,
    const Contravariant4& point,
    const FluidSample& sample) {
    const double norm = metric_inner_product(
        metric.covariant(point),
        sample.four_velocity,
        sample.four_velocity);
    const double error = std::fabs(norm + 1.0);
    maximum_velocity_norm_error =
        std::max(maximum_velocity_norm_error, error);
    check(
        name,
        sample.valid &&
            sample.four_velocity.v[0] > 0.0 &&
            error < 1.0e-10,
        error);
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

template <typename Callable>
void check_domain_error(
    const char* name,
    Callable&& callable) {
    try {
        callable();
        check(name, false);
    } catch (const std::domain_error&) {
        check(name, true);
    }
}

} // namespace

int main() {
    constexpr double pi = 3.14159265358979323846;
    const KerrBoyerLindquistMetric metric(2.0, 0.5);
    const AnalyticCircularDiskFluid disk(
        disk_config(2.0, 0.5, OrbitSense::Prograde));
    const double inner = kerr_isco_radius(
        metric, OrbitSense::Prograde);
    check_near(
        "disk defaults inner radius to matching ISCO",
        disk.inner_radius(),
        inner,
        2.0e-15);
    check_near(
        "disk retains configured outer radius",
        disk.outer_radius(),
        60.0,
        0.0);

    const FluidSample inner_sample =
        disk.sample(metric, bl_point(inner, 0.5 * pi));
    check("disk inner edge is valid material", inner_sample.valid);
    check_near(
        "zero-torque inner edge has zero temperature",
        inner_sample.temperature,
        0.0,
        0.0);
    check_near(
        "inner edge density follows profile",
        inner_sample.density,
        4.0,
        2.0e-14);

    const double radius = 2.0 * inner;
    const FluidSample profile =
        disk.sample(metric, bl_point(radius, 0.5 * pi));
    const double expected_density =
        4.0 * std::pow(2.0, -1.5);
    const double expected_flux_shape =
        std::pow(2.0, -3.0) *
        (1.0 - std::sqrt(0.5));
    const double expected_temperature =
        10.0 * std::pow(expected_flux_shape, 0.25);
    check_near(
        "disk density uses literal radial power",
        profile.density,
        expected_density,
        2.0e-14);
    check_near(
        "disk temperature uses literal zero-torque shape",
        profile.temperature,
        expected_temperature,
        2.0e-14);
    check_unit_velocity(
        "disk BL velocity is future unit timelike",
        metric,
        bl_point(radius, 0.5 * pi),
        profile);

    check(
        "disk rejects radius below inner edge as vacuum",
        !disk.sample(
                 metric,
                 bl_point(
                     std::nextafter(
                         inner, -std::numeric_limits<double>::infinity()),
                     0.5 * pi))
             .valid);
    const double sub_photon_radius =
        metric.outer_horizon_radius() + 0.1;
    const Contravariant4 sub_photon_bl =
        bl_point(sub_photon_radius, 0.5 * pi);
    check(
        "disk treats sub-photon-orbit BL point as vacuum",
        !disk.sample(metric, sub_photon_bl).valid);
    const KerrChartTransform sub_photon_transform(2.0, 0.5);
    check(
        "disk treats sub-photon-orbit KS point as vacuum",
        !disk.sample(
                 KerrSchildCartesianMetric(2.0, 0.5),
                 sub_photon_transform.position_to_kerr_schild(
                     sub_photon_bl))
             .valid);
    check(
        "disk rejects radius above outer edge as vacuum",
        !disk.sample(
                 metric,
                 bl_point(61.0, 0.5 * pi))
             .valid);
    check(
        "disk rejects points away from surface as vacuum",
        !disk.sample(
                 metric,
                 bl_point(radius, 0.5 * pi + 1.0e-6))
             .valid);

    for (const double spin : {-0.7, 0.0, 0.6}) {
        for (const OrbitSense sense :
             {OrbitSense::Prograde, OrbitSense::Retrograde}) {
            const KerrBoyerLindquistMetric bl_metric(1.0, spin);
            const KerrSchildCartesianMetric ks_metric(1.0, spin);
            const AnalyticCircularDiskFluid fixture_disk(
                disk_config(1.0, spin, sense));
            const double fixture_radius =
                std::max(
                    fixture_disk.inner_radius() + 1.0,
                    8.0);
            const Contravariant4 fixture_bl =
                bl_point(fixture_radius, 0.5 * pi);
            const KerrChartTransform transform(1.0, spin);
            const Contravariant4 fixture_ks =
                transform.position_to_kerr_schild(fixture_bl);
            const FluidSample bl_sample =
                fixture_disk.sample(bl_metric, fixture_bl);
            const FluidSample ks_sample =
                fixture_disk.sample(ks_metric, fixture_ks);

            check(
                "common BL/KS disk samples are valid",
                bl_sample.valid && ks_sample.valid);
            check_near(
                "common BL/KS disk density agrees",
                ks_sample.density,
                bl_sample.density,
                1.0e-12);
            check_near(
                "common BL/KS disk temperature agrees",
                ks_sample.temperature,
                bl_sample.temperature,
                1.0e-12);
            check_unit_velocity(
                "fixture BL velocity is future unit timelike",
                bl_metric,
                fixture_bl,
                bl_sample);
            check_unit_velocity(
                "fixture KS velocity is future unit timelike",
                ks_metric,
                fixture_ks,
                ks_sample);

            const Mat4 jacobian =
                transform
                    .boyer_lindquist_to_kerr_schild_jacobian(
                        fixture_bl);
            const Contravariant4 expected_ks_velocity =
                apply_jacobian(
                    jacobian, bl_sample.four_velocity);
            const double chart_error =
                maximum_component_error(
                    ks_sample.four_velocity,
                    expected_ks_velocity);
            maximum_chart_velocity_error =
                std::max(
                    maximum_chart_velocity_error,
                    chart_error);
            check(
                "KS disk velocity uses full chart Jacobian",
                chart_error < 1.0e-10,
                chart_error);
        }
    }

    const MinkowskiMetric minkowski;
    check_domain_error(
        "disk rejects unsupported metric",
        [&] {
            (void)disk.sample(
                minkowski,
                Contravariant4{
                    Vec4{{0.0, radius, 0.0, 0.0}}});
        });
    const KerrBoyerLindquistMetric mismatched_metric(
        2.0, 0.4);
    check_domain_error(
        "disk rejects mismatched Kerr parameters",
        [&] {
            (void)disk.sample(
                mismatched_metric,
                bl_point(radius, 0.5 * pi));
        });

    const AnalyticOpticallyThinTorus torus(torus_config());
    const KerrBoyerLindquistMetric torus_bl_metric(1.0, 0.6);
    const KerrSchildCartesianMetric torus_ks_metric(1.0, 0.6);
    const Contravariant4 torus_center =
        bl_point(8.0, 0.5 * pi);
    const FluidSample torus_center_sample =
        torus.sample(torus_bl_metric, torus_center);
    check_near(
        "torus center density equals configured scale",
        torus_center_sample.density,
        5.0,
        2.0e-14);
    check_near(
        "torus center temperature equals configured scale",
        torus_center_sample.temperature,
        7.0,
        2.0e-14);
    check_unit_velocity(
        "torus center velocity is future unit timelike",
        torus_bl_metric,
        torus_center,
        torus_center_sample);

    const double one_width_shape = std::exp(-0.5);
    const FluidSample one_width =
        torus.sample(
            torus_bl_metric,
            bl_point(10.0, 0.5 * pi));
    check_near(
        "torus radial Gaussian is literal",
        one_width.density,
        5.0 * one_width_shape,
        2.0e-14);
    check_near(
        "torus temperature power is literal",
        one_width.temperature,
        7.0 * std::pow(one_width_shape, 0.4),
        2.0e-14);

    const FluidSample torus_north =
        torus.sample(
            torus_bl_metric,
            bl_point(8.0, 0.5 * pi - 0.1));
    const FluidSample torus_south =
        torus.sample(
            torus_bl_metric,
            bl_point(8.0, 0.5 * pi + 0.1));
    check_near(
        "torus density is north-south symmetric",
        torus_north.density,
        torus_south.density,
        2.0e-14);
    check_near(
        "torus temperature is north-south symmetric",
        torus_north.temperature,
        torus_south.temperature,
        2.0e-14);
    check(
        "torus compact cutoff returns vacuum",
        !torus.sample(
                  torus_bl_metric,
                  bl_point(20.0, 0.5 * pi))
             .valid);

    const Contravariant4 torus_bl =
        bl_point(8.5, 0.5 * pi - 0.05);
    const KerrChartTransform torus_transform(1.0, 0.6);
    const Contravariant4 torus_ks =
        torus_transform.position_to_kerr_schild(torus_bl);
    const FluidSample torus_bl_sample =
        torus.sample(torus_bl_metric, torus_bl);
    const FluidSample torus_ks_sample =
        torus.sample(torus_ks_metric, torus_ks);
    check(
        "common BL/KS torus samples are valid",
        torus_bl_sample.valid && torus_ks_sample.valid);
    check_near(
        "common BL/KS torus density agrees",
        torus_ks_sample.density,
        torus_bl_sample.density,
        1.0e-12);
    check_near(
        "common BL/KS torus temperature agrees",
        torus_ks_sample.temperature,
        torus_bl_sample.temperature,
        1.0e-12);
    check_unit_velocity(
        "torus BL velocity is future unit timelike",
        torus_bl_metric,
        torus_bl,
        torus_bl_sample);
    check_unit_velocity(
        "torus KS velocity is future unit timelike",
        torus_ks_metric,
        torus_ks,
        torus_ks_sample);
    const Contravariant4 expected_torus_ks_velocity =
        apply_jacobian(
            torus_transform
                .boyer_lindquist_to_kerr_schild_jacobian(
                    torus_bl),
            torus_bl_sample.four_velocity);
    const double torus_chart_error =
        maximum_component_error(
            torus_ks_sample.four_velocity,
            expected_torus_ks_velocity);
    maximum_chart_velocity_error =
        std::max(
            maximum_chart_velocity_error,
            torus_chart_error);
    check(
        "KS torus velocity uses full chart Jacobian",
        torus_chart_error < 1.0e-10,
        torus_chart_error);

    check_domain_error(
        "torus rejects unsupported metric",
        [&] {
            (void)torus.sample(
                minkowski,
                Contravariant4{
                    Vec4{{0.0, 8.0, 0.0, 0.0}}});
        });
    check_domain_error(
        "torus rejects mismatched Kerr parameters",
        [&] {
            (void)torus.sample(
                KerrBoyerLindquistMetric(1.0, 0.5),
                torus_center);
        });

    const double nan =
        std::numeric_limits<double>::quiet_NaN();
    check_invalid_argument(
        "disk rejects zero mass",
        [&] {
            auto config =
                disk_config(2.0, 0.5, OrbitSense::Prograde);
            config.mass_M = 0.0;
            (void)AnalyticCircularDiskFluid(config);
        });
    check_invalid_argument(
        "disk rejects extremal spin",
        [&] {
            auto config =
                disk_config(2.0, 0.5, OrbitSense::Prograde);
            config.spin_chi = 1.0;
            (void)AnalyticCircularDiskFluid(config);
        });
    check_invalid_argument(
        "disk rejects non-finite outer radius",
        [&] {
            auto config =
                disk_config(2.0, 0.5, OrbitSense::Prograde);
            config.outer_radius_M = nan;
            (void)AnalyticCircularDiskFluid(config);
        });
    check_invalid_argument(
        "disk rejects reversed radial bounds",
        [&] {
            auto config =
                disk_config(2.0, 0.5, OrbitSense::Prograde);
            config.inner_radius_M = 20.0;
            config.outer_radius_M = 10.0;
            (void)AnalyticCircularDiskFluid(config);
        });
    check_invalid_argument(
        "disk rejects negative density scale",
        [&] {
            auto config =
                disk_config(2.0, 0.5, OrbitSense::Prograde);
            config.density_scale = -1.0;
            (void)AnalyticCircularDiskFluid(config);
        });
    check_invalid_argument(
        "disk rejects zero temperature scale",
        [&] {
            auto config =
                disk_config(2.0, 0.5, OrbitSense::Prograde);
            config.temperature_scale = 0.0;
            (void)AnalyticCircularDiskFluid(config);
        });
    check_invalid_argument(
        "disk rejects negative density power",
        [&] {
            auto config =
                disk_config(2.0, 0.5, OrbitSense::Prograde);
            config.density_power = -0.1;
            (void)AnalyticCircularDiskFluid(config);
        });
    check_invalid_argument(
        "disk rejects zero surface tolerance",
        [&] {
            auto config =
                disk_config(2.0, 0.5, OrbitSense::Prograde);
            config.surface_height_tolerance = 0.0;
            (void)AnalyticCircularDiskFluid(config);
        });
    check_invalid_argument(
        "disk rejects unknown orbit sense",
        [&] {
            auto config =
                disk_config(2.0, 0.5, OrbitSense::Prograde);
            config.sense = static_cast<OrbitSense>(99);
            (void)AnalyticCircularDiskFluid(config);
        });
    check_invalid_argument(
        "torus rejects zero radial width",
        [&] {
            auto config = torus_config();
            config.radial_width_M = 0.0;
            (void)AnalyticOpticallyThinTorus(config);
        });
    check_invalid_argument(
        "torus rejects zero angular width",
        [&] {
            auto config = torus_config();
            config.angular_width = 0.0;
            (void)AnalyticOpticallyThinTorus(config);
        });
    check_invalid_argument(
        "torus rejects invalid cutoff fraction",
        [&] {
            auto config = torus_config();
            config.density_cutoff_fraction = 1.0;
            (void)AnalyticOpticallyThinTorus(config);
        });
    check_invalid_argument(
        "torus rejects negative temperature power",
        [&] {
            auto config = torus_config();
            config.temperature_power = -0.1;
            (void)AnalyticOpticallyThinTorus(config);
        });

    std::cout.precision(17);
    std::cout << "  max_bl_ks_velocity_error="
              << maximum_chart_velocity_error
              << " max_velocity_norm_error="
              << maximum_velocity_norm_error
              << '\n';
    std::cout << "\n=== Results: " << passed
              << " passed, " << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
