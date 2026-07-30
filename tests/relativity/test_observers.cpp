#include "solar/relativity/kerr_bl_metric.h"
#include "solar/relativity/minkowski_metric.h"
#include "solar/relativity/observer.h"
#include "solar/relativity/spacetime_algebra.h"

#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>

using namespace solar::relativity;

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

void check_near(
    const std::string& name,
    double actual,
    double expected,
    double tolerance) {
    check(name, std::fabs(actual - expected) <= tolerance);
}

double spatial_determinant(const Tetrad& tetrad) {
    const auto& a = tetrad.basis[1].v;
    const auto& b = tetrad.basis[2].v;
    const auto& c = tetrad.basis[3].v;
    return a[1] * (b[2] * c[3] - b[3] * c[2]) -
           a[2] * (b[1] * c[3] - b[3] * c[1]) +
           a[3] * (b[1] * c[2] - b[2] * c[1]);
}

} // namespace

int main() {
    const MinkowskiMetric minkowski;
    const Contravariant4 origin{
        Vec4{{0.0, 0.0, 0.0, 0.0}}};

    const ObserverResult stationary =
        make_static_observer(minkowski, origin);
    check("Minkowski static observer exists", bool(stationary));
    check(
        "successful observer has no error",
        stationary.error == ObserverError::None);
    check_near(
        "static tetrad orthonormal",
        tetrad_orthonormality_error(
            minkowski, *stationary.frame),
        0.0,
        1.0e-15);
    check_near(
        "static observer time leg",
        stationary.frame->tetrad.basis[0].v[0],
        1.0,
        0.0);
    check(
        "static observer is right handed",
        spatial_determinant(stationary.frame->tetrad) > 0.0);

    const Vec4 local{{1.2, -0.5, 0.25, 2.0}};
    const Contravariant4 coordinate =
        tetrad_to_coordinate(stationary.frame->tetrad, local);
    const Vec4 local_round_trip = coordinate_to_tetrad(
        minkowski.covariant(origin),
        stationary.frame->tetrad,
        coordinate);
    for (std::size_t component = 0; component < 4; ++component) {
        check_near(
            "static tetrad local round trip",
            local_round_trip[component],
            local[component],
            1.0e-15);
    }

    const double gamma = 1.25;
    const Contravariant4 boosted_velocity{
        Vec4{{gamma, 0.75, 0.0, 0.0}}};
    const LookAtAttitude attitude{
        Contravariant4{
            Vec4{{0.0, 0.0, 0.0, -2.0}}},
        Contravariant4{
            Vec4{{0.0, 0.0, 3.0, 0.0}}},
    };
    const ObserverResult look_at =
        make_look_at_observer(
            minkowski,
            origin,
            boosted_velocity,
            attitude);
    check("boosted look-at observer exists", bool(look_at));
    check(
        "boosted tetrad meets gate",
        tetrad_orthonormality_error(
            minkowski, *look_at.frame) < 2.0e-15);
    check(
        "look direction is preserved",
        metric_inner_product(
            minkowski.covariant(origin),
            look_at.frame->tetrad.basis[1],
            attitude.look_direction) > 0.0);
    check(
        "up reference is preserved",
        metric_inner_product(
            minkowski.covariant(origin),
            look_at.frame->tetrad.basis[2],
            attitude.up_reference) > 0.0);
    check(
        "look-at observer is right handed",
        spatial_determinant(look_at.frame->tetrad) > 0.0);

    const Vec4 boosted_local{{0.8, -1.0, 0.4, 0.2}};
    const Contravariant4 boosted_coordinate =
        tetrad_to_coordinate(
            look_at.frame->tetrad, boosted_local);
    const Vec4 boosted_round_trip = coordinate_to_tetrad(
        minkowski.covariant(origin),
        look_at.frame->tetrad,
        boosted_coordinate);
    for (std::size_t component = 0; component < 4; ++component) {
        check_near(
            "boosted tetrad local round trip",
            boosted_round_trip[component],
            boosted_local[component],
            2.0e-15);
    }

    const std::array<Contravariant4, 3> coordinate_seeds{{
        Contravariant4{Vec4{{0.0, 1.0, 0.0, 0.0}}},
        Contravariant4{Vec4{{0.0, 0.0, 1.0, 0.0}}},
        Contravariant4{Vec4{{0.0, 0.0, 0.0, 1.0}}},
    }};
    const ObserverResult arbitrary =
        make_arbitrary_observer(
            minkowski,
            origin,
            boosted_velocity,
            coordinate_seeds);
    check("arbitrary boosted observer exists", bool(arbitrary));
    check(
        "arbitrary boosted observer meets gate",
        tetrad_orthonormality_error(
            minkowski, *arbitrary.frame) < 2.0e-15);

    const Contravariant4 non_unit_velocity{
        Vec4{{1.0, 0.5, 0.0, 0.0}}};
    check(
        "non-unit four-velocity rejected",
        make_look_at_observer(
            minkowski,
            origin,
            non_unit_velocity,
            attitude).error ==
            ObserverError::FourVelocityNotUnitTimelike);

    const LookAtAttitude parallel_attitude{
        attitude.look_direction,
        attitude.look_direction,
    };
    check(
        "parallel look and up rejected",
        make_look_at_observer(
            minkowski,
            origin,
            boosted_velocity,
            parallel_attitude).error ==
            ObserverError::DegenerateSpatialSeed);

    auto degenerate_seeds = coordinate_seeds;
    degenerate_seeds[1] = degenerate_seeds[0];
    check(
        "degenerate arbitrary seed rejected",
        make_arbitrary_observer(
            minkowski,
            origin,
            boosted_velocity,
            degenerate_seeds).error ==
            ObserverError::DegenerateSpatialSeed);

    Contravariant4 non_finite_position = origin;
    non_finite_position.v[1] =
        std::numeric_limits<double>::quiet_NaN();
    check(
        "non-finite observer position rejected",
        make_static_observer(
            minkowski, non_finite_position).error ==
            ObserverError::NonFiniteInput);

    const KerrBoyerLindquistMetric kerr(1.0, 0.8);
    const Contravariant4 ergosphere_point{
        Vec4{{0.0, 1.8, 1.5707963267948966, 0.0}}};
    check(
        "ergosphere test point is in metric domain",
        kerr.valid_point(ergosphere_point));
    check(
        "static observer rejected inside ergosphere",
        make_static_observer(
            kerr, ergosphere_point).error ==
            ObserverError::StaticWorldlineNotTimelike);

    const Contravariant4 horizon_point{
        Vec4{{
            0.0,
            kerr.outer_horizon_radius(),
            1.5707963267948966,
            0.0,
        }}};
    check(
        "invalid metric point rejected before construction",
        make_static_observer(
            kerr, horizon_point).error ==
            ObserverError::InvalidMetricPoint);

    std::cout << "\n=== Results: " << passed
              << " passed, " << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
