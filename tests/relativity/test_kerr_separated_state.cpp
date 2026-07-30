#include "solar/relativity/hamiltonian.h"
#include "solar/relativity/kerr_bl_metric.h"

#include "../../src/relativity/kerr_separated_state.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace solar::relativity;

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

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

PhaseSpaceState constrained_state(
    const KerrBoyerLindquistMetric& metric,
    GeodesicKind kind,
    double radius,
    double theta,
    double energy,
    double lz,
    double p_theta,
    double radial_sign) {
    PhaseSpaceState state{
        0.25,
        Contravariant4{
            Vec4{{1.5, radius, theta, -0.75}}},
        Covariant4{
            Vec4{{-energy, 0.0, p_theta, lz}}},
    };
    const Mat4 inverse = metric.contravariant(state.x);
    double non_radial_twice_hamiltonian = 0.0;
    for (std::size_t row = 0; row < 4; ++row) {
        if (row == 1) {
            continue;
        }
        for (std::size_t column = 0; column < 4; ++column) {
            if (column == 1) {
                continue;
            }
            non_radial_twice_hamiltonian +=
                inverse[row][column] *
                state.p.v[row] *
                state.p.v[column];
        }
    }
    const double radial_squared =
        (2.0 * hamiltonian_target(kind) -
         non_radial_twice_hamiltonian) /
        inverse[1][1];
    if (!(radial_squared >= 0.0) ||
        !std::isfinite(radial_squared)) {
        throw std::domain_error(
            "test fixture has no real radial momentum");
    }
    state.p.v[1] =
        std::copysign(std::sqrt(radial_squared), radial_sign);
    return state;
}

template <typename Exception, typename Function>
void check_throws(
    const std::string& name,
    Function&& function) {
    try {
        function();
        check(name, false);
    } catch (const Exception&) {
        check(name, true);
    }
}

} // namespace

int main() {
    const KerrBoyerLindquistMetric kerr(1.0, 0.5);
    const PhaseSpaceState initial = constrained_state(
        kerr,
        GeodesicKind::Null,
        8.0,
        1.1,
        1.0,
        2.0,
        1.0,
        -1.0);
    check(
        "round-trip fixture satisfies Hamiltonian constraint",
        hamiltonian_constraint_error(
            kerr, initial, GeodesicKind::Null) <
            1.0e-14);

    const auto separated =
        detail::initialize_kerr_separated_state(
            kerr,
            initial,
            GeodesicKind::Null,
            1.0e-12,
            1.0e-10,
            1.0);
    const PhaseSpaceState rebuilt =
        detail::reconstruct_kerr_phase_space(
            kerr, separated.constants, separated.state);

    check_near(
        "round-trip affine",
        rebuilt.affine,
        initial.affine,
        1.0e-14);
    for (std::size_t component = 0;
         component < 4;
         ++component) {
        check_near(
            "round-trip coordinate " +
                std::to_string(component),
            rebuilt.x.v[component],
            initial.x.v[component],
            1.0e-12);
        check_near(
            "round-trip momentum " +
                std::to_string(component),
            rebuilt.p.v[component],
            initial.p.v[component],
            1.0e-11);
    }
    check(
        "inward radial direction",
        separated.state.radial_direction ==
            detail::SeparatedDirection::Negative);
    check(
        "mu direction follows negative sine theta-dot",
        separated.state.polar_direction ==
            detail::SeparatedDirection::Negative);

    const PhaseSpaceState equatorial = constrained_state(
        kerr,
        GeodesicKind::Null,
        8.0,
        kPi / 2.0,
        1.0,
        2.0,
        0.0,
        -1.0);
    const auto equatorial_separated =
        detail::initialize_kerr_separated_state(
            kerr,
            equatorial,
            GeodesicKind::Null,
            1.0e-12,
            1.0e-10,
            1.0);
    check(
        "equatorial polar motion is locked",
        equatorial_separated.state.polar_direction ==
            detail::SeparatedDirection::Locked);

    const KerrBoyerLindquistMetric schwarzschild(1.0, 0.0);
    const double outer_turn_radius = 4.4533631938113558;
    const PhaseSpaceState simple_root{
        0.0,
        Contravariant4{
            Vec4{{0.0, outer_turn_radius, kPi / 2.0, 0.0}}},
        Covariant4{
            Vec4{{-1.0, 0.0, 0.0, 6.0}}},
    };
    check(
        "simple-root fixture satisfies Hamiltonian constraint",
        hamiltonian_constraint_error(
            schwarzschild,
            simple_root,
            GeodesicKind::Null) <
            1.0e-14);
    const auto root_forward =
        detail::initialize_kerr_separated_state(
            schwarzschild,
            simple_root,
            GeodesicKind::Null,
            1.0e-12,
            1.0e-10,
            1.0);
    check(
        "positive integration leaves outer root outward",
        root_forward.state.radial_direction ==
            detail::SeparatedDirection::Positive);
    const auto root_backward =
        detail::initialize_kerr_separated_state(
            schwarzschild,
            simple_root,
            GeodesicKind::Null,
            1.0e-12,
            1.0e-10,
            -1.0);
    check(
        "negative integration leaves outer root on past branch",
        root_backward.state.radial_direction ==
            detail::SeparatedDirection::Negative);

    const PhaseSpaceState critical_root{
        0.0,
        Contravariant4{
            Vec4{{0.0, 3.0, kPi / 2.0, 0.0}}},
        Covariant4{
            Vec4{{-1.0, 0.0, 0.0, std::sqrt(27.0)}}},
    };
    check_throws<detail::KerrSeparatedCriticalInitialState>(
        "spherical photon double root is explicit",
        [&] {
            static_cast<void>(
                detail::initialize_kerr_separated_state(
                    schwarzschild,
                    critical_root,
                    GeodesicKind::Null,
                    1.0e-12,
                    1.0e-10,
                    1.0));
        });

    PhaseSpaceState inconsistent = initial;
    inconsistent.p.v[1] = 0.0;
    check_throws<std::domain_error>(
        "inconsistent Hamiltonian and radial potential rejected",
        [&] {
            static_cast<void>(
                detail::initialize_kerr_separated_state(
                    kerr,
                    inconsistent,
                    GeodesicKind::Null,
                    1.0e-12,
                    1.0e-10,
                    1.0));
        });

    detail::KerrSeparatedState axis_state =
        equatorial_separated.state;
    axis_state.values[detail::kMu] = 1.0;
    check_throws<std::domain_error>(
        "BL polar axis reconstruction rejected",
        [&] {
            static_cast<void>(
                detail::reconstruct_kerr_phase_space(
                    kerr,
                    equatorial_separated.constants,
                    axis_state));
        });

    check_throws<std::invalid_argument>(
        "zero integration direction rejected",
        [&] {
            static_cast<void>(
                detail::initialize_kerr_separated_state(
                    kerr,
                    initial,
                    GeodesicKind::Null,
                    1.0e-12,
                    1.0e-10,
                    0.0));
        });

    std::cout << "\n=== Results: " << passed
              << " passed, " << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
