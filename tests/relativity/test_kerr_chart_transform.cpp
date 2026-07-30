#include "solar/constants.h"
#include "solar/relativity/hamiltonian.h"
#include "solar/relativity/kerr_bl_metric.h"
#include "solar/relativity/kerr_chart_transform.h"
#include "solar/relativity/kerr_schild_metric.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace solar::relativity;

namespace {

int passed = 0;
int failed = 0;

void check(
    const char* name,
    bool condition,
    double diagnostic = 0.0) {
    std::cout << (condition ? "  PASS: " : "  FAIL: ")
              << name << " (error=" << diagnostic << ")\n";
    condition ? ++passed : ++failed;
}

double wrapped_difference(
    double left,
    double right) {
    return std::remainder(
        left - right,
        2.0 * solar::constants::PI);
}

double matrix_identity_error(const Mat4& matrix) {
    double maximum = 0.0;
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            const double expected = row == column ? 1.0 : 0.0;
            maximum = std::max(
                maximum,
                std::fabs(matrix[row][column] - expected));
        }
    }
    return maximum;
}

double position_error(
    const Contravariant4& actual,
    const Contravariant4& expected) {
    return std::max(
        {std::fabs(actual.v[0] - expected.v[0]),
         std::fabs(actual.v[1] - expected.v[1]),
         std::fabs(actual.v[2] - expected.v[2]),
         std::fabs(wrapped_difference(
             actual.v[3], expected.v[3]))});
}

double momentum_error(
    const Covariant4& actual,
    const Covariant4& expected) {
    double maximum = 0.0;
    for (std::size_t component = 0; component < 4; ++component) {
        maximum = std::max(
            maximum,
            std::fabs(
                actual.v[component] -
                expected.v[component]));
    }
    return maximum;
}

double pairing(
    const Covariant4& covector,
    const Vec4& vector) {
    double result = 0.0;
    for (std::size_t component = 0; component < 4; ++component) {
        result += covector.v[component] * vector[component];
    }
    return result;
}

double recovered_tilde_phi(
    const Contravariant4& ks,
    double radius,
    double spin_a) {
    return std::atan2(
        radius * ks.v[2] - spin_a * ks.v[1],
        radius * ks.v[1] + spin_a * ks.v[2]);
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

struct TransformFixture {
    double mass;
    double spin_chi;
    Contravariant4 bl_position;
};

} // namespace

int main() {
    const std::vector<TransformFixture> fixtures{
        {1.0, 0.7,
         Contravariant4{
             Vec4{{2.0, 8.0, 1.1, 0.4}}}},
        {3.0, -0.8,
         Contravariant4{
             Vec4{{-5.0, 18.0, 2.0, -2.2}}}},
        {2.0, 0.0,
         Contravariant4{
             Vec4{{1.0, 12.0, 0.8, 2.6}}}},
    };

    double maximum_position_error = 0.0;
    double maximum_jacobian_error = 0.0;
    double maximum_momentum_error = 0.0;
    double maximum_pairing_error = 0.0;
    double maximum_hamiltonian_error = 0.0;

    for (const TransformFixture& fixture : fixtures) {
        const KerrChartTransform transform(
            fixture.mass, fixture.spin_chi);
        const KerrBoyerLindquistMetric bl_metric(
            fixture.mass, fixture.spin_chi);
        const KerrSchildCartesianMetric ks_metric(
            fixture.mass, fixture.spin_chi);

        const Contravariant4 ks_position =
            transform.position_to_kerr_schild(
                fixture.bl_position);
        const Contravariant4 bl_round_trip =
            transform.position_to_boyer_lindquist(
                ks_position);
        maximum_position_error = std::max(
            maximum_position_error,
            position_error(
                bl_round_trip,
                fixture.bl_position));

        const Mat4 forward =
            transform.boyer_lindquist_to_kerr_schild_jacobian(
                fixture.bl_position);
        const Mat4 reverse =
            transform.kerr_schild_to_boyer_lindquist_jacobian(
                ks_position);
        maximum_jacobian_error = std::max(
            {maximum_jacobian_error,
             matrix_identity_error(
                 multiply(reverse, forward)),
             matrix_identity_error(
                 multiply(forward, reverse))});

        const PhaseSpaceState bl_state{
            4.5,
            fixture.bl_position,
            Covariant4{
                Vec4{{-1.2, 0.3, -0.2, 2.0}}},
        };
        const PhaseSpaceState ks_state =
            transform.state_to_kerr_schild(bl_state);
        const PhaseSpaceState state_round_trip =
            transform.state_to_boyer_lindquist(ks_state);
        maximum_position_error = std::max(
            maximum_position_error,
            position_error(
                state_round_trip.x,
                bl_state.x));
        maximum_momentum_error = std::max(
            maximum_momentum_error,
            momentum_error(
                state_round_trip.p,
                bl_state.p));
        check("state transform preserves affine",
              ks_state.affine == bl_state.affine &&
                  state_round_trip.affine == bl_state.affine);

        const Vec4 bl_vector{{1.1, -0.4, 0.2, 0.7}};
        const Vec4 ks_vector = multiply(forward, bl_vector);
        maximum_pairing_error = std::max(
            maximum_pairing_error,
            std::fabs(
                pairing(bl_state.p, bl_vector) -
                pairing(ks_state.p, ks_vector)));
        maximum_hamiltonian_error = std::max(
            maximum_hamiltonian_error,
            std::fabs(
                hamiltonian(bl_metric, bl_state) -
                hamiltonian(ks_metric, ks_state)));
    }

    check("BL/KS position round trip",
          maximum_position_error < 1.0e-10,
          maximum_position_error);
    check("BL/KS Jacobians are mutual inverses",
          maximum_jacobian_error < 5.0e-12,
          maximum_jacobian_error);
    check("BL/KS covariant momentum round trip",
          maximum_momentum_error < 1.0e-10,
          maximum_momentum_error);
    check("covector-vector pairing is invariant",
          maximum_pairing_error < 1.0e-10,
          maximum_pairing_error);
    check("Hamiltonian is invariant",
          maximum_hamiltonian_error < 1.0e-10,
          maximum_hamiltonian_error);

    const double mass = 1.0;
    const double spin_chi = 0.6;
    const double spin_a = mass * spin_chi;
    const KerrChartTransform sign_transform(mass, spin_chi);
    const KerrSchildCartesianMetric sign_metric(
        mass, spin_chi);
    const double radius = 6.0;
    const double step = 1.0e-5;
    const Contravariant4 bl_minus{
        Vec4{{0.0, radius - step, 1.2, 0.3}}};
    const Contravariant4 bl_plus{
        Vec4{{0.0, radius + step, 1.2, 0.3}}};
    const Contravariant4 ks_minus =
        sign_transform.position_to_kerr_schild(bl_minus);
    const Contravariant4 ks_plus =
        sign_transform.position_to_kerr_schild(bl_plus);
    const double d_t_dr =
        (ks_plus.v[0] - ks_minus.v[0]) /
        (2.0 * step);
    const double d_phi_dr =
        wrapped_difference(
            recovered_tilde_phi(
                ks_plus, radius + step, spin_a),
            recovered_tilde_phi(
                ks_minus, radius - step, spin_a)) /
        (2.0 * step);
    const double delta =
        radius * radius - 2.0 * mass * radius +
        spin_a * spin_a;
    check("ingoing KS time differential has positive sign",
          std::fabs(d_t_dr - 2.0 * mass * radius / delta) <
              2.0e-10,
          std::fabs(
              d_t_dr - 2.0 * mass * radius / delta));
    check("ingoing KS azimuth differential has positive sign",
          std::fabs(d_phi_dr - spin_a / delta) <
              2.0e-10,
          std::fabs(d_phi_dr - spin_a / delta));

    const double unsafe_radius =
        sign_metric.outer_horizon_radius() + 1.0e-4;
    const Contravariant4 unsafe_bl{
        Vec4{{0.0, unsafe_radius, 1.0, 0.2}}};
    check_domain_error(
        "overlap boundary rejected",
        [&] {
            (void)sign_transform.position_to_kerr_schild(
                unsafe_bl);
        });
    const Contravariant4 unsafe_ks =
        KerrChartTransform(mass, spin_chi, 1.0e-8)
            .position_to_kerr_schild(unsafe_bl);
    check_domain_error(
        "inverse overlap boundary rejected",
        [&] {
            (void)sign_transform.position_to_boyer_lindquist(
                unsafe_ks);
        });

    Contravariant4 axis_bl{
        Vec4{{0.0, 8.0, 0.0, 0.0}}};
    check_domain_error(
        "BL polar axis rejected",
        [&] {
            (void)sign_transform.position_to_kerr_schild(
                axis_bl);
        });
    const Contravariant4 axis_ks{
        Vec4{{0.0, 0.0, 0.0, 8.0}}};
    check_domain_error(
        "KS polar axis rejected",
        [&] {
            (void)sign_transform.position_to_boyer_lindquist(
                axis_ks);
        });
    Contravariant4 non_finite_bl = fixtures[0].bl_position;
    non_finite_bl.v[3] =
        std::numeric_limits<double>::quiet_NaN();
    check_domain_error(
        "non-finite BL point rejected",
        [&] {
            (void)sign_transform.position_to_kerr_schild(
                non_finite_bl);
        });
    PhaseSpaceState non_finite_state{
        0.0,
        fixtures[0].bl_position,
        Covariant4{
            Vec4{{-1.0, 0.0, 0.0,
                  std::numeric_limits<double>::infinity()}}},
    };
    check_domain_error(
        "non-finite canonical momentum rejected",
        [&] {
            (void)sign_transform.state_to_kerr_schild(
                non_finite_state);
        });

    check_invalid_argument(
        "transform zero mass rejected",
        [] { (void)KerrChartTransform(0.0, 0.0); });
    check_invalid_argument(
        "transform extremal spin rejected",
        [] { (void)KerrChartTransform(1.0, 1.0); });
    check_invalid_argument(
        "transform non-finite spin rejected",
        [] {
            (void)KerrChartTransform(
                1.0,
                std::numeric_limits<double>::quiet_NaN());
        });
    check_invalid_argument(
        "transform zero overlap margin rejected",
        [] { (void)KerrChartTransform(1.0, 0.0, 0.0); });

    std::cout << "\n=== Results: " << passed << " passed, "
              << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
