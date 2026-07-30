#include "solar/relativity/kerr_schild_metric.h"

#include <algorithm>
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

void check_near(
    const char* name,
    double actual,
    double expected,
    double tolerance) {
    check(name, std::isfinite(actual) &&
                    std::fabs(actual - expected) <= tolerance);
}

Contravariant4 from_oblate(
    double time,
    double radius,
    double theta,
    double azimuth,
    double spin_a) {
    const double sine = std::sin(theta);
    const double cosine = std::cos(theta);
    return Contravariant4{Vec4{{
        time,
        (radius * std::cos(azimuth) -
         spin_a * std::sin(azimuth)) * sine,
        (radius * std::sin(azimuth) +
         spin_a * std::cos(azimuth)) * sine,
        radius * cosine,
    }}};
}

double inverse_identity_error(
    const Mat4& covariant,
    const Mat4& contravariant) {
    const Mat4 product = multiply(covariant, contravariant);
    double maximum = 0.0;
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            const double expected = row == column ? 1.0 : 0.0;
            maximum = std::max(
                maximum,
                std::fabs(product[row][column] - expected));
        }
    }
    return maximum;
}

void check_metric_point(
    const char* validity_name,
    const KerrSchildCartesianMetric& metric,
    const Contravariant4& point,
    double inverse_gate) {
    check(validity_name, metric.valid_point(point));
    const double radius = metric.radial_coordinate(point);
    const double a = metric.spin_length();
    const double x = point.v[1];
    const double y = point.v[2];
    const double z = point.v[3];
    const double rho_squared = x * x + y * y + z * z;
    const double residual =
        radius * radius * radius * radius -
        (rho_squared - a * a) * radius * radius -
        a * a * z * z;
    const double residual_scale = std::max(
        {1.0,
         radius * radius * radius * radius,
         std::fabs((rho_squared - a * a) *
                   radius * radius),
         std::fabs(a * a * z * z)});
    check("implicit radius quartic residual",
          std::fabs(residual) <
              2.0e-14 * residual_scale);

    const double r2 = radius * radius;
    const double l0 = 1.0;
    const double lx = (radius * x + a * y) / (r2 + a * a);
    const double ly = (radius * y - a * x) / (r2 + a * a);
    const double lz = z / radius;
    check("KS one-form is Minkowski null",
          std::fabs(-l0 * l0 + lx * lx + ly * ly + lz * lz) <
              5.0e-13);

    const Mat4 covariant = metric.covariant(point);
    const Mat4 contravariant = metric.contravariant(point);
    check("KS covariant metric finite", all_finite(covariant));
    check("KS inverse metric finite", all_finite(contravariant));
    bool covariant_symmetric = true;
    bool contravariant_symmetric = true;
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            covariant_symmetric &=
                covariant[row][column] ==
                covariant[column][row];
            contravariant_symmetric &=
                contravariant[row][column] ==
                contravariant[column][row];
        }
    }
    check("KS covariant metric symmetric",
          covariant_symmetric);
    check("KS inverse metric symmetric",
          contravariant_symmetric);
    check("KS analytic inverse identity",
          inverse_identity_error(covariant, contravariant) <
              inverse_gate);
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
    const KerrSchildCartesianMetric metric(2.0, 0.6);
    check("KS chart",
          metric.chart() == Chart::KerrSchildCartesian);
    check("KS name", metric.name() == "kerr-schild-cartesian");
    check_near("KS mass", metric.mass(), 2.0, 0.0);
    check_near("KS dimensionless spin",
               metric.spin_chi(), 0.6, 0.0);
    check_near("KS spin length",
               metric.spin_length(), 1.2, 0.0);
    check_near("KS outer horizon",
               metric.outer_horizon_radius(), 3.6, 1.0e-15);
    check_near("KS inner horizon",
               metric.inner_horizon_radius(), 0.4, 2.0e-16);

    const Contravariant4 exterior =
        from_oblate(3.0, 8.0, 1.1, 0.4, 1.2);
    const Contravariant4 axis =
        from_oblate(-2.0, 5.0, 0.0, 0.0, 1.2);
    const Contravariant4 horizon =
        from_oblate(0.0, 3.6, 1.2, -0.3, 1.2);
    const Contravariant4 inside =
        from_oblate(1.0, 1.5, 0.9, 0.8, 1.2);
    check_metric_point(
        "ordinary exterior accepted",
        metric, exterior, 5.0e-13);
    check_metric_point(
        "polar axis accepted",
        metric, axis, 5.0e-13);
    check_metric_point(
        "outer horizon accepted",
        metric, horizon, 1.0e-10);
    check_metric_point(
        "positive-radius interior accepted",
        metric, inside, 5.0e-13);

    const KerrSchildCartesianMetric schwarzschild(1.0, 0.0);
    const Contravariant4 cartesian{
        Vec4{{4.0, 3.0, 4.0, 0.0}}};
    const Mat4 schwarzschild_covariant =
        schwarzschild.covariant(cartesian);
    const double spatial_l[3] = {0.6, 0.8, 0.0};
    double schwarzschild_max_error = 0.0;
    for (std::size_t row = 0; row < 4; ++row) {
        const double l_row =
            row == 0 ? 1.0 : spatial_l[row - 1];
        for (std::size_t column = 0; column < 4; ++column) {
            const double l_column =
                column == 0 ? 1.0 : spatial_l[column - 1];
            const double eta =
                row == column ? (row == 0 ? -1.0 : 1.0) : 0.0;
            schwarzschild_max_error = std::max(
                schwarzschild_max_error,
                std::fabs(
                    schwarzschild_covariant[row][column] -
                    (eta + 0.4 * l_row * l_column)));
        }
    }
    check("Schwarzschild Cartesian KS metric",
          schwarzschild_max_error < 2.0e-16);

    const Vec3 gradient =
        schwarzschild.radial_coordinate_gradient(cartesian);
    check_near("Schwarzschild dr/dx", gradient[0], 0.6, 2.0e-16);
    check_near("Schwarzschild dr/dy", gradient[1], 0.8, 2.0e-16);
    check_near("Schwarzschild dr/dz", gradient[2], 0.0, 0.0);

    const auto derivatives =
        schwarzschild.contravariant_derivatives(cartesian);
    bool derivatives_finite = true;
    bool derivatives_symmetric = true;
    for (const auto& derivative : derivatives) {
        derivatives_finite &= all_finite(derivative);
        for (std::size_t row = 0; row < 4; ++row) {
            for (std::size_t column = 0; column < 4; ++column) {
                derivatives_symmetric &=
                    derivative[row][column] ==
                    derivative[column][row];
            }
        }
    }
    check("KS inverse derivatives finite",
          derivatives_finite);
    check("KS inverse derivatives symmetric",
          derivatives_symmetric);
    Mat4 zero{};
    check("stationary derivative is exactly zero",
          derivatives[0] == zero);

    const PhaseSpaceState state{
        0.0,
        exterior,
        Covariant4{Vec4{{-2.5, 0.3, -0.7, 0.4}}},
    };
    check_near("Cartesian stationary energy",
               kerr_schild_stationary_energy(state),
               2.5, 0.0);
    check_near("Cartesian axial angular momentum",
               kerr_schild_axial_angular_momentum(state),
               exterior.v[1] * -0.7 -
                   exterior.v[2] * 0.3,
               0.0);

    const Contravariant4 ring{
        Vec4{{0.0, metric.spin_length(), 0.0, 0.0}}};
    const Contravariant4 disk{
        Vec4{{0.0, 0.5 * metric.spin_length(), 0.0, 0.0}}};
    Contravariant4 non_finite = exterior;
    non_finite.v[2] =
        std::numeric_limits<double>::quiet_NaN();
    check("Kerr ring rejected", !metric.valid_point(ring));
    check("zero-radius oblate disk rejected",
          !metric.valid_point(disk));
    check("non-finite point rejected",
          !metric.valid_point(non_finite));
    check("Schwarzschild origin rejected",
          !schwarzschild.valid_point(Contravariant4{}));
    try {
        (void)metric.covariant(ring);
        check("invalid metric evaluation throws", false);
    } catch (const std::domain_error&) {
        check("invalid metric evaluation throws", true);
    }

    check_invalid_argument(
        "zero mass rejected",
        [] { (void)KerrSchildCartesianMetric(0.0, 0.0); });
    check_invalid_argument(
        "non-finite mass rejected",
        [] {
            (void)KerrSchildCartesianMetric(
                std::numeric_limits<double>::infinity(), 0.0);
        });
    check_invalid_argument(
        "extremal spin rejected",
        [] { (void)KerrSchildCartesianMetric(1.0, 1.0); });
    check_invalid_argument(
        "negative extremal spin rejected",
        [] { (void)KerrSchildCartesianMetric(1.0, -1.0); });
    check_invalid_argument(
        "non-finite spin rejected",
        [] {
            (void)KerrSchildCartesianMetric(
                1.0,
                std::numeric_limits<double>::quiet_NaN());
        });
    check_invalid_argument(
        "overflowing horizon scale rejected",
        [] {
            (void)KerrSchildCartesianMetric(
                std::numeric_limits<double>::max(), 0.5);
        });

    std::cout << "\n=== Results: " << passed << " passed, "
              << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
