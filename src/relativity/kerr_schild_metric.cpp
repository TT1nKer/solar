#include "solar/relativity/kerr_schild_metric.h"

#include "kerr_schild_fields.h"

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace solar::relativity {
namespace {

using DoubleFields = detail::KerrSchildFields<double>;

std::array<double, 4> components(
    const Contravariant4& x) {
    return {
        x.v[0],
        x.v[1],
        x.v[2],
        x.v[3],
    };
}

Mat4 to_matrix(
    const std::array<std::array<double, 4>, 4>& input) {
    Mat4 result{};
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            result[row][column] = input[row][column];
        }
    }
    return result;
}

bool all_fields_finite(
    const DoubleFields& fields) noexcept {
    if (!std::isfinite(fields.radius) ||
        !std::isfinite(fields.scalar_h)) {
        return false;
    }
    for (const double component : fields.null_covector) {
        if (!std::isfinite(component)) {
            return false;
        }
    }
    return all_finite(to_matrix(fields.covariant)) &&
           all_finite(to_matrix(fields.contravariant));
}

void require_valid(
    const KerrSchildCartesianMetric& metric,
    const Contravariant4& x) {
    if (!metric.valid_point(x)) {
        throw std::domain_error(
            "point is outside the positive-radius Cartesian "
            "Kerr-Schild domain");
    }
}

} // namespace

KerrSchildCartesianMetric::KerrSchildCartesianMetric(
    double mass_M,
    double spin_chi)
    : mass_M_(mass_M),
      spin_chi_(spin_chi),
      spin_a_M_(mass_M * spin_chi),
      outer_horizon_M_(0.0),
      inner_horizon_M_(0.0),
      minimum_radius_M_(
          64.0 * std::numeric_limits<double>::epsilon() *
          mass_M) {
    if (!std::isfinite(mass_M_) || mass_M_ <= 0.0) {
        throw std::invalid_argument(
            "Kerr mass must be finite and positive");
    }
    if (!std::isfinite(spin_chi_) ||
        std::fabs(spin_chi_) >= 1.0) {
        throw std::invalid_argument(
            "Kerr dimensionless spin must be finite with abs(chi)<1");
    }
    if (!std::isfinite(spin_a_M_) ||
        !std::isfinite(minimum_radius_M_) ||
        minimum_radius_M_ <= 0.0) {
        throw std::invalid_argument(
            "Kerr derived length scales must remain finite and positive");
    }

    const double horizon_factor =
        std::sqrt(1.0 - spin_chi_ * spin_chi_);
    outer_horizon_M_ =
        mass_M_ * (1.0 + horizon_factor);
    inner_horizon_M_ =
        mass_M_ * spin_chi_ * spin_chi_ /
        (1.0 + horizon_factor);
    if (!std::isfinite(outer_horizon_M_) ||
        !std::isfinite(inner_horizon_M_)) {
        throw std::invalid_argument(
            "Kerr horizon radii must remain finite");
    }
}

Mat4 KerrSchildCartesianMetric::covariant(
    const Contravariant4& x) const {
    require_valid(*this, x);
    return to_matrix(
        detail::evaluate_kerr_schild_fields(
            components(x), mass_M_, spin_a_M_)
            .covariant);
}

Mat4 KerrSchildCartesianMetric::contravariant(
    const Contravariant4& x) const {
    require_valid(*this, x);
    return to_matrix(
        detail::evaluate_kerr_schild_fields(
            components(x), mass_M_, spin_a_M_)
            .contravariant);
}

std::array<Mat4, 4>
KerrSchildCartesianMetric::contravariant_derivatives(
    const Contravariant4& x) const {
    require_valid(*this, x);
    std::array<Dual4, 4> dual_coordinates{};
    for (std::size_t coordinate = 0;
         coordinate < 4;
         ++coordinate) {
        dual_coordinates[coordinate] =
            Dual4::variable(x.v[coordinate], coordinate);
    }
    const auto fields =
        detail::evaluate_kerr_schild_fields(
            dual_coordinates, mass_M_, spin_a_M_);

    std::array<Mat4, 4> result{};
    for (std::size_t coordinate = 0;
         coordinate < 4;
         ++coordinate) {
        for (std::size_t row = 0; row < 4; ++row) {
            for (std::size_t column = 0;
                 column < 4;
                 ++column) {
                result[coordinate][row][column] =
                    fields.contravariant[row][column]
                        .derivative[coordinate];
            }
        }
        if (!all_finite(result[coordinate])) {
            throw std::domain_error(
                "Kerr-Schild inverse derivative is non-finite");
        }
    }
    return result;
}

bool KerrSchildCartesianMetric::valid_point(
    const Contravariant4& x) const noexcept {
    if (!x.v.all_finite()) {
        return false;
    }
    try {
        const DoubleFields fields =
            detail::evaluate_kerr_schild_fields(
                components(x), mass_M_, spin_a_M_);
        return fields.radius > minimum_radius_M_ &&
               all_fields_finite(fields);
    } catch (...) {
        return false;
    }
}

double KerrSchildCartesianMetric::radial_coordinate(
    const Contravariant4& x) const {
    require_valid(*this, x);
    return detail::evaluate_kerr_schild_fields(
        components(x), mass_M_, spin_a_M_)
        .radius;
}

Vec3 KerrSchildCartesianMetric::radial_coordinate_gradient(
    const Contravariant4& x) const {
    require_valid(*this, x);
    const double radius = radial_coordinate(x);
    const double radius_squared = radius * radius;
    const double spin_squared = spin_a_M_ * spin_a_M_;
    const double rho_squared =
        x.v[1] * x.v[1] +
        x.v[2] * x.v[2] +
        x.v[3] * x.v[3];
    const double denominator =
        2.0 * radius_squared - rho_squared + spin_squared;
    const double vertical_denominator =
        radius * denominator;
    if (!std::isfinite(denominator) ||
        denominator == 0.0 ||
        !std::isfinite(vertical_denominator) ||
        vertical_denominator == 0.0) {
        throw std::domain_error(
            "Kerr-Schild radial gradient is singular");
    }
    const Vec3 result{{
        x.v[1] * radius / denominator,
        x.v[2] * radius / denominator,
        x.v[3] * (radius_squared + spin_squared) /
            vertical_denominator,
    }};
    if (!result.all_finite()) {
        throw std::domain_error(
            "Kerr-Schild radial gradient is non-finite");
    }
    return result;
}

double kerr_schild_stationary_energy(
    const PhaseSpaceState& state) {
    return -state.p.v[0];
}

double kerr_schild_axial_angular_momentum(
    const PhaseSpaceState& state) {
    return state.x.v[1] * state.p.v[2] -
           state.x.v[2] * state.p.v[1];
}

} // namespace solar::relativity
