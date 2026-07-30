#pragma once

#include "solar/relativity/dual4.h"

#include <array>
#include <cmath>
#include <cstddef>

namespace solar::relativity::detail {

template <typename Scalar>
struct KerrSchildFields {
    Scalar radius;
    Scalar scalar_h;
    std::array<Scalar, 4> null_covector;
    std::array<std::array<Scalar, 4>, 4> covariant;
    std::array<std::array<Scalar, 4>, 4> contravariant;
};

inline double scalar_value(double value) noexcept {
    return value;
}

inline double scalar_value(const Dual4& value) noexcept {
    return value.value;
}

template <typename Scalar>
KerrSchildFields<Scalar> evaluate_kerr_schild_fields(
    const std::array<Scalar, 4>& coordinates,
    double mass,
    double spin_a) {
    using std::sqrt;

    const Scalar& x = coordinates[1];
    const Scalar& y = coordinates[2];
    const Scalar& z = coordinates[3];
    const double spin_squared = spin_a * spin_a;
    const Scalar rho_squared = x * x + y * y + z * z;
    const Scalar q = rho_squared - spin_squared;
    const Scalar discriminant = sqrt(
        q * q + 4.0 * spin_squared * z * z);

    const Scalar radius_squared =
        scalar_value(q) >= 0.0
            ? 0.5 * (q + discriminant)
            : 2.0 * spin_squared * z * z /
                  (discriminant - q);
    const Scalar radius = sqrt(radius_squared);
    const Scalar z_over_radius = z / radius;
    const Scalar scalar_h =
        mass * radius /
        (radius_squared +
         spin_squared * z_over_radius * z_over_radius);
    const Scalar radial_spin_sum =
        radius_squared + spin_squared;

    std::array<Scalar, 4> null_covector{};
    null_covector[0] = Scalar{1.0};
    null_covector[1] =
        (radius * x + spin_a * y) / radial_spin_sum;
    null_covector[2] =
        (radius * y - spin_a * x) / radial_spin_sum;
    null_covector[3] = z_over_radius;

    std::array<Scalar, 4> null_contravariant =
        null_covector;
    null_contravariant[0] = -null_contravariant[0];

    KerrSchildFields<Scalar> result{};
    result.radius = radius;
    result.scalar_h = scalar_h;
    result.null_covector = null_covector;
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = row;
             column < 4;
             ++column) {
            const double eta =
                row == column
                    ? (row == 0 ? -1.0 : 1.0)
                    : 0.0;
            const Scalar covariant =
                Scalar{eta} +
                2.0 * scalar_h *
                    null_covector[row] *
                    null_covector[column];
            const Scalar contravariant =
                Scalar{eta} -
                2.0 * scalar_h *
                    null_contravariant[row] *
                    null_contravariant[column];
            result.covariant[row][column] = covariant;
            result.covariant[column][row] = covariant;
            result.contravariant[row][column] =
                contravariant;
            result.contravariant[column][row] =
                contravariant;
        }
    }
    return result;
}

} // namespace solar::relativity::detail
