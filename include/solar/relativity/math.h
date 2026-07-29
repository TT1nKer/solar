#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace solar::relativity {

template <std::size_t N>
struct Vec {
    std::array<double, N> data{};

    double& operator[](std::size_t index) { return data[index]; }
    const double& operator[](std::size_t index) const { return data[index]; }

    bool all_finite() const noexcept {
        for (const double component : data) {
            if (!std::isfinite(component)) {
                return false;
            }
        }
        return true;
    }
};

using Vec3 = Vec<3>;
using Vec4 = Vec<4>;
using Mat4 = std::array<std::array<double, 4>, 4>;

template <std::size_t N>
Vec<N> operator+(const Vec<N>& left, const Vec<N>& right) {
    Vec<N> result;
    for (std::size_t index = 0; index < N; ++index) {
        result[index] = left[index] + right[index];
    }
    return result;
}

template <std::size_t N>
Vec<N> operator-(const Vec<N>& left, const Vec<N>& right) {
    Vec<N> result;
    for (std::size_t index = 0; index < N; ++index) {
        result[index] = left[index] - right[index];
    }
    return result;
}

template <std::size_t N>
Vec<N> operator-(const Vec<N>& value) {
    Vec<N> result;
    for (std::size_t index = 0; index < N; ++index) {
        result[index] = -value[index];
    }
    return result;
}

template <std::size_t N>
Vec<N> operator*(const Vec<N>& value, double scalar) {
    Vec<N> result;
    for (std::size_t index = 0; index < N; ++index) {
        result[index] = value[index] * scalar;
    }
    return result;
}

template <std::size_t N>
Vec<N> operator*(double scalar, const Vec<N>& value) {
    return value * scalar;
}

template <std::size_t N>
Vec<N> operator/(const Vec<N>& value, double scalar) {
    if (!std::isfinite(scalar) || scalar == 0.0) {
        throw std::domain_error("vector divisor must be finite and nonzero");
    }
    return value * (1.0 / scalar);
}

template <std::size_t N>
double max_norm(const Vec<N>& value) noexcept {
    double maximum = 0.0;
    for (const double component : value.data) {
        maximum = std::max(maximum, std::fabs(component));
    }
    return maximum;
}

inline bool all_finite(const Mat4& matrix) noexcept {
    for (const auto& row : matrix) {
        for (const double component : row) {
            if (!std::isfinite(component)) {
                return false;
            }
        }
    }
    return true;
}

inline Vec4 multiply(const Mat4& matrix, const Vec4& vector) {
    Vec4 result;
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            result[row] += matrix[row][column] * vector[column];
        }
    }
    return result;
}

inline Mat4 multiply(const Mat4& left, const Mat4& right) {
    Mat4 result{};
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            for (std::size_t inner = 0; inner < 4; ++inner) {
                result[row][column] +=
                    left[row][inner] * right[inner][column];
            }
        }
    }
    return result;
}

/**
 * Inverts a finite 4x4 matrix using Gauss-Jordan elimination with partial
 * pivoting.
 *
 * Failure:
 *   Throws std::domain_error for a non-finite matrix or a pivot indistinguish-
 *   able from zero at double precision. No regularization is applied.
 *
 * Validation:
 *   tests/relativity/test_math.cpp
 */
inline Mat4 inverse(const Mat4& matrix) {
    if (!all_finite(matrix)) {
        throw std::domain_error("cannot invert a non-finite matrix");
    }

    std::array<std::array<double, 8>, 4> augmented{};
    double matrix_scale = 0.0;
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            augmented[row][column] = matrix[row][column];
            matrix_scale =
                std::max(matrix_scale, std::fabs(matrix[row][column]));
        }
        augmented[row][row + 4] = 1.0;
    }

    const double pivot_floor =
        16.0 * std::numeric_limits<double>::epsilon() * matrix_scale;
    if (matrix_scale == 0.0) {
        throw std::domain_error("cannot invert a singular matrix");
    }

    for (std::size_t column = 0; column < 4; ++column) {
        std::size_t pivot_row = column;
        for (std::size_t candidate = column + 1; candidate < 4; ++candidate) {
            if (std::fabs(augmented[candidate][column]) >
                std::fabs(augmented[pivot_row][column])) {
                pivot_row = candidate;
            }
        }

        const double pivot = augmented[pivot_row][column];
        if (!std::isfinite(pivot) || std::fabs(pivot) <= pivot_floor) {
            throw std::domain_error(
                "cannot invert a singular or ill-conditioned matrix");
        }
        if (pivot_row != column) {
            std::swap(augmented[pivot_row], augmented[column]);
        }

        const double selected_pivot = augmented[column][column];
        for (double& component : augmented[column]) {
            component /= selected_pivot;
        }

        for (std::size_t row = 0; row < 4; ++row) {
            if (row == column) {
                continue;
            }
            const double factor = augmented[row][column];
            for (std::size_t entry = 0; entry < 8; ++entry) {
                augmented[row][entry] -=
                    factor * augmented[column][entry];
            }
        }
    }

    Mat4 result{};
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            result[row][column] = augmented[row][column + 4];
        }
    }
    if (!all_finite(result)) {
        throw std::domain_error("matrix inversion produced a non-finite result");
    }
    return result;
}

/**
 * Contracts two contravariant component arrays with diag(-1,+1,+1,+1).
 *
 * This explicitly named operation is the only metric-free spacetime
 * contraction in Phase 0A. General contractions require a supplied metric.
 */
inline double minkowski_dot_minus_plus_plus_plus(
    const Vec4& left, const Vec4& right) noexcept {
    return -left[0] * right[0] + left[1] * right[1] +
           left[2] * right[2] + left[3] * right[3];
}

} // namespace solar::relativity
