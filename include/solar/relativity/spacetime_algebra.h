#pragma once

#include "solar/relativity/types.h"

#include <cstddef>

namespace solar::relativity {

inline double covector_vector_pairing(
    const Covariant4& covector,
    const Contravariant4& vector) noexcept {
    double result = 0.0;
    for (std::size_t component = 0; component < 4; ++component) {
        result += covector.v[component] * vector.v[component];
    }
    return result;
}

inline Covariant4 lower_index(
    const Mat4& covariant,
    const Contravariant4& vector) noexcept {
    Covariant4 result;
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            result.v[row] +=
                covariant[row][column] * vector.v[column];
        }
    }
    return result;
}

inline Contravariant4 raise_index(
    const Mat4& contravariant,
    const Covariant4& covector) noexcept {
    Contravariant4 result;
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            result.v[row] +=
                contravariant[row][column] * covector.v[column];
        }
    }
    return result;
}

inline double metric_inner_product(
    const Mat4& covariant,
    const Contravariant4& left,
    const Contravariant4& right) noexcept {
    return covector_vector_pairing(
        lower_index(covariant, left), right);
}

} // namespace solar::relativity
