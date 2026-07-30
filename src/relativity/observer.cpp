#include "solar/relativity/observer.h"

#include "solar/relativity/spacetime_algebra.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace solar::relativity {
namespace {

bool finite_tetrad(const Tetrad& tetrad) {
    for (const Contravariant4& leg : tetrad.basis) {
        if (!leg.v.all_finite()) {
            return false;
        }
    }
    return true;
}

} // namespace

double tetrad_orthonormality_error(
    const Metric& metric,
    const ObserverFrame& observer) {
    if (!observer.x.v.all_finite() ||
        !finite_tetrad(observer.tetrad) ||
        !metric.valid_point(observer.x)) {
        return std::numeric_limits<double>::infinity();
    }

    Mat4 covariant;
    try {
        covariant = metric.covariant(observer.x);
    } catch (const std::domain_error&) {
        return std::numeric_limits<double>::infinity();
    }

    double maximum = 0.0;
    for (std::size_t left = 0; left < 4; ++left) {
        for (std::size_t right = 0; right < 4; ++right) {
            const double expected =
                left == right ? (left == 0 ? -1.0 : 1.0)
                              : 0.0;
            const double actual = metric_inner_product(
                covariant,
                observer.tetrad.basis[left],
                observer.tetrad.basis[right]);
            if (!std::isfinite(actual)) {
                return std::numeric_limits<double>::infinity();
            }
            maximum = std::max(
                maximum, std::fabs(actual - expected));
        }
    }
    return maximum;
}

Contravariant4 tetrad_to_coordinate(
    const Tetrad& tetrad,
    const Vec4& local_components) noexcept {
    Contravariant4 coordinate;
    for (std::size_t leg = 0; leg < 4; ++leg) {
        for (std::size_t component = 0;
             component < 4;
             ++component) {
            coordinate.v[component] +=
                local_components[leg] *
                tetrad.basis[leg].v[component];
        }
    }
    return coordinate;
}

Vec4 coordinate_to_tetrad(
    const Mat4& covariant,
    const Tetrad& tetrad,
    const Contravariant4& coordinate_vector) noexcept {
    Vec4 local;
    local[0] = -metric_inner_product(
        covariant, tetrad.basis[0], coordinate_vector);
    for (std::size_t leg = 1; leg < 4; ++leg) {
        local[leg] = metric_inner_product(
            covariant, tetrad.basis[leg], coordinate_vector);
    }
    return local;
}

} // namespace solar::relativity
