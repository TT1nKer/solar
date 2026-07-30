#include "solar/relativity/minkowski_metric.h"

#include <stdexcept>

namespace solar::relativity {
namespace {

Mat4 minkowski_matrix() {
    Mat4 metric{};
    metric[0][0] = -1.0;
    metric[1][1] = 1.0;
    metric[2][2] = 1.0;
    metric[3][3] = 1.0;
    return metric;
}

void require_valid(const MinkowskiMetric& metric, const Contravariant4& x) {
    if (!metric.valid_point(x)) {
        throw std::domain_error(
            "Minkowski coordinates must contain only finite values");
    }
}

} // namespace

const char* chart_name(Chart chart) noexcept {
    switch (chart) {
        case Chart::MinkowskiCartesian:
            return "minkowski-cartesian";
        case Chart::BoyerLindquist:
            return "boyer-lindquist";
        case Chart::KerrSchildCartesian:
            return "kerr-schild-cartesian";
    }
    return "unknown";
}

Mat4 MinkowskiMetric::covariant(const Contravariant4& x) const {
    require_valid(*this, x);
    return minkowski_matrix();
}

Mat4 MinkowskiMetric::contravariant(const Contravariant4& x) const {
    require_valid(*this, x);
    return minkowski_matrix();
}

std::array<Mat4, 4>
MinkowskiMetric::contravariant_derivatives(
    const Contravariant4& x) const {
    require_valid(*this, x);
    return {};
}

bool MinkowskiMetric::valid_point(const Contravariant4& x) const noexcept {
    return x.v.all_finite();
}

} // namespace solar::relativity
