#include "thin_disk_geometry.h"

#include "solar/relativity/kerr_bl_metric.h"
#include "solar/relativity/kerr_chart_transform.h"
#include "solar/relativity/kerr_schild_metric.h"
#include "solar/relativity/spacetime_algebra.h"

#include <cmath>
#include <stdexcept>

namespace solar::relativity::detail {
namespace {

constexpr double surface_geometry_tolerance = 1.0e-10;

Contravariant4 north_normal_in_boyer_lindquist(
    const KerrBoyerLindquistMetric& metric,
    const Contravariant4& position) {
    const Mat4 covariant = metric.covariant(position);
    if (!std::isfinite(covariant[2][2]) ||
        covariant[2][2] <= 0.0) {
        throw std::domain_error(
            "thin-disk BL normal is not defined");
    }
    Contravariant4 normal;
    normal.v[2] = -1.0 / std::sqrt(covariant[2][2]);
    return normal;
}

} // namespace

ThinDiskSurfaceGeometry evaluate_thin_disk_surface_geometry(
    const Metric& metric,
    const Contravariant4& position) {
    if (const auto* kerr_bl =
            dynamic_cast<
                const KerrBoyerLindquistMetric*>(&metric)) {
        return ThinDiskSurfaceGeometry{
            position.v[1],
            north_normal_in_boyer_lindquist(
                *kerr_bl, position),
        };
    }

    const auto* kerr_schild =
        dynamic_cast<
            const KerrSchildCartesianMetric*>(&metric);
    if (kerr_schild == nullptr) {
        throw std::domain_error(
            "thin-disk normal requires a supported Kerr chart");
    }
    const KerrChartTransform transform(
        kerr_schild->mass(),
        kerr_schild->spin_chi());
    const Contravariant4 bl_position =
        transform.position_to_boyer_lindquist(position);
    const KerrBoyerLindquistMetric bl_metric(
        kerr_schild->mass(),
        kerr_schild->spin_chi());
    const Contravariant4 bl_normal =
        north_normal_in_boyer_lindquist(
            bl_metric, bl_position);
    return ThinDiskSurfaceGeometry{
        bl_position.v[1],
        Contravariant4{multiply(
            transform
                .boyer_lindquist_to_kerr_schild_jacobian(
                    bl_position),
            bl_normal.v)},
    };
}

bool valid_thin_disk_surface_geometry(
    const Metric& metric,
    const Contravariant4& position,
    const Contravariant4& normal,
    const Contravariant4& emitter_four_velocity) {
    if (!normal.v.all_finite()) {
        return false;
    }
    const Mat4 covariant = metric.covariant(position);
    const double normal_norm = metric_inner_product(
        covariant, normal, normal);
    const double orthogonality = metric_inner_product(
        covariant, normal, emitter_four_velocity);
    return std::isfinite(normal_norm) &&
           std::isfinite(orthogonality) &&
           std::fabs(normal_norm - 1.0) <=
               surface_geometry_tolerance &&
           std::fabs(orthogonality) <=
               surface_geometry_tolerance;
}

} // namespace solar::relativity::detail
