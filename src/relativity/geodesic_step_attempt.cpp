#include "geodesic_step_attempt.h"

#include <limits>
#include <stdexcept>
#include <utility>

namespace solar::relativity::detail {

GeodesicStepAttempt attempt_geodesic_step(
    const Metric& metric,
    const PhaseSpaceState& current,
    double attempted_step,
    const numerics::Dopri5Config<8>& config) {
    bool invalid_metric_point = false;
    bool non_finite_stage = false;
    const auto packed_rhs =
        [&metric, &invalid_metric_point, &non_finite_stage](
            double affine,
            const numerics::StateN<8>& packed) {
            numerics::StateN<8> derivative{};
            PhaseSpaceState state;
            try {
                state = unpack_phase_space(affine, packed);
            } catch (const std::domain_error&) {
                non_finite_stage = true;
                derivative.fill(
                    std::numeric_limits<double>::quiet_NaN());
                return derivative;
            }
            if (!metric.valid_point(state.x)) {
                invalid_metric_point = true;
                derivative.fill(
                    std::numeric_limits<double>::quiet_NaN());
                return derivative;
            }
            try {
                const PhaseSpaceDerivative rhs =
                    HamiltonGeodesicRhs(metric)(state);
                for (std::size_t component = 0;
                     component < 4;
                     ++component) {
                    derivative[component] =
                        rhs.dx.v[component];
                    derivative[component + 4] =
                        rhs.dp.v[component];
                }
            } catch (const std::domain_error&) {
                non_finite_stage = true;
                derivative.fill(
                    std::numeric_limits<double>::quiet_NaN());
            }
            return derivative;
        };

    auto step = numerics::dopri5_step(
        pack_phase_space(current),
        current.affine,
        attempted_step,
        packed_rhs,
        config);
    return GeodesicStepAttempt{
        std::move(step),
        invalid_metric_point,
        non_finite_stage,
    };
}

} // namespace solar::relativity::detail
