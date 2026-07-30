#include "solar/relativity/kerr_chart_transform.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace solar::relativity {
namespace {

void require_finite_state(
    const PhaseSpaceState& state,
    const char* chart_name) {
    if (!std::isfinite(state.affine) ||
        !state.x.v.all_finite() ||
        !state.p.v.all_finite()) {
        throw std::domain_error(
            std::string(chart_name) +
            " transform state must be finite");
    }
}

Covariant4 boyer_lindquist_covector_to_kerr_schild(
    const Mat4& inverse_jacobian,
    const Covariant4& source) {
    Covariant4 result{};
    for (std::size_t target = 0; target < 4; ++target) {
        for (std::size_t source_component = 0;
             source_component < 4;
             ++source_component) {
            result.v[target] +=
                inverse_jacobian[source_component][target] *
                source.v[source_component];
        }
    }
    if (!result.v.all_finite()) {
        throw std::domain_error(
            "Kerr chart covector transform is non-finite");
    }
    return result;
}

Covariant4 kerr_schild_covector_to_boyer_lindquist(
    const Mat4& forward_jacobian,
    const Covariant4& source) {
    Covariant4 result{};
    for (std::size_t target = 0; target < 4; ++target) {
        for (std::size_t source_component = 0;
             source_component < 4;
             ++source_component) {
            result.v[target] +=
                forward_jacobian[source_component][target] *
                source.v[source_component];
        }
    }
    if (!result.v.all_finite()) {
        throw std::domain_error(
            "inverse Kerr chart covector transform is non-finite");
    }
    return result;
}

} // namespace

PhaseSpaceState KerrChartTransform::state_to_kerr_schild(
    const PhaseSpaceState& boyer_lindquist) const {
    require_finite_state(
        boyer_lindquist,
        "Boyer-Lindquist");
    const Contravariant4 position =
        position_to_kerr_schild(boyer_lindquist.x);
    const Mat4 forward =
        boyer_lindquist_to_kerr_schild_jacobian(
            boyer_lindquist.x);
    const Mat4 inverse_jacobian = inverse(forward);
    return PhaseSpaceState{
        boyer_lindquist.affine,
        position,
        boyer_lindquist_covector_to_kerr_schild(
            inverse_jacobian,
            boyer_lindquist.p),
    };
}

PhaseSpaceState KerrChartTransform::state_to_boyer_lindquist(
    const PhaseSpaceState& kerr_schild) const {
    require_finite_state(kerr_schild, "Kerr-Schild");
    const Contravariant4 position =
        position_to_boyer_lindquist(kerr_schild.x);
    const Mat4 forward =
        boyer_lindquist_to_kerr_schild_jacobian(position);
    return PhaseSpaceState{
        kerr_schild.affine,
        position,
        kerr_schild_covector_to_boyer_lindquist(
            forward,
            kerr_schild.p),
    };
}

} // namespace solar::relativity
