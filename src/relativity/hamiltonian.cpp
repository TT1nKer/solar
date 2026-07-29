#include "solar/relativity/hamiltonian.h"

#include <cmath>
#include <stdexcept>

namespace solar::relativity {
namespace {

void require_finite_state(const PhaseSpaceState& state) {
    if (!std::isfinite(state.affine) ||
        !state.x.v.all_finite() ||
        !state.p.v.all_finite()) {
        throw std::domain_error(
            "Hamiltonian phase-space state must be finite");
    }
}

Mat4 validated_inverse(
    const Metric& metric,
    const PhaseSpaceState& state) {
    require_finite_state(state);
    if (!metric.valid_point(state.x)) {
        throw std::domain_error(
            "Hamiltonian state is outside the metric domain");
    }
    const Mat4 inverse = metric.contravariant(state.x);
    if (!all_finite(inverse)) {
        throw std::domain_error(
            "metric returned a non-finite inverse");
    }
    return inverse;
}

double contract_hamiltonian(
    const Mat4& inverse,
    const Covariant4& momentum) {
    double value = 0.0;
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            value += 0.5 * inverse[row][column] *
                     momentum.v[row] * momentum.v[column];
        }
    }
    if (!std::isfinite(value)) {
        throw std::domain_error(
            "Hamiltonian contraction is non-finite");
    }
    return value;
}

} // namespace

double hamiltonian(
    const Metric& metric,
    const PhaseSpaceState& state) {
    return contract_hamiltonian(
        validated_inverse(metric, state), state.p);
}

double hamiltonian_constraint_error(
    const Metric& metric,
    const PhaseSpaceState& state,
    GeodesicKind kind) {
    const Mat4 inverse = validated_inverse(metric, state);
    const double value = contract_hamiltonian(inverse, state.p);
    const double target = hamiltonian_target(kind);

    double absolute_scale_sum = 0.0;
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            absolute_scale_sum += std::fabs(
                inverse[row][column] *
                state.p.v[row] * state.p.v[column]);
        }
    }
    const double denominator =
        1.0 + std::fabs(target) + 0.5 * absolute_scale_sum;
    const double error = std::fabs(value - target) / denominator;
    if (!std::isfinite(error)) {
        throw std::domain_error(
            "Hamiltonian constraint error is non-finite");
    }
    return error;
}

HamiltonGeodesicRhs::HamiltonGeodesicRhs(
    const Metric& metric) noexcept
    : metric_(&metric) {}

PhaseSpaceDerivative HamiltonGeodesicRhs::operator()(
    const PhaseSpaceState& state) const {
    const Mat4 inverse = validated_inverse(*metric_, state);
    const auto derivatives =
        metric_->contravariant_derivatives(state.x);
    for (const Mat4& derivative : derivatives) {
        if (!all_finite(derivative)) {
            throw std::domain_error(
                "metric returned a non-finite inverse derivative");
        }
    }

    PhaseSpaceDerivative result{};
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            result.dx.v[row] +=
                inverse[row][column] * state.p.v[column];
        }
    }
    for (std::size_t coordinate = 0; coordinate < 4; ++coordinate) {
        for (std::size_t row = 0; row < 4; ++row) {
            for (std::size_t column = 0; column < 4; ++column) {
                result.dp.v[coordinate] -=
                    0.5 * derivatives[coordinate][row][column] *
                    state.p.v[row] * state.p.v[column];
            }
        }
    }
    if (!result.dx.v.all_finite() ||
        !result.dp.v.all_finite()) {
        throw std::domain_error(
            "Hamiltonian derivative is non-finite");
    }
    return result;
}

numerics::StateN<8> pack_phase_space(
    const PhaseSpaceState& state) {
    require_finite_state(state);
    numerics::StateN<8> packed{};
    for (std::size_t component = 0; component < 4; ++component) {
        packed[component] = state.x.v[component];
        packed[component + 4] = state.p.v[component];
    }
    return packed;
}

PhaseSpaceState unpack_phase_space(
    double affine,
    const numerics::StateN<8>& packed) {
    if (!std::isfinite(affine)) {
        throw std::domain_error(
            "affine parameter must be finite");
    }
    PhaseSpaceState state{};
    state.affine = affine;
    for (std::size_t component = 0; component < 4; ++component) {
        state.x.v[component] = packed[component];
        state.p.v[component] = packed[component + 4];
    }
    require_finite_state(state);
    return state;
}

} // namespace solar::relativity
