#include "kerr_separated_step.h"

#include "kerr_separated_potentials.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace solar::relativity::detail {
namespace {

double sign_of(SeparatedDirection direction) noexcept {
    return static_cast<double>(
        static_cast<int>(direction));
}

void fill_nan(numerics::StateN<5>& derivative) noexcept {
    derivative.fill(
        std::numeric_limits<double>::quiet_NaN());
}

bool finite_state(const numerics::StateN<5>& state) noexcept {
    return std::all_of(
        state.begin(),
        state.end(),
        [](double value) {
            return std::isfinite(value);
        });
}

} // namespace

KerrSeparatedStepAttempt attempt_kerr_separated_step(
    const KerrBoyerLindquistMetric& metric,
    const KerrConstants& constants,
    const KerrSeparatedState& current,
    double current_mino,
    double attempted_step,
    const KerrSeparatedConfig& config) {
    bool radial_forbidden = false;
    bool polar_forbidden = false;
    bool invalid_metric_point = false;
    bool non_finite_stage = false;
    double first_radial_forbidden_coordinate =
        std::numeric_limits<double>::quiet_NaN();
    double first_polar_forbidden_coordinate =
        std::numeric_limits<double>::quiet_NaN();
    const KerrSeparatedPotentials potentials(
        metric.mass(), metric.spin_length(), constants);

    const auto rhs =
        [&](double, const numerics::StateN<5>& state) {
            numerics::StateN<5> derivative{};
            if (!finite_state(state)) {
                non_finite_stage = true;
                fill_nan(derivative);
                return derivative;
            }

            const double radius = state[kRadius];
            const double mu = state[kMu];
            if (std::fabs(mu) >= 1.0) {
                invalid_metric_point = true;
                fill_nan(derivative);
                return derivative;
            }

            KerrSeparatedPotentialValues values;
            try {
                values = potentials.evaluate(radius, mu);
            } catch (const std::domain_error&) {
                non_finite_stage = true;
                fill_nan(derivative);
                return derivative;
            }

            if (current.radial_direction !=
                    SeparatedDirection::Locked &&
                values.radial < 0.0) {
                radial_forbidden = true;
                if (!std::isfinite(
                        first_radial_forbidden_coordinate)) {
                    first_radial_forbidden_coordinate = radius;
                }
            }
            if (current.polar_direction !=
                    SeparatedDirection::Locked &&
                values.polar < 0.0) {
                polar_forbidden = true;
                if (!std::isfinite(
                        first_polar_forbidden_coordinate)) {
                    first_polar_forbidden_coordinate = mu;
                }
            }
            if (radial_forbidden || polar_forbidden) {
                fill_nan(derivative);
                return derivative;
            }

            const double theta = std::acos(mu);
            const Contravariant4 x{
                Vec4{{
                    state[kTime],
                    radius,
                    theta,
                    state[kAzimuth],
                }}};
            if (!metric.valid_point(x)) {
                invalid_metric_point = true;
                fill_nan(derivative);
                return derivative;
            }

            const double one_minus_mu_sq =
                1.0 - mu * mu;
            const double radius_sq = radius * radius;
            const double spin = metric.spin_length();
            const double spin_sq = spin * spin;
            const double radial_momentum =
                constants.E * (radius_sq + spin_sq) -
                spin * constants.Lz;
            derivative[kTime] =
                (radius_sq + spin_sq) *
                    radial_momentum / values.delta +
                spin *
                    (constants.Lz -
                     spin * constants.E *
                         one_minus_mu_sq);
            derivative[kRadius] =
                sign_of(current.radial_direction) *
                std::sqrt(values.radial);
            derivative[kMu] =
                sign_of(current.polar_direction) *
                std::sqrt(values.polar);
            derivative[kAzimuth] =
                spin * radial_momentum / values.delta +
                constants.Lz / one_minus_mu_sq -
                spin * constants.E;
            derivative[kAffine] = values.sigma;

            if (!finite_state(derivative)) {
                non_finite_stage = true;
                fill_nan(derivative);
            }
            return derivative;
        };

    auto step = numerics::dopri5_step(
        current.values,
        current_mino,
        attempted_step,
        rhs,
        config.dopri5);
    return KerrSeparatedStepAttempt{
        std::move(step),
        radial_forbidden,
        polar_forbidden,
        invalid_metric_point,
        non_finite_stage,
        first_radial_forbidden_coordinate,
        first_polar_forbidden_coordinate,
    };
}

} // namespace solar::relativity::detail
