#pragma once

#include "solar/relativity/dual4.h"

#include <array>
#include <cmath>

namespace solar::relativity::detail {

template <typename Scalar>
Scalar kerr_ingoing_time_offset(
    const Scalar& radius,
    double mass,
    double outer_horizon,
    double inner_horizon) {
    using std::log;
    const double horizon_gap =
        outer_horizon - inner_horizon;
    return
        (2.0 * mass * outer_horizon / horizon_gap) *
            log(radius - outer_horizon) -
        (2.0 * mass * inner_horizon / horizon_gap) *
            log(radius - inner_horizon);
}

template <typename Scalar>
Scalar kerr_ingoing_azimuth_offset(
    const Scalar& radius,
    double spin_a,
    double outer_horizon,
    double inner_horizon) {
    using std::log;
    const double horizon_gap =
        outer_horizon - inner_horizon;
    return
        (spin_a / horizon_gap) *
        (log(radius - outer_horizon) -
         log(radius - inner_horizon));
}

template <typename Scalar>
std::array<Scalar, 4>
boyer_lindquist_to_kerr_schild_position(
    const std::array<Scalar, 4>& boyer_lindquist,
    double mass,
    double spin_a,
    double outer_horizon,
    double inner_horizon) {
    using std::cos;
    using std::sin;

    const Scalar& time = boyer_lindquist[0];
    const Scalar& radius = boyer_lindquist[1];
    const Scalar& theta = boyer_lindquist[2];
    const Scalar& phi = boyer_lindquist[3];
    const Scalar tilde_phi =
        phi + kerr_ingoing_azimuth_offset(
                  radius,
                  spin_a,
                  outer_horizon,
                  inner_horizon);
    const Scalar sine_theta = sin(theta);
    return {
        time + kerr_ingoing_time_offset(
                   radius,
                   mass,
                   outer_horizon,
                   inner_horizon),
        (radius * cos(tilde_phi) -
         spin_a * sin(tilde_phi)) * sine_theta,
        (radius * sin(tilde_phi) +
         spin_a * cos(tilde_phi)) * sine_theta,
        radius * cos(theta),
    };
}

} // namespace solar::relativity::detail
