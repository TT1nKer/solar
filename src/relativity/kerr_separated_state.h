#pragma once

#include "kerr_separated_potentials.h"

#include "solar/numerics/dopri5.h"
#include "solar/relativity/kerr_bl_metric.h"

#include <cstddef>
#include <stdexcept>

namespace solar::relativity::detail {

using KerrMinoState = numerics::StateN<5>;

inline constexpr std::size_t kTime = 0;
inline constexpr std::size_t kRadius = 1;
inline constexpr std::size_t kMu = 2;
inline constexpr std::size_t kAzimuth = 3;
inline constexpr std::size_t kAffine = 4;

enum class SeparatedDirection : int {
    Negative = -1,
    Locked = 0,
    Positive = 1,
};

struct KerrSeparatedState {
    KerrMinoState values;
    SeparatedDirection radial_direction;
    SeparatedDirection polar_direction;
};

struct KerrSeparatedInitialState {
    KerrSeparatedState state;
    KerrConstants constants;
};

class KerrSeparatedCriticalInitialState
    : public std::domain_error {
public:
    using std::domain_error::domain_error;
};

KerrSeparatedInitialState initialize_kerr_separated_state(
    const KerrBoyerLindquistMetric& metric,
    const PhaseSpaceState& initial,
    GeodesicKind kind,
    double normalized_potential_tolerance,
    double normalized_critical_tolerance,
    double integration_direction);

PhaseSpaceState reconstruct_kerr_phase_space(
    const KerrBoyerLindquistMetric& metric,
    const KerrConstants& constants,
    const KerrSeparatedState& state);

} // namespace solar::relativity::detail
