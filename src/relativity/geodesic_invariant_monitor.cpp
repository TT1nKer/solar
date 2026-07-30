#include "geodesic_invariant_monitor.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <utility>

namespace solar::relativity::detail {
namespace {

double normalized_invariant_error(
    double current,
    double initial) {
    return std::fabs(current - initial) /
           std::max(1.0, std::fabs(initial));
}

} // namespace

GeodesicInvariantMonitor::GeodesicInvariantMonitor(
    const PhaseSpaceState& initial,
    const GeodesicIntegrationConfig& config)
    : monitor_energy_(config.monitor_energy),
      monitor_lz_(config.monitor_lz),
      carter_evaluator_(config.carter_evaluator),
      initial_state_(initial),
      initial_energy_(-initial.p.v[0]),
      initial_lz_(initial.p.v[3]) {}

bool GeodesicInvariantMonitor::evaluate_carter(
    const PhaseSpaceState& state,
    double& value,
    std::string& failure_message) const {
    try {
        value = carter_evaluator_(state);
    } catch (const std::exception& error) {
        failure_message =
            std::string("Carter evaluator failed: ") +
            error.what();
        return false;
    } catch (...) {
        failure_message =
            "Carter evaluator failed with a non-standard exception";
        return false;
    }
    if (!std::isfinite(value)) {
        failure_message =
            "Carter evaluator returned a non-finite value";
        return false;
    }
    return true;
}

bool GeodesicInvariantMonitor::initialize(
    IntegrationDiagnostics& diagnostics,
    std::string& failure_message) {
    if (monitor_energy_) {
        diagnostics.max_energy_rel_error = 0.0;
    }
    if (monitor_lz_) {
        diagnostics.max_lz_rel_error = 0.0;
    }
    if (!carter_evaluator_) {
        return true;
    }
    if (!evaluate_carter(
            initial_state_,
            initial_carter_,
            failure_message)) {
        return false;
    }
    diagnostics.max_carter_rel_error = 0.0;
    diagnostics.max_carter_abs_error = 0.0;
    return true;
}

bool GeodesicInvariantMonitor::update(
    const PhaseSpaceState& state,
    IntegrationDiagnostics& diagnostics,
    std::string& failure_message) const {
    const double energy_error = monitor_energy_
        ? normalized_invariant_error(
              -state.p.v[0], initial_energy_)
        : 0.0;
    const double lz_error = monitor_lz_
        ? normalized_invariant_error(
              state.p.v[3], initial_lz_)
        : 0.0;

    double carter = 0.0;
    if (carter_evaluator_ &&
        !evaluate_carter(state, carter, failure_message)) {
        return false;
    }

    if (monitor_energy_) {
        diagnostics.max_energy_rel_error = std::max(
            diagnostics.max_energy_rel_error, energy_error);
    }
    if (monitor_lz_) {
        diagnostics.max_lz_rel_error = std::max(
            diagnostics.max_lz_rel_error, lz_error);
    }
    if (carter_evaluator_) {
        const double absolute =
            std::fabs(carter - initial_carter_);
        diagnostics.max_carter_rel_error = std::max(
            diagnostics.max_carter_rel_error,
            absolute /
                std::max(1.0, std::fabs(initial_carter_)));
        diagnostics.max_carter_abs_error = std::max(
            diagnostics.max_carter_abs_error, absolute);
    }
    return true;
}

} // namespace solar::relativity::detail
