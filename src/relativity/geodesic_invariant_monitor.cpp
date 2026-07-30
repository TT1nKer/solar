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
      stationary_energy_evaluator_(
          config.stationary_energy_evaluator),
      axial_angular_momentum_evaluator_(
          config.axial_angular_momentum_evaluator),
      initial_state_(initial),
      initial_energy_(-initial.p.v[0]),
      initial_lz_(initial.p.v[3]) {}

bool GeodesicInvariantMonitor::evaluate(
    const char* label,
    const InvariantEvaluator& evaluator,
    const PhaseSpaceState& state,
    double default_value,
    double& value,
    std::string& failure_message) const {
    if (!evaluator) {
        value = default_value;
        return true;
    }
    try {
        value = evaluator(state);
    } catch (const std::exception& error) {
        failure_message =
            std::string(label) + " evaluator failed: " +
            error.what();
        return false;
    } catch (...) {
        failure_message =
            std::string(label) +
            " evaluator failed with a non-standard exception";
        return false;
    }
    if (!std::isfinite(value)) {
        failure_message =
            std::string(label) +
            " evaluator returned a non-finite value";
        return false;
    }
    return true;
}

bool GeodesicInvariantMonitor::initialize(
    IntegrationDiagnostics& diagnostics,
    std::string& failure_message) {
    if (monitor_energy_) {
        if (!evaluate(
                "stationary energy",
                stationary_energy_evaluator_,
                initial_state_,
                -initial_state_.p.v[0],
                initial_energy_,
                failure_message)) {
            return false;
        }
        diagnostics.max_energy_rel_error = 0.0;
    }
    if (monitor_lz_) {
        if (!evaluate(
                "axial angular momentum",
                axial_angular_momentum_evaluator_,
                initial_state_,
                initial_state_.p.v[3],
                initial_lz_,
                failure_message)) {
            return false;
        }
        diagnostics.max_lz_rel_error = 0.0;
    }
    if (!carter_evaluator_) {
        return true;
    }
    if (!evaluate(
            "Carter",
            carter_evaluator_,
            initial_state_,
            0.0,
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
    double energy = -state.p.v[0];
    if (monitor_energy_ &&
        !evaluate(
            "stationary energy",
            stationary_energy_evaluator_,
            state,
            energy,
            energy,
            failure_message)) {
        return false;
    }
    const double energy_error =
        monitor_energy_
            ? normalized_invariant_error(
                  energy, initial_energy_)
            : 0.0;

    double lz = state.p.v[3];
    if (monitor_lz_ &&
        !evaluate(
            "axial angular momentum",
            axial_angular_momentum_evaluator_,
            state,
            lz,
            lz,
            failure_message)) {
        return false;
    }
    const double lz_error =
        monitor_lz_
            ? normalized_invariant_error(
                  lz, initial_lz_)
            : 0.0;

    double carter = 0.0;
    if (carter_evaluator_ &&
        !evaluate(
            "Carter",
            carter_evaluator_,
            state,
            0.0,
            carter,
            failure_message)) {
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
