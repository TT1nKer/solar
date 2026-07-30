#include "solar/relativity/geodesic_integrator.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace solar::relativity;

namespace {

int passed = 0;
int failed = 0;

void check(const char* name, bool condition) {
    std::cout << (condition ? "  PASS: " : "  FAIL: ") << name << '\n';
    condition ? ++passed : ++failed;
}

class DomainLimitedMinkowskiMetric final : public Metric {
public:
    Chart chart() const noexcept override {
        return Chart::MinkowskiCartesian;
    }

    std::string name() const override {
        return "domain-limited-minkowski-test";
    }

    Mat4 covariant(const Contravariant4& x) const override {
        require_valid(x);
        return diagonal();
    }

    Mat4 contravariant(const Contravariant4& x) const override {
        require_valid(x);
        return diagonal();
    }

    std::array<Mat4, 4> contravariant_derivatives(
        const Contravariant4& x) const override {
        require_valid(x);
        return {};
    }

    bool valid_point(const Contravariant4& x) const noexcept override {
        return x.v.all_finite() && x.v[1] < 1.0;
    }

private:
    static Mat4 diagonal() {
        Mat4 metric{};
        metric[0][0] = -1.0;
        metric[1][1] = 1.0;
        metric[2][2] = 1.0;
        metric[3][3] = 1.0;
        return metric;
    }

    void require_valid(const Contravariant4& x) const {
        if (!valid_point(x)) {
            throw std::domain_error(
                "point is outside the test metric domain");
        }
    }
};

PhaseSpaceState photon() {
    return PhaseSpaceState{
        0.0,
        Contravariant4{Vec4{{0.0, 0.0, 0.0, 0.0}}},
        Covariant4{Vec4{{-1.0, 1.0, 0.0, 0.0}}},
    };
}

GeodesicIntegrationConfig config() {
    auto result = GeodesicIntegrationConfig::cpu_reference(
        GeodesicKind::Null, 1.0, 2.0, 2.0, 2.0);
    result.min_step = 1.0e-2;
    result.max_rejections_per_step = 8;
    result.max_total_steps = 1000;
    return result;
}

} // namespace

int main() {
    const DomainLimitedMinkowskiMetric metric;
    const GeodesicIntegrator integrator(metric);

    const auto domain_result =
        integrator.integrate(photon(), config());
    check("invalid trial domain is reported explicitly",
          domain_result.diagnostics.reason ==
              TerminationReason::InvalidMetricPoint);
    check("invalid trial domain causes rejected steps",
          domain_result.diagnostics.rejected_steps > 0);
    check("smaller valid steps are accepted before domain stop",
          domain_result.diagnostics.accepted_steps > 0);
    check("domain stop remains inside the valid chart",
          metric.valid_point(domain_result.final_state.x));
    check("BL-style invalid domain is not called a horizon",
          domain_result.diagnostics.reason !=
              TerminationReason::HorizonCrossing);

    const GeodesicEvent initial_event{
        "initial boundary",
        [](const PhaseSpaceState& state) {
            return state.x.v[1];
        },
        EventDirection::Any,
        TerminationReason::UserEvent,
        1.0e-12,
    };
    const auto initial_event_result = integrator.integrate(
        photon(), config(), {initial_event});
    check("initial any-direction event terminates immediately",
          initial_event_result.diagnostics.reason ==
              TerminationReason::UserEvent);
    check("initial event performs no rejected trial",
          initial_event_result.diagnostics.rejected_steps == 0);
    check("initial event performs no accepted trial",
          initial_event_result.diagnostics.accepted_steps == 0);

    const GeodesicEvent invalid_directed_event{
        "missing directed event function",
        {},
        EventDirection::Increasing,
        TerminationReason::UserEvent,
        1.0e-12,
    };
    const auto invalid_event_contract = integrator.integrate(
        photon(), config(), {invalid_directed_event});
    check("invalid directed event contract fails before integration",
          invalid_event_contract.diagnostics.reason ==
              TerminationReason::EventRootFailure);
    check("invalid directed event contract performs no rejected trial",
          invalid_event_contract.diagnostics.rejected_steps == 0);
    check("invalid directed event contract performs no accepted trial",
          invalid_event_contract.diagnostics.accepted_steps == 0);

    int premature_invariant_calls = 0;
    auto invalid_event_monitor_config = config();
    invalid_event_monitor_config.carter_evaluator =
        [&premature_invariant_calls](
            const PhaseSpaceState&) {
            ++premature_invariant_calls;
            return 0.0;
        };
    const auto invalid_event_with_monitor = integrator.integrate(
        photon(),
        invalid_event_monitor_config,
        {invalid_directed_event});
    check("invalid event still reports contract failure with monitor",
          invalid_event_with_monitor.diagnostics.reason ==
              TerminationReason::EventRootFailure);
    check("invalid event is validated before invariant callback",
          premature_invariant_calls == 0);

    GeodesicEvent unknown_direction_event = initial_event;
    unknown_direction_event.direction =
        static_cast<EventDirection>(99);
    const auto unknown_direction_contract = integrator.integrate(
        photon(), config(), {unknown_direction_event});
    check("unknown event direction fails before integration",
          unknown_direction_contract.diagnostics.reason ==
              TerminationReason::EventRootFailure);
    check("unknown event direction performs no rejected trial",
          unknown_direction_contract.diagnostics.rejected_steps == 0);
    check("unknown event direction performs no accepted trial",
          unknown_direction_contract.diagnostics.accepted_steps == 0);

    GeodesicEvent non_finite_event{
        "non-finite internal event",
        [](const PhaseSpaceState& state) {
            if (state.x.v[1] > 0.1 &&
                state.x.v[1] < 0.9) {
                return std::numeric_limits<double>::quiet_NaN();
            }
            return state.x.v[1] - 0.3;
        },
        EventDirection::Any,
        TerminationReason::UserEvent,
        1.0e-12,
    };
    auto event_config = config();
    event_config.initial_step = 0.9;
    event_config.max_step = 0.9;
    event_config.max_affine = 0.9;
    event_config.min_step = 1.0e-12;
    const auto event_failure = integrator.integrate(
        photon(), event_config, {non_finite_event});
    check("event root failure is explicit",
          event_failure.diagnostics.reason ==
              TerminationReason::EventRootFailure);
    check("failed event step is not accepted",
          event_failure.diagnostics.accepted_steps == 0);

    std::cout << "\n=== Results: " << passed << " passed, "
              << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
