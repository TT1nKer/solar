#include "solar/relativity/geodesic_integrator.h"
#include "solar/relativity/kerr_constants.h"
#include "solar/relativity/kerr_shadow.h"
#include "solar/relativity/local_initialization.h"
#include "solar/relativity/observer.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

using namespace solar::relativity;

namespace {

struct ShadowRayBenchmark {
    bool all_rays_classified = true;
    bool all_rays_future_directed = true;
    double left_edge = std::numeric_limits<double>::quiet_NaN();
    double right_edge = std::numeric_limits<double>::quiet_NaN();
    double max_hamiltonian_error = 0.0;
    double max_carter_relative_error = 0.0;
    std::string failure;
};

std::optional<bool> classify_shadow_ray_termination(
    const GeodesicIntegrationResult& ray) {
    if (!ray.event) {
        return std::nullopt;
    }
    if (ray.event->event_index == 0 &&
        ray.diagnostics.reason ==
            TerminationReason::UserEvent) {
        return true;
    }
    if (ray.event->event_index == 1 &&
        ray.diagnostics.reason ==
            TerminationReason::Escaped) {
        return false;
    }
    return std::nullopt;
}

class ShadowRayClassifier {
public:
    ShadowRayClassifier(
        const KerrBoyerLindquistMetric& metric,
        const ObserverFrame& observer,
        double observer_radius)
        : metric_(metric),
          observer_(observer),
          observer_radius_(observer_radius),
          integrator_(metric) {
        const double inner_radius =
            metric_.outer_horizon_radius() + 1.0e-3;
        const double escape_radius = 1.1 * observer_radius_;
        events_ = {
            GeodesicEvent{
                "BL capture proxy",
                [inner_radius](const PhaseSpaceState& state) {
                    return state.x.v[1] - inner_radius;
                },
                EventDirection::Decreasing,
                TerminationReason::UserEvent,
                1.0e-10,
            },
            GeodesicEvent{
                "outer escape",
                [escape_radius](const PhaseSpaceState& state) {
                    return state.x.v[1] - escape_radius;
                },
                EventDirection::Increasing,
                TerminationReason::Escaped,
                1.0e-10,
            },
        };
    }

    std::optional<bool> is_captured(
        double alpha,
        ShadowRayBenchmark& benchmark) const {
        const double local_phi = -alpha / observer_radius_;
        const double radial_squared =
            1.0 - local_phi * local_phi;
        if (!std::isfinite(radial_squared) ||
            radial_squared < 0.0) {
            return fail(
                benchmark,
                "screen direction is outside the local sky");
        }

        const InitialStateResult initial =
            initialize_local_photon(
                metric_,
                observer_,
                Vec3{{
                    std::sqrt(radial_squared),
                    0.0,
                    local_phi,
                }});
        if (!initial) {
            return fail(
                benchmark,
                "local photon initialization failed: " +
                    initial.message);
        }
        benchmark.all_rays_future_directed &=
            std::isfinite(initial.measured_frequency) &&
            initial.measured_frequency > 0.0;

        GeodesicIntegrationConfig config =
            GeodesicIntegrationConfig::cpu_reference(
                GeodesicKind::Null,
                metric_.mass(),
                -0.5,
                2.0,
                4000.0);
        config.carter_evaluator =
            [this](const PhaseSpaceState& state) {
                return evaluate_kerr_constants(
                    metric_, state, GeodesicKind::Null).Q;
            };
        const GeodesicIntegrationResult ray =
            integrator_.integrate(
                *initial.state, config, events_);

        if (!std::isfinite(
                ray.diagnostics.max_constraint_error) ||
            !std::isfinite(
                ray.diagnostics.max_carter_rel_error)) {
            return fail(
                benchmark,
                "ray produced non-finite invariant diagnostics");
        }
        benchmark.max_hamiltonian_error = std::max(
            benchmark.max_hamiltonian_error,
            ray.diagnostics.max_constraint_error);
        benchmark.max_carter_relative_error = std::max(
            benchmark.max_carter_relative_error,
            ray.diagnostics.max_carter_rel_error);

        const std::optional<bool> classification =
            classify_shadow_ray_termination(ray);
        if (classification) {
            return classification;
        }
        return fail(
            benchmark,
            "ray terminated without an accepted classification: " +
                ray.diagnostics.message);
    }

private:
    static std::optional<bool> fail(
        ShadowRayBenchmark& benchmark,
        std::string message) {
        benchmark.all_rays_classified = false;
        if (benchmark.failure.empty()) {
            benchmark.failure = std::move(message);
        }
        return std::nullopt;
    }

    const KerrBoyerLindquistMetric& metric_;
    const ObserverFrame& observer_;
    double observer_radius_;
    GeodesicIntegrator integrator_;
    std::vector<GeodesicEvent> events_;
};

std::optional<double> locate_shadow_edge(
    const ShadowRayClassifier& classifier,
    double escaped_alpha,
    double captured_alpha,
    ShadowRayBenchmark& benchmark) {
    const std::optional<bool> escaped =
        classifier.is_captured(escaped_alpha, benchmark);
    const std::optional<bool> captured =
        classifier.is_captured(captured_alpha, benchmark);
    if (!escaped || !captured || *escaped || !*captured) {
        benchmark.all_rays_classified = false;
        if (benchmark.failure.empty()) {
            benchmark.failure =
                "shadow-edge bracket does not straddle the boundary";
        }
        return std::nullopt;
    }

    constexpr double screen_tolerance = 1.0e-3;
    while (std::fabs(escaped_alpha - captured_alpha) >
           screen_tolerance) {
        const double midpoint =
            0.5 * (escaped_alpha + captured_alpha);
        const std::optional<bool> midpoint_captured =
            classifier.is_captured(midpoint, benchmark);
        if (!midpoint_captured) {
            return std::nullopt;
        }
        if (*midpoint_captured) {
            captured_alpha = midpoint;
        } else {
            escaped_alpha = midpoint;
        }
    }
    return 0.5 * (escaped_alpha + captured_alpha);
}

ShadowRayBenchmark run_shadow_ray_benchmark() {
    ShadowRayBenchmark benchmark;
    constexpr double observer_radius = 1000.0;
    constexpr double half_pi = 1.5707963267948966;
    const KerrBoyerLindquistMetric metric(1.0, 0.5);
    const ObserverResult observer = make_zamo_observer(
        metric,
        Contravariant4{
            Vec4{{0.0, observer_radius, half_pi, 0.0}}});
    if (!observer) {
        benchmark.all_rays_classified = false;
        benchmark.failure =
            "distant ZAMO construction failed: " +
            observer.message;
        return benchmark;
    }

    const ShadowRayClassifier classifier(
        metric, *observer.frame, observer_radius);
    const std::optional<double> left = locate_shadow_edge(
        classifier, -8.0, 0.0, benchmark);
    const std::optional<double> right = locate_shadow_edge(
        classifier, 8.0, 0.0, benchmark);
    if (left) {
        benchmark.left_edge = *left;
    }
    if (right) {
        benchmark.right_edge = *right;
    }
    return benchmark;
}

void check(
    const std::string& name,
    bool condition,
    int& passed,
    int& failed) {
    std::cout << (condition ? "  PASS: " : "  FAIL: ")
              << name << "\n";
    condition ? ++passed : ++failed;
}

} // namespace

int main() {
    int passed = 0;
    int failed = 0;
    constexpr double half_pi = 1.5707963267948966;
    const KerrBoyerLindquistMetric metric(1.0, 0.5);
    const std::vector<ShadowCriticalPoint> analytic_curve =
        bardeen_shadow_curve(metric, half_pi, 65);
    const auto analytic_edges = std::minmax_element(
        analytic_curve.begin(),
        analytic_curve.end(),
        [](const ShadowCriticalPoint& left,
           const ShadowCriticalPoint& right) {
            return left.alpha < right.alpha;
        });
    const ShadowRayBenchmark benchmark =
        run_shadow_ray_benchmark();
    GeodesicIntegrationResult invalid_metric_ray{};
    invalid_metric_ray.diagnostics.reason =
        TerminationReason::InvalidMetricPoint;

    check(
        "CPU shadow rays all reach explicit classification events",
        benchmark.all_rays_classified,
        passed,
        failed);
    check(
        "CPU shadow rays remain future-directed",
        benchmark.all_rays_future_directed,
        passed,
        failed);
    check(
        "CPU left shadow edge converges to Bardeen edge",
        std::isfinite(benchmark.left_edge) &&
            std::fabs(
                benchmark.left_edge -
                analytic_edges.first->alpha) < 3.0e-2,
        passed,
        failed);
    check(
        "CPU right shadow edge converges to Bardeen edge",
        std::isfinite(benchmark.right_edge) &&
            std::fabs(
                benchmark.right_edge -
                analytic_edges.second->alpha) < 3.0e-2,
        passed,
        failed);
    check(
        "CPU shadow Hamiltonian gate",
        benchmark.max_hamiltonian_error < 1.0e-10,
        passed,
        failed);
    check(
        "CPU shadow Carter gate",
        benchmark.max_carter_relative_error < 1.0e-10,
        passed,
        failed);
    check(
        "invalid metric termination is never classified as capture",
        !classify_shadow_ray_termination(
             invalid_metric_ray).has_value(),
        passed,
        failed);

    if (!benchmark.failure.empty()) {
        std::cerr << "  shadow benchmark failure: "
                  << benchmark.failure << "\n";
    }
    std::cout << std::setprecision(17)
              << "  numerical_shadow_left="
              << benchmark.left_edge
              << " analytic_shadow_left="
              << analytic_edges.first->alpha
              << " numerical_shadow_right="
              << benchmark.right_edge
              << " analytic_shadow_right="
              << analytic_edges.second->alpha
              << " shadow_max_constraint="
              << benchmark.max_hamiltonian_error
              << " shadow_carter_rel="
              << benchmark.max_carter_relative_error
              << "\n\n=== Results: " << passed
              << " passed, " << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
