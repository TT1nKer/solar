#include "solar/relativity/kerr_orbits.h"
#include "solar/relativity/kerr_shadow.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace solar::relativity;

namespace {

int passed = 0;
int failed = 0;

void check(const std::string& name, bool condition) {
    if (condition) {
        std::cout << "  PASS: " << name << "\n";
        ++passed;
    } else {
        std::cerr << "  FAIL: " << name << "\n";
        ++failed;
    }
}

void check_near(
    const std::string& name,
    double actual,
    double expected,
    double tolerance) {
    const bool within_tolerance =
        std::isfinite(actual) &&
        std::fabs(actual - expected) <= tolerance;
    if (!within_tolerance) {
        std::cerr << std::setprecision(17)
                  << "    actual=" << actual
                  << " expected=" << expected
                  << " error=" << std::fabs(actual - expected)
                  << " tolerance=" << tolerance << "\n";
    }
    check(name, within_tolerance);
}

} // namespace

int main() {
    constexpr double half_pi = 1.5707963267948966;
    constexpr std::size_t samples = 65;

    const KerrBoyerLindquistMetric schwarzschild(2.0, 0.0);
    const std::vector<ShadowCriticalPoint> schwarzschild_curve =
        bardeen_shadow_curve(
            schwarzschild, half_pi, samples);
    check(
        "Schwarzschild curve has closed-curve sample count",
        schwarzschild_curve.size() == 2 * samples - 2);
    for (const ShadowCriticalPoint& point :
         schwarzschild_curve) {
        check_near(
            "Schwarzschild shadow circle",
            point.alpha * point.alpha +
                point.beta * point.beta,
            108.0,
            5.0e-13);
        check_near(
            "Schwarzschild photon radius",
            point.photon_radius,
            6.0,
            0.0);
    }

    const KerrBoyerLindquistMetric tiny_spin(1.0, 1.0e-8);
    const std::vector<ShadowCriticalPoint> tiny_curve =
        bardeen_shadow_curve(
            tiny_spin, half_pi, 9);
    for (const ShadowCriticalPoint& point : tiny_curve) {
        check_near(
            "tiny-spin branch uses stable Schwarzschild limit",
            point.alpha * point.alpha +
                point.beta * point.beta,
            27.0,
            2.0e-13);
        check_near(
            "tiny-spin photon radius uses Schwarzschild limit",
            point.photon_radius,
            3.0,
            0.0);
    }

    const KerrBoyerLindquistMetric positive_spin(1.0, 0.5);
    const std::vector<ShadowCriticalPoint> positive_curve =
        bardeen_shadow_curve(
            positive_spin, half_pi, samples);
    check(
        "equatorial Kerr curve keeps every physical sample",
        positive_curve.size() == 2 * samples - 2);
    const auto alpha_extrema = std::minmax_element(
        positive_curve.begin(),
        positive_curve.end(),
        [](const ShadowCriticalPoint& left,
           const ShadowCriticalPoint& right) {
            return left.alpha < right.alpha;
        });
    check_near(
        "positive-spin left shadow edge",
        alpha_extrema.first->alpha,
        -4.096266658713869,
        4.0e-15);
    check_near(
        "positive-spin right shadow edge",
        alpha_extrema.second->alpha,
        6.138155724715452,
        6.0e-15);

    const double prograde_photon =
        kerr_equatorial_photon_radius(
            positive_spin, OrbitSense::Prograde);
    const double retrograde_photon =
        kerr_equatorial_photon_radius(
            positive_spin, OrbitSense::Retrograde);
    for (std::size_t index = 0;
         index < positive_curve.size();
         ++index) {
        const ShadowCriticalPoint& point =
            positive_curve[index];
        check(
            "Kerr shadow point is finite",
            std::isfinite(point.alpha) &&
                std::isfinite(point.beta) &&
                std::isfinite(point.photon_radius));
        check(
            "Kerr shadow photon radius is physical",
            point.photon_radius >= prograde_photon &&
                point.photon_radius <= retrograde_photon);
    }
    for (std::size_t upper = 1;
         upper + 1 < samples;
         ++upper) {
        const std::size_t reflected =
            2 * samples - 2 - upper;
        check_near(
            "Kerr shadow reflected alpha",
            positive_curve[reflected].alpha,
            positive_curve[upper].alpha,
            0.0);
        check_near(
            "Kerr shadow reflected beta",
            positive_curve[reflected].beta,
            -positive_curve[upper].beta,
            0.0);
    }

    constexpr double inclined_view = 1.0471975511965976;
    constexpr std::size_t inclined_samples = 33;
    const std::vector<ShadowCriticalPoint> inclined_curve =
        bardeen_shadow_curve(
            positive_spin,
            inclined_view,
            inclined_samples);
    check(
        "inclined Kerr curve includes both visible tips",
        inclined_curve.size() == 2 * inclined_samples - 2);
    if (inclined_curve.size() ==
        2 * inclined_samples - 2) {
        check_near(
            "inclined left tip photon radius",
            inclined_curve.front().photon_radius,
            2.408224077239731,
            2.0e-15);
        check_near(
            "inclined left tip alpha",
            inclined_curve.front().alpha,
            -4.230999071965547,
            4.0e-15);
        check_near(
            "inclined left tip beta",
            inclined_curve.front().beta,
            0.0,
            0.0);
        check_near(
            "inclined right tip photon radius",
            inclined_curve[inclined_samples - 1]
                .photon_radius,
            3.4421204900209483,
            3.0e-15);
        check_near(
            "inclined right tip alpha",
            inclined_curve[inclined_samples - 1].alpha,
            6.003824294846052,
            6.0e-15);
        check_near(
            "inclined right tip beta",
            inclined_curve[inclined_samples - 1].beta,
            0.0,
            0.0);
    }

    const std::vector<ShadowCriticalPoint>
        minimum_inclined_curve =
            bardeen_shadow_curve(
                positive_spin, inclined_view, 2);
    check(
        "two inclined samples still return both physical tips",
        minimum_inclined_curve.size() == 2);

    bool near_equatorial_resolved = false;
    try {
        const std::vector<ShadowCriticalPoint>
            near_equatorial_curve =
                bardeen_shadow_curve(
                    positive_spin,
                    half_pi + 1.0e-10,
                    2);
        near_equatorial_resolved =
            near_equatorial_curve.size() == 2;
    } catch (const std::domain_error&) {
    }
    check(
        "near-equatorial visible tips remain resolvable",
        near_equatorial_resolved);

    constexpr double near_axis_view = 1.0e-3;
    const std::vector<ShadowCriticalPoint> near_axis_curve =
        bardeen_shadow_curve(
            positive_spin, near_axis_view, 2);
    check(
        "near-axis Kerr curve resolves its narrow visible interval",
        near_axis_curve.size() == 2);
    if (near_axis_curve.size() == 2) {
        check_near(
            "near-axis left tip photon radius",
            near_axis_curve[0].photon_radius,
            2.882606664595796,
            3.0e-15);
        check_near(
            "near-axis right tip photon radius",
            near_axis_curve[1].photon_radius,
            2.883828927101564,
            3.0e-15);
        check_near(
            "near-axis left tip alpha",
            near_axis_curve[0].alpha,
            -5.119500178989135,
            7.0e-13);
        check_near(
            "near-axis right tip alpha",
            near_axis_curve[1].alpha,
            5.121562190711024,
            7.0e-13);
    }

    const KerrBoyerLindquistMetric negative_spin(1.0, -0.5);
    const std::vector<ShadowCriticalPoint> negative_curve =
        bardeen_shadow_curve(
            negative_spin, half_pi, samples);
    check(
        "signed-spin curves have matching sample count",
        negative_curve.size() == positive_curve.size());
    for (std::size_t index = 0;
         index < positive_curve.size();
         ++index) {
        check_near(
            "negative spin mirrors alpha",
            negative_curve[index].alpha,
            -positive_curve[index].alpha,
            5.0e-15);
        check_near(
            "negative spin preserves beta",
            negative_curve[index].beta,
            positive_curve[index].beta,
            5.0e-15);
    }

    const KerrBoyerLindquistMetric scaled_spin(3.0, 0.5);
    const std::vector<ShadowCriticalPoint> scaled_curve =
        bardeen_shadow_curve(
            scaled_spin, half_pi, samples);
    check(
        "mass-scaled curve has matching sample count",
        scaled_curve.size() == positive_curve.size());
    for (std::size_t index = 0;
         index < positive_curve.size();
         ++index) {
        check_near(
            "shadow alpha scales with mass",
            scaled_curve[index].alpha,
            3.0 * positive_curve[index].alpha,
            2.0e-14);
        check_near(
            "shadow beta scales with mass",
            scaled_curve[index].beta,
            3.0 * positive_curve[index].beta,
            2.0e-14);
        check_near(
            "spherical photon radius scales with mass",
            scaled_curve[index].photon_radius,
            3.0 * positive_curve[index].photon_radius,
            3.0e-15);
    }

    for (const double invalid_inclination : {
             0.0,
             3.1415926535897932,
             std::numeric_limits<double>::quiet_NaN(),
         }) {
        try {
            static_cast<void>(bardeen_shadow_curve(
                positive_spin,
                invalid_inclination,
                samples));
            check("invalid shadow inclination rejected", false);
        } catch (const std::invalid_argument&) {
            check("invalid shadow inclination rejected", true);
        }
    }
    try {
        static_cast<void>(bardeen_shadow_curve(
            positive_spin, half_pi, 1));
        check("insufficient shadow samples rejected", false);
    } catch (const std::invalid_argument&) {
        check("insufficient shadow samples rejected", true);
    }
    try {
        const KerrBoyerLindquistMetric huge_metric(
            std::numeric_limits<double>::max() / 2.0,
            0.5);
        static_cast<void>(bardeen_shadow_curve(
            huge_metric, half_pi, samples));
        check("overflowing shadow output rejected", false);
    } catch (const std::overflow_error&) {
        check("overflowing shadow output rejected", true);
    }

    std::cout << "\n=== Results: " << passed
              << " passed, " << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
