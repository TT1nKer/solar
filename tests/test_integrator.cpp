#include "solar/integrator.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

int passed = 0;
int failed = 0;

void check(const char* name, bool condition) {
    std::cout << (condition ? "  PASS: " : "  FAIL: ") << name << '\n';
    condition ? ++passed : ++failed;
}

void check_near(const char* name, double actual, double expected,
                double tolerance) {
    check(name, std::isfinite(actual) &&
                    std::fabs(actual - expected) <= tolerance);
}

template <typename Exception, typename Action>
void check_throws(const char* name, Action action) {
    try {
        action();
        check(name, false);
    } catch (const Exception&) {
        check(name, true);
    }
}

} // namespace

int main() {
    const std::vector<double> initial{1.0, -2.0};
    const auto rhs =
        [](double time, const std::vector<double>& state) {
            return std::vector<double>{
                state[0] + time,
                -0.5 * state[1] + 2.0 * time,
            };
        };

    const auto result = solar::dopri5_generic_step(
        initial, 0.3, 0.1, rhs, 1.0e-9, 1.0e-9);

    check("generic result preserves dimension", result.y.size() == 2);
    check_near("generic first component",
               result.y[0], 1.1418931121666667, 3.0e-15);
    check_near("generic second component",
               result.y[1], -1.834098762375, 3.0e-15);
    check_near("generic maximum normalized error",
               result.error, 8.3354999831097523, 2.0e-13);
    check_near("generic next step",
               result.dt_next, 0.058891983053898568, 2.0e-15);
    check_near("generic used step", result.dt_used, 0.1, 0.0);
    check("generic strict step rejected", !result.accepted);

    check_throws<std::invalid_argument>(
        "generic empty state rejected", [&] {
            (void)solar::dopri5_generic_step(
                {}, 0.0, 0.1, rhs, 1.0e-9, 1.0e-9);
        });
    check_throws<std::invalid_argument>(
        "generic non-positive tolerance rejected", [&] {
            (void)solar::dopri5_generic_step(
                initial, 0.0, 0.1, rhs, 0.0, 1.0e-9);
        });
    check_throws<std::invalid_argument>(
        "generic step must advance time", [&] {
            (void)solar::dopri5_generic_step(
                initial, 1.0e20, 1.0, rhs, 1.0e-9, 1.0e-9);
        });
    check_throws<std::domain_error>(
        "generic non-finite state rejected", [&] {
            const std::vector<double> invalid{
                std::numeric_limits<double>::quiet_NaN(),
                0.0,
            };
            (void)solar::dopri5_generic_step(
                invalid, 0.0, 0.1, rhs, 1.0e-9, 1.0e-9);
        });

    std::cout << "\n=== Results: " << passed << " passed, "
              << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
