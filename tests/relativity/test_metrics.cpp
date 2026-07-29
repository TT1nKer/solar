#include "solar/relativity/minkowski_metric.h"

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

void check_near(const char* name, double actual, double expected,
                double tolerance) {
    check(name, std::fabs(actual - expected) <= tolerance);
}

double derivative_max_abs(const std::array<Mat4, 4>& derivatives) {
    double maximum = 0.0;
    for (const auto& derivative : derivatives) {
        for (const auto& row : derivative) {
            for (const double component : row) {
                maximum = std::max(maximum, std::fabs(component));
            }
        }
    }
    return maximum;
}

} // namespace

int main() {
    const MinkowskiMetric metric;
    const Contravariant4 x{Vec4{{2.0, -3.0, 4.0, 5.0}}};

    check("Minkowski chart", metric.chart() == Chart::MinkowskiCartesian);
    check("Minkowski name", metric.name() == "minkowski");
    check("finite point accepted", metric.valid_point(x));

    const Mat4 covariant = metric.covariant(x);
    const Mat4 contravariant = metric.contravariant(x);
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            const double expected =
                row == column ? (row == 0 ? -1.0 : 1.0) : 0.0;
            check_near("Minkowski covariant component",
                       covariant[row][column], expected, 0.0);
            check_near("Minkowski inverse component",
                       contravariant[row][column], expected, 0.0);
        }
    }

    check_near("Minkowski derivatives are exactly zero",
               derivative_max_abs(metric.contravariant_derivatives(x)),
               0.0, 0.0);

    Contravariant4 invalid = x;
    invalid.v[2] = std::numeric_limits<double>::quiet_NaN();
    check("non-finite Minkowski point rejected", !metric.valid_point(invalid));
    try {
        (void)metric.covariant(invalid);
        check("invalid Minkowski evaluation throws", false);
    } catch (const std::domain_error&) {
        check("invalid Minkowski evaluation throws", true);
    }

    std::cout << "\n=== Results: " << passed << " passed, " << failed
              << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
