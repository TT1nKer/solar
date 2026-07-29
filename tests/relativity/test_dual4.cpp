#include "solar/relativity/dual4.h"

#include <cmath>
#include <iostream>
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
                double tolerance = 1.0e-14) {
    check(name, std::fabs(actual - expected) <= tolerance);
}

template <typename Operation>
void check_domain_error(const char* name, Operation operation) {
    try {
        operation();
        check(name, false);
    } catch (const std::domain_error&) {
        check(name, true);
    } catch (...) {
        check(name, false);
    }
}

} // namespace

int main() {
    const Dual4 x = Dual4::variable(2.0, 0);
    const Dual4 y = Dual4::variable(3.0, 1);
    const Dual4 f = x * x * y + sin(y);

    check_near("multivariable value", f.value, 12.0 + std::sin(3.0));
    check_near("multivariable df/dx", f.derivative[0], 12.0);
    check_near("multivariable df/dy", f.derivative[1],
               4.0 + std::cos(3.0));
    check_near("unused derivative remains zero", f.derivative[2], 0.0);

    const Dual4 quotient = x / y;
    check_near("quotient value", quotient.value, 2.0 / 3.0);
    check_near("quotient dx", quotient.derivative[0], 1.0 / 3.0);
    check_near("quotient dy", quotient.derivative[1], -2.0 / 9.0);

    const Dual4 positive = Dual4::variable(4.0, 2);
    const Dual4 chain = log(sqrt(positive));
    check_near("sqrt-log value", chain.value, std::log(2.0));
    check_near("sqrt-log derivative", chain.derivative[2], 0.125);

    const Dual4 zero_root = sqrt(Dual4{0.0});
    check_near("constant zero square root value", zero_root.value, 0.0);
    check_near("constant zero square root derivative",
               zero_root.derivative[0], 0.0);

    const Dual4 angle = atan2(y, x);
    check_near("atan2 value", angle.value, std::atan2(3.0, 2.0));
    check_near("atan2 dx", angle.derivative[0], -3.0 / 13.0);
    check_near("atan2 dy", angle.derivative[1], 2.0 / 13.0);

    const Dual4 trig = sin(x) * sin(x) + cos(x) * cos(x);
    check_near("trigonometric identity value", trig.value, 1.0);
    check_near("trigonometric identity derivative", trig.derivative[0], 0.0);
    check("finite Dual4 detected", trig.all_finite());

    try {
        (void)Dual4::variable(1.0, 4);
        check("invalid variable index rejected", false);
    } catch (const std::out_of_range&) {
        check("invalid variable index rejected", true);
    }

    check_domain_error("Dual4 zero division rejected", [x] {
        (void)(x / 0.0);
    });
    check_domain_error("negative square root rejected", [] {
        (void)sqrt(Dual4{-1.0});
    });
    check_domain_error("singular square-root derivative rejected", [] {
        (void)sqrt(Dual4::variable(0.0, 0));
    });
    check_domain_error("non-positive logarithm rejected", [] {
        (void)log(Dual4{0.0});
    });
    check_domain_error("atan2 origin rejected", [] {
        (void)atan2(Dual4{0.0}, Dual4{0.0});
    });

    std::cout << "\n=== Results: " << passed << " passed, " << failed
              << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
