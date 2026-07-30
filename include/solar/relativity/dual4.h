#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace solar::relativity {

/**
 * A scalar value with forward derivatives in four coordinate directions.
 *
 * The derivative array order matches spacetime coordinates x^0 through x^3.
 * Arithmetic applies exact first-order product, quotient, and chain rules.
 * Domain errors are reported rather than clamped.
 *
 * Validation:
 *   tests/relativity/test_dual4.cpp
 */
struct Dual4 {
    double value = 0.0;
    std::array<double, 4> derivative{};

    static Dual4 variable(double value, std::size_t coordinate_index) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("Dual4 variable value must be finite");
        }
        if (coordinate_index >= 4) {
            throw std::out_of_range(
                "Dual4 coordinate index must be in [0,3]");
        }
        Dual4 result{value};
        result.derivative[coordinate_index] = 1.0;
        return result;
    }

    bool all_finite() const noexcept {
        if (!std::isfinite(value)) {
            return false;
        }
        for (const double partial : derivative) {
            if (!std::isfinite(partial)) {
                return false;
            }
        }
        return true;
    }
};

inline Dual4 operator+(const Dual4& left, const Dual4& right) {
    Dual4 result{left.value + right.value};
    for (std::size_t index = 0; index < 4; ++index) {
        result.derivative[index] =
            left.derivative[index] + right.derivative[index];
    }
    return result;
}

inline Dual4 operator-(const Dual4& left, const Dual4& right) {
    Dual4 result{left.value - right.value};
    for (std::size_t index = 0; index < 4; ++index) {
        result.derivative[index] =
            left.derivative[index] - right.derivative[index];
    }
    return result;
}

inline Dual4 operator-(const Dual4& input) {
    Dual4 result{-input.value};
    for (std::size_t index = 0; index < 4; ++index) {
        result.derivative[index] = -input.derivative[index];
    }
    return result;
}

inline Dual4 operator*(const Dual4& left, const Dual4& right) {
    Dual4 result{left.value * right.value};
    for (std::size_t index = 0; index < 4; ++index) {
        result.derivative[index] =
            left.derivative[index] * right.value +
            left.value * right.derivative[index];
    }
    return result;
}

inline Dual4 operator/(const Dual4& numerator, const Dual4& denominator) {
    if (!std::isfinite(denominator.value) || denominator.value == 0.0) {
        throw std::domain_error(
            "Dual4 denominator value must be finite and nonzero");
    }
    Dual4 result{numerator.value / denominator.value};
    for (std::size_t index = 0; index < 4; ++index) {
        result.derivative[index] =
            (numerator.derivative[index] -
             result.value * denominator.derivative[index]) /
            denominator.value;
    }
    if (!result.all_finite()) {
        throw std::domain_error(
            "Dual4 quotient produced a non-finite value or derivative");
    }
    return result;
}

inline Dual4 operator+(const Dual4& left, double right) {
    return left + Dual4{right};
}

inline Dual4 operator+(double left, const Dual4& right) {
    return Dual4{left} + right;
}

inline Dual4 operator-(const Dual4& left, double right) {
    return left - Dual4{right};
}

inline Dual4 operator-(double left, const Dual4& right) {
    return Dual4{left} - right;
}

inline Dual4 operator*(const Dual4& left, double right) {
    return left * Dual4{right};
}

inline Dual4 operator*(double left, const Dual4& right) {
    return Dual4{left} * right;
}

inline Dual4 operator/(const Dual4& left, double right) {
    return left / Dual4{right};
}

inline Dual4 operator/(double left, const Dual4& right) {
    return Dual4{left} / right;
}

inline Dual4 sqrt(const Dual4& input) {
    if (input.value < 0.0 || !std::isfinite(input.value)) {
        throw std::domain_error(
            "Dual4 square root requires a finite non-negative value");
    }
    if (input.value == 0.0) {
        for (const double partial : input.derivative) {
            if (partial != 0.0) {
                throw std::domain_error(
                    "Dual4 square-root derivative is singular at zero");
            }
        }
        return Dual4{0.0};
    }
    const double root = std::sqrt(input.value);
    Dual4 result{root};
    for (std::size_t index = 0; index < 4; ++index) {
        result.derivative[index] =
            input.derivative[index] / (2.0 * root);
    }
    return result;
}

inline Dual4 sin(const Dual4& input) {
    Dual4 result{std::sin(input.value)};
    const double slope = std::cos(input.value);
    for (std::size_t index = 0; index < 4; ++index) {
        result.derivative[index] = slope * input.derivative[index];
    }
    return result;
}

inline Dual4 cos(const Dual4& input) {
    Dual4 result{std::cos(input.value)};
    const double slope = -std::sin(input.value);
    for (std::size_t index = 0; index < 4; ++index) {
        result.derivative[index] = slope * input.derivative[index];
    }
    return result;
}

inline Dual4 log(const Dual4& input) {
    if (input.value <= 0.0 || !std::isfinite(input.value)) {
        throw std::domain_error(
            "Dual4 logarithm requires a finite positive value");
    }
    Dual4 result{std::log(input.value)};
    for (std::size_t index = 0; index < 4; ++index) {
        result.derivative[index] =
            input.derivative[index] / input.value;
    }
    return result;
}

inline Dual4 atan2(const Dual4& y, const Dual4& x) {
    const double radius_squared = x.value * x.value + y.value * y.value;
    if (!std::isfinite(radius_squared) || radius_squared == 0.0) {
        throw std::domain_error(
            "Dual4 atan2 requires a finite point away from the origin");
    }
    Dual4 result{std::atan2(y.value, x.value)};
    for (std::size_t index = 0; index < 4; ++index) {
        result.derivative[index] =
            (x.value * y.derivative[index] -
             y.value * x.derivative[index]) /
            radius_squared;
    }
    return result;
}

} // namespace solar::relativity
