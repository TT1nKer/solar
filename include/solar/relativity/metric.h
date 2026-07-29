#pragma once

#include "solar/relativity/types.h"

#include <array>
#include <string>

namespace solar::relativity {

enum class Chart {
    MinkowskiCartesian,
    BoyerLindquist,
    KerrSchildCartesian,
};

const char* chart_name(Chart chart) noexcept;

/**
 * Defines a fixed-background spacetime metric.
 *
 * Signature:
 *   (-,+,+,+)
 *
 * Coordinates:
 *   Reported by chart(); component x^0 is always first.
 *
 * Units:
 *   Geometrized G=c=1.
 *
 * Model:
 *   Analytic fixed background; test particles/photons do not backreact.
 *
 * Failure:
 *   Matrix methods throw std::domain_error when x is outside valid_point().
 */
class Metric {
public:
    virtual ~Metric() = default;

    virtual Chart chart() const noexcept = 0;
    virtual std::string name() const = 0;

    virtual Mat4 covariant(const Contravariant4& x) const = 0;
    virtual Mat4 contravariant(const Contravariant4& x) const = 0;

    // result[coordinate_mu][alpha][beta]
    virtual std::array<Mat4, 4>
    contravariant_derivatives(const Contravariant4& x) const = 0;

    virtual bool valid_point(const Contravariant4& x) const noexcept = 0;
};

} // namespace solar::relativity
