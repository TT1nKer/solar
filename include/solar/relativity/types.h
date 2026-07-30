#pragma once

#include "solar/relativity/math.h"

#include <cstdint>
#include <limits>

namespace solar::relativity {

struct Contravariant4 {
    Vec4 v;
};

struct Covariant4 {
    Vec4 v;
};

enum class GeodesicKind {
    Null,
    TimelikeUnitMass,
};

constexpr bool is_valid_geodesic_kind(GeodesicKind kind) noexcept {
    return kind == GeodesicKind::Null ||
           kind == GeodesicKind::TimelikeUnitMass;
}

constexpr double hamiltonian_target(GeodesicKind kind) noexcept {
    return kind == GeodesicKind::Null
               ? 0.0
               : kind == GeodesicKind::TimelikeUnitMass
                     ? -0.5
                     : std::numeric_limits<double>::quiet_NaN();
}

constexpr bool has_proper_time(GeodesicKind kind) noexcept {
    return kind == GeodesicKind::TimelikeUnitMass;
}

struct PhaseSpaceState {
    double affine = 0.0;
    Contravariant4 x;
    Covariant4 p;
};

/**
 * Records a sampled Hamiltonian geodesic state and its invariants.
 *
 * Signature:
 *   (-,+,+,+)
 *
 * Units:
 *   Geometrized G=c=1.
 *
 * Model:
 *   Fixed background spacetime; test particle or photon. A null sample has no
 *   proper time, represented by quiet NaN rather than coordinate time.
 */
struct GeodesicSample {
    double affine = 0.0;
    double proper_time = std::numeric_limits<double>::quiet_NaN();
    Contravariant4 x;
    Covariant4 p;
    Contravariant4 tangent;
    double hamiltonian = 0.0;
    double energy = 0.0;
    double lz = 0.0;
    double carter_q = 0.0;
    std::uint32_t flags = 0;
};

} // namespace solar::relativity
