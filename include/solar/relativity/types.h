#pragma once

namespace solar::relativity {

// Schwarzschild coordinates use geometric units: G = c = 1.
// Coordinates t and r are lengths in km; theta and phi are radians.
struct FourVector {
    double t = 0.0;
    double r = 0.0;
    double theta = 0.0;
    double phi = 0.0;
};

struct GeodesicState {
    FourVector position;
    FourVector tangent; // dx^mu / d(lambda), with lambda measured in km
};

enum class GeodesicKind {
    Timelike,
    Null,
};

struct ConservedQuantities {
    double energy = 0.0;
    double angular_momentum_z = 0.0;
};

} // namespace solar::relativity
