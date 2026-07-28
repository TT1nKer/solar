#include "solar/gravity.h"
#include "solar/constants.h"
#include <cmath>

namespace solar {

void NewtonianGravity::compute(
    const std::vector<Body>& bodies,
    double /*time*/,
    std::vector<Vec3>& acc) const
{
    size_t n = bodies.size();

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            Vec3 rij = bodies[j].state.pos - bodies[i].state.pos;
            double dist_sq = rij.norm_sq();
            double dist = std::sqrt(dist_sq);
            double dist_cubed = dist_sq * dist;

            Vec3 force_dir = rij / dist_cubed;

            acc[i] += force_dir * bodies[j].mu;
            acc[j] -= force_dir * bodies[i].mu;
        }
    }
}

double NewtonianGravity::potential_energy(
    const std::vector<Body>& bodies,
    double /*time*/) const
{
    double pe = 0.0;
    size_t n = bodies.size();

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            Vec3 rij = bodies[j].state.pos - bodies[i].state.pos;
            double dist = rij.norm();
            // Dynamics uses mu directly, so the matching inertial mass is mu/G.
            pe -= bodies[i].mu * bodies[j].mu / (constants::G * dist);
        }
    }

    return pe;
}

} // namespace solar
