#include "solar/gr_correction.h"
#include "solar/constants.h"
#include <cmath>

namespace solar {

GRCorrection::GRCorrection(int central_body_index)
    : central_index_(central_body_index) {}

void GRCorrection::compute(
    const std::vector<Body>& bodies,
    double /*time*/,
    std::vector<Vec3>& acc) const
{
    size_t n = bodies.size();
    if (central_index_ < 0 || central_index_ >= static_cast<int>(n)) return;

    double mu = bodies[central_index_].mu;
    double c2 = constants::C_LIGHT * constants::C_LIGHT;

    // 1PN Schwarzschild correction for each body orbiting the central body
    // a_GR = (mu / (r^2 * c^2)) * [ (4*mu/r - v^2) * r_hat + 4*(r_hat . v) * v ]
    // Reference: Soffel et al. (2003), Einstein-Infeld-Hoffmann equations

    for (size_t i = 0; i < n; ++i) {
        if (static_cast<int>(i) == central_index_) continue;

        Vec3 r = bodies[i].state.pos - bodies[central_index_].state.pos;
        Vec3 v = bodies[i].state.vel - bodies[central_index_].state.vel;

        double r_mag = r.norm();
        double r2 = r_mag * r_mag;
        double v2 = v.norm_sq();
        double rdotv = r.dot(v);

        double coeff = mu / (r2 * c2);

        // Radial and velocity terms
        double radial_factor = (4.0 * mu / r_mag) - v2;
        double velocity_factor = 4.0 * rdotv / r_mag;

        Vec3 r_hat = r / r_mag;

        acc[i] += (r_hat * radial_factor + v * velocity_factor) * coeff;
    }
}

} // namespace solar
