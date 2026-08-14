#include "solar/dynamics/pn_gravity.h"

#include "solar/constants.h"

#include <cmath>

namespace solar {
namespace dynamics {

PostNewtonianGravity::PostNewtonianGravity(double speed_of_light_km_s)
    : speed_of_light_km_s_(speed_of_light_km_s) {}

void PostNewtonianGravity::compute(
    const std::vector<Body>& bodies,
    double /*time*/,
    std::vector<Vec3>& acc) const {
    const double c = speed_of_light_km_s_ > 0.0
        ? speed_of_light_km_s_
        : constants::C_LIGHT;  // km/s
    const double c_sq = c * c;

    for (std::size_t i = 0; i < bodies.size(); ++i) {
        const Vec3& position = bodies[i].state.pos;
        const Vec3& velocity = bodies[i].state.vel;
        const double v_sq = velocity.norm_sq();

        for (std::size_t j = 0; j < bodies.size(); ++j) {
            if (i == j) continue;
            Vec3 direction = bodies[j].state.pos - position;
            const double r = direction.norm();
            if (r == 0.0) continue;
            direction = direction * (1.0 / r);

            const double mu_j = bodies[j].mu;
            const double newtonian = mu_j / (r * r);
            const double compactness = mu_j / (r * c_sq);
            const double v_dot_n = velocity.dot(direction);

            // n points from the source toward the accelerated body in
            // the standard formula; here direction points the other way,
            // so both terms flip sign.
            acc[i] += direction *
                      (newtonian *
                       (1.0 - 4.0 * compactness - v_sq / c_sq));
            acc[i] -= velocity * (4.0 * newtonian * v_dot_n / c_sq);
        }
    }
}

} // namespace dynamics
} // namespace solar
