#include "solar/dynamics/pn_gravity.h"

#include "solar/constants.h"
#include "solar/dynamics/barnes_hut_octree.h"

#include <cmath>
#include <cstddef>
#include <vector>

namespace solar {
namespace dynamics {

namespace {

double smoothstep(double x) {
    return x * x * (3.0 - 2.0 * x);
}

} // namespace

PostNewtonianGravity::PostNewtonianGravity(double speed_of_light_km_s)
    : PostNewtonianGravity(Config{}, speed_of_light_km_s) {}

PostNewtonianGravity::PostNewtonianGravity(Config config,
                                           double speed_of_light_km_s)
    : config_(config), speed_of_light_km_s_(speed_of_light_km_s) {}

void PostNewtonianGravity::compute(
    const std::vector<Body>& bodies,
    double /*time*/,
    std::vector<Vec3>& acc) const {
    if (bodies.empty()) return;
    const double c = speed_of_light_km_s_ > 0.0
        ? speed_of_light_km_s_
        : constants::C_LIGHT;  // km/s
    const double c_sq = c * c;
    const octree_detail::Tree tree(bodies);
    const double softening_sq = config_.softening_km * config_.softening_km;
    const double theta = config_.opening_angle;

    for (std::size_t i = 0; i < bodies.size(); ++i) {
        // Newtonian field and potential of the OTHER bodies, from the
        // same monopole walk as BarnesHutGravity (same self-exclusion,
        // same opening-angle policy).
        Vec3 field{};
        double potential = 0.0;
        octree_detail::monopole_field(tree, bodies, i, softening_sq, theta,
                                      field, potential);

        const Vec3& velocity = bodies[i].state.vel;
        const double epsilon = -potential / c_sq;  // G M_enc / (r c^2)
        double fade = 1.0;
        if (config_.fade_hi > 0.0) {
            fade = epsilon < config_.fade_hi
                ? 1.0 - smoothstep(epsilon / config_.fade_hi)
                : 0.0;
        }
        const double v_sq = velocity.norm_sq();
        const double g_dot_v = field.dot(velocity);
        acc[i] += (field * (1.0 + 4.0 * potential / c_sq - v_sq / c_sq) -
                   velocity * (4.0 * g_dot_v / c_sq)) * fade;
    }
}

} // namespace dynamics
} // namespace solar