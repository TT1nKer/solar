#include "solar/spherical_harmonics.h"
#include "solar/ephemeris.h"
#include "solar/constants.h"
#include <cmath>

namespace solar {

SphericalHarmonics::SphericalHarmonics(std::vector<GravityField> fields)
    : fields_(std::move(fields)) {}

SphericalHarmonics SphericalHarmonics::solar_system_defaults(
    const std::vector<Body>& bodies, int max_degree)
{
    std::vector<GravityField> fields;

    auto add = [&](const char* name, double mu, double R_eq, int iau_id,
                   std::vector<double> Jn) {
        const Body* b = find_body(bodies, name);
        if (!b) return;
        int idx = static_cast<int>(b - bodies.data());
        int deg = std::min(max_degree, static_cast<int>(Jn.size()) + 1);
        Jn.resize(deg - 1);  // Jn[0]=J2, so size = max_degree - 1
        fields.push_back({idx, mu, R_eq, deg, std::move(Jn), iau_id});
    };

    using namespace constants;

    // Earth: J2-J6
    add("Earth", MU_EARTH, 6378.137, 3,
        {1.08263e-3, -2.532e-6, -1.620e-6, -2.273e-7, 5.407e-7});

    // Mars: J2-J4
    add("Mars", MU_MARS, 3396.2, 4,
        {1.96045e-3, 3.145e-5, -1.538e-5});

    // Jupiter: J2-J4
    add("Jupiter", MU_JUPITER, 71492.0, -1,
        {1.4736e-2, -5.80e-4, -3.40e-5});

    // Saturn: J2-J4
    add("Saturn", MU_SATURN, 60268.0, -1,
        {1.6298e-2, -1.00e-3, -1.00e-4});

    // Moon: J2 only (tesseral would need body-fixed frame)
    add("Moon", MU_MOON, 1738.1, 10,
        {2.033e-4});

    return SphericalHarmonics(std::move(fields));
}

// Zonal harmonics acceleration.
// Potential: V_Jn = mu*Jn*(Re/r)^n * Pn(u) / r  where u = z/r
// Acceleration: a = -grad(V_Jn)
// Derivation gives:
//   Q = mu * Jn * (Re/r)^n / r^3
//   F_xy = (n+1)*Pn(u) + u*Pn'(u)     [dimensionless]
//   a_x = Q * F_xy * x
//   a_y = Q * F_xy * y
//   a_z = Q * [(n+1)*Pn(u)*z - Pn'(u)*(r^2-z^2)/r]
Vec3 SphericalHarmonics::zonal_acceleration(const GravityField& field,
                                             const Vec3& r) const {
    double x = r.x, y = r.y, z = r.z;
    double r_mag = r.norm();

    double mu = field.mu;
    double Re = field.R_eq;
    double ax = 0.0, ay = 0.0, az = 0.0;

    double u = z / r_mag; // sin(latitude)
    double Re_over_r = Re / r_mag;

    // Legendre recursion: P_0=1, P_1=u
    double P_nm2 = 1.0;
    double P_nm1 = u;
    double Re_r_n = Re_over_r * Re_over_r; // starts at (Re/r)^2 for J2

    for (size_t i = 0; i < field.Jn.size(); ++i) {
        int n = static_cast<int>(i) + 2;
        double Jn = field.Jn[i];

        // P_n via recursion
        double Pn = (static_cast<double>(2*n - 1) * u * P_nm1
                    - static_cast<double>(n - 1) * P_nm2) / n;

        if (std::fabs(Jn) > 1e-30) {
            // dPn/du using: dPn = n*(P_{n-1} - u*Pn) / (1 - u^2)
            // Note: 1-u^2 = cos^2(lat)
            double cos2 = 1.0 - u * u;
            double dPn;
            if (cos2 > 1e-20) {
                dPn = static_cast<double>(n) * (P_nm1 - u * Pn) / cos2;
            } else {
                dPn = 0.5 * n * (n + 1) * (u > 0 ? 1.0 : -1.0);
            }

            double Q = mu * Jn * Re_r_n / (r_mag * r_mag * r_mag);
            double F_xy = static_cast<double>(n + 1) * Pn + u * dPn;

            ax += Q * F_xy * x;
            ay += Q * F_xy * y;
            az += Q * (static_cast<double>(n + 1) * Pn * z
                      - dPn * (r_mag * r_mag - z * z) / r_mag);
        }

        P_nm2 = P_nm1;
        P_nm1 = Pn;
        Re_r_n *= Re_over_r;
    }

    return {ax, ay, az};
}

void SphericalHarmonics::compute(
    const std::vector<Body>& bodies,
    double /*time*/,
    std::vector<Vec3>& acc) const
{
    size_t n = bodies.size();

    for (const auto& field : fields_) {
        int j = field.body_index;
        if (j < 0 || j >= static_cast<int>(n)) continue;

        for (size_t i = 0; i < n; ++i) {
            if (static_cast<int>(i) == j) continue;

            Vec3 r = bodies[i].state.pos - bodies[j].state.pos;
            Vec3 a = zonal_acceleration(field, r);
            acc[i] += a;
        }
    }
}

} // namespace solar
