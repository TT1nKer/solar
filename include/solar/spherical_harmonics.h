#pragma once

#include "solar/force_model.h"
#include <vector>

namespace solar {

// Gravity field specification for one body
struct GravityField {
    int body_index;         // index in bodies vector
    double mu;              // GM (km^3/s^2)
    double R_eq;            // reference equatorial radius (km)
    int max_degree;         // highest degree (N in Jn)
    std::vector<double> Jn; // zonal harmonics: Jn[0]=J2, Jn[1]=J3, ...
    int iau_body_id;        // for body-fixed frame (3=Earth, 4=Mars, 10=Moon)
};

// Spherical harmonics gravity perturbation (zonal Jn expansion)
// Generalizes J2Perturbation to arbitrary degree
class SphericalHarmonics : public ForceModel {
public:
    explicit SphericalHarmonics(std::vector<GravityField> fields);

    // Factory: standard solar system gravity fields
    static SphericalHarmonics solar_system_defaults(
        const std::vector<Body>& bodies, int max_degree = 6);

    std::string name() const override { return "spherical_harmonics"; }

    void compute(
        const std::vector<Body>& bodies,
        double time,
        std::vector<Vec3>& acc) const override;

private:
    std::vector<GravityField> fields_;

    // Compute zonal harmonics acceleration on a body
    // r: position relative to the oblate body (inertial frame, pole=z)
    Vec3 zonal_acceleration(const GravityField& field, const Vec3& r) const;
};

} // namespace solar
