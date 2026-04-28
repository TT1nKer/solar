#include "solar/force_model.h"

namespace solar {

double ForceModel::potential_energy(
    const std::vector<Body>& /*bodies*/,
    double /*time*/) const {
    return 0.0;
}

} // namespace solar
