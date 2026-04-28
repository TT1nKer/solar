#pragma once

#include "solar/body.h"
#include <vector>

namespace solar {

// Append moon bodies to an existing bodies vector.
// Planets must already be loaded (sets parent_index by scanning for parent name).
// Moon elements are relative to their parent planet.
void load_moons(std::vector<Body>& bodies);

} // namespace solar
