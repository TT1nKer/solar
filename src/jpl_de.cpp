#include "solar/jpl_de.h"
#include "solar/constants.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <cctype>

namespace solar {

// Global instance pointer
static JPLEphemeris* g_de_ephemeris = nullptr;

JPLEphemeris* global_de_ephemeris() { return g_de_ephemeris; }
void set_global_de_ephemeris(JPLEphemeris* eph) { g_de_ephemeris = eph; }

// Convert Fortran D-notation to standard E-notation
static double parse_d_float(const std::string& s) {
    std::string t = s;
    for (auto& c : t) {
        if (c == 'D' || c == 'd') c = 'E';
    }
    return std::stod(t);
}

// Trim whitespace
static std::string trim(const std::string& s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    size_t end = s.find_last_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    return s.substr(start, end - start + 1);
}

bool JPLEphemeris::load_ascii(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "DE: cannot open " << path << "\n";
        return false;
    }

    if (!parse_header(file)) {
        std::cerr << "DE: failed to parse header\n";
        return false;
    }

    if (!parse_data(file)) {
        std::cerr << "DE: failed to parse data records\n";
        return false;
    }

    // Extract key constants
    if (constants_.count("EMRAT")) emrat_ = constants_["EMRAT"];
    if (constants_.count("AU")) au_km_ = constants_["AU"];

    loaded_ = true;
    return true;
}

bool JPLEphemeris::parse_header(std::istream& in) {
    std::string line;
    int ncoeff_from_header = 0;

    // Read KSIZE/NCOEFF line
    while (std::getline(in, line)) {
        if (line.find("NCOEFF") != std::string::npos) {
            std::istringstream iss(line);
            std::string kw;
            int ksize;
            iss >> kw >> ksize >> kw >> ncoeff_from_header;
            ncoeff_ = ncoeff_from_header;
            break;
        }
    }

    // Parse groups
    while (std::getline(in, line)) {
        std::string trimmed = trim(line);

        if (trimmed.find("GROUP") != std::string::npos && trimmed.find("1030") != std::string::npos) {
            // Group 1030: start JD, end JD, record span
            std::getline(in, line); // blank
            std::getline(in, line);
            std::istringstream iss(line);
            std::string s1, s2, s3;
            iss >> s1 >> s2 >> s3;
            jd_start_ = parse_d_float(s1);
            jd_end_ = parse_d_float(s2);
            record_span_ = parse_d_float(s3);
        }
        else if (trimmed.find("GROUP") != std::string::npos && trimmed.find("1040") != std::string::npos) {
            // Group 1040: constant names
            std::getline(in, line); // blank
            std::getline(in, line);
            int nconst = std::stoi(trim(line));
            const_names_.clear();
            while (static_cast<int>(const_names_.size()) < nconst) {
                std::getline(in, line);
                std::istringstream iss(line);
                std::string name;
                while (iss >> name) {
                    const_names_.push_back(name);
                }
            }
        }
        else if (trimmed.find("GROUP") != std::string::npos && trimmed.find("1041") != std::string::npos) {
            // Group 1041: constant values
            std::getline(in, line); // blank
            std::getline(in, line);
            int nconst = std::stoi(trim(line));
            const_values_.clear();
            while (static_cast<int>(const_values_.size()) < nconst) {
                std::getline(in, line);
                std::istringstream iss(line);
                std::string val;
                while (iss >> val) {
                    const_values_.push_back(parse_d_float(val));
                }
            }
            // Build constants map
            for (size_t i = 0; i < const_names_.size() && i < const_values_.size(); ++i) {
                constants_[const_names_[i]] = const_values_[i];
            }
        }
        else if (trimmed.find("GROUP") != std::string::npos && trimmed.find("1050") != std::string::npos) {
            // Group 1050: coefficient pointers (3 rows of integers)
            std::getline(in, line); // blank

            // Row 1: offsets
            std::vector<int> offsets, ncoeffs, nsubints;
            std::getline(in, line);
            { std::istringstream iss(line); int v; while (iss >> v) offsets.push_back(v); }

            // Row 2: num coefficients
            std::getline(in, line);
            { std::istringstream iss(line); int v; while (iss >> v) ncoeffs.push_back(v); }

            // Row 3: num subintervals
            std::getline(in, line);
            { std::istringstream iss(line); int v; while (iss >> v) nsubints.push_back(v); }

            int nc = std::min({static_cast<int>(offsets.size()),
                               static_cast<int>(ncoeffs.size()),
                               static_cast<int>(nsubints.size()),
                               MAX_COMPONENTS});
            num_components_ = nc;
            for (int i = 0; i < nc; ++i) {
                components_[i].offset = offsets[i];
                components_[i].num_coeffs = ncoeffs[i];
                components_[i].num_subintervals = nsubints[i];
            }
        }
        else if (trimmed.find("GROUP") != std::string::npos && trimmed.find("1070") != std::string::npos) {
            // End of header, data follows
            break;
        }
    }

    return ncoeff_ > 0 && record_span_ > 0;
}

bool JPLEphemeris::parse_data(std::istream& in) {
    std::string line;
    records_.clear();
    data_jd_start_ = 0;

    while (std::getline(in, line)) {
        std::string trimmed = trim(line);
        if (trimmed.empty()) continue;

        // Each record starts with: record_number  ncoeff
        // Then the coefficients follow
        std::istringstream first_iss(trimmed);
        int rec_num, nc;
        if (!(first_iss >> rec_num >> nc)) continue;
        if (nc != ncoeff_) continue; // sanity check

        // Read ncoeff_ doubles
        std::vector<double> coeffs;
        coeffs.reserve(ncoeff_);

        while (static_cast<int>(coeffs.size()) < ncoeff_) {
            if (!std::getline(in, line)) break;
            std::istringstream iss(line);
            std::string token;
            while (iss >> token && static_cast<int>(coeffs.size()) < ncoeff_) {
                coeffs.push_back(parse_d_float(token));
            }
        }

        if (static_cast<int>(coeffs.size()) == ncoeff_) {
            if (records_.empty()) {
                data_jd_start_ = coeffs[0]; // first two doubles are jd_start, jd_end of record
            }
            records_.push_back(std::move(coeffs));
        }
    }

    return !records_.empty();
}

double JPLEphemeris::get_constant(const std::string& name) const {
    auto it = constants_.find(name);
    if (it != constants_.end()) return it->second;
    return 0.0;
}

void JPLEphemeris::eval_chebyshev(int component, double jd, Vec3& pos, Vec3& vel) const {
    const auto& ci = components_[component];
    if (ci.num_coeffs == 0 || ci.num_subintervals == 0) {
        pos = vel = {0, 0, 0};
        return;
    }

    // Find record
    int rec_idx = static_cast<int>((jd - data_jd_start_) / record_span_);
    if (rec_idx < 0) rec_idx = 0;
    if (rec_idx >= static_cast<int>(records_.size())) rec_idx = static_cast<int>(records_.size()) - 1;

    const auto& rec = records_[rec_idx];
    double rec_start = rec[0];
    double rec_end = rec[1];

    // Find subinterval
    double subinterval_span = (rec_end - rec_start) / ci.num_subintervals;
    int subint = static_cast<int>((jd - rec_start) / subinterval_span);
    if (subint >= ci.num_subintervals) subint = ci.num_subintervals - 1;
    if (subint < 0) subint = 0;

    double subint_start = rec_start + subint * subinterval_span;

    // Normalize time to [-1, 1] within subinterval
    double tau = 2.0 * (jd - subint_start) / subinterval_span - 1.0;

    // Pointer into coefficients (0-based, offset is 1-based)
    int base = ci.offset - 1 + subint * ci.num_coeffs * 3;
    int nc = ci.num_coeffs;

    double p[3], v[3];

    for (int coord = 0; coord < 3; ++coord) {
        int idx = base + coord * nc;

        // --- Position: Clenshaw on T_n(tau) ---
        // f(tau) = sum_{i=0}^{N-1} c_i * T_i(tau)
        double w0 = 0.0, w1 = 0.0;
        for (int k = nc - 1; k >= 1; --k) {
            double w_new = 2.0 * tau * w0 - w1 + rec[idx + k];
            w1 = w0;
            w0 = w_new;
        }
        p[coord] = tau * w0 - w1 + rec[idx];

        // --- Velocity: Clenshaw on U_{n-1}(tau) for derivative ---
        // f'(tau) = sum_{i=1}^{N-1} i * c_i * U_{i-1}(tau)
        //         = sum_{j=0}^{N-2} (j+1) * c_{j+1} * U_j(tau)
        // Clenshaw for U: b0 = b1 = 0
        //   b_j = 2*tau*b_{j+1} - b_{j+2} + a_j  where a_j = (j+1)*c_{j+1}
        double b0 = 0.0, b1 = 0.0;
        for (int j = nc - 2; j >= 1; --j) {
            double b_new = 2.0 * tau * b0 - b1 + static_cast<double>(j + 1) * rec[idx + j + 1];
            b1 = b0;
            b0 = b_new;
        }
        // j=0 term: a_0 = 1*c_1, U_0 = 1
        v[coord] = 2.0 * tau * b0 - b1 + rec[idx + 1];

        // Scale: d/dtau → d/dt (days to seconds)
        v[coord] *= 2.0 / (subinterval_span * constants::DAY);
    }

    pos = {p[0], p[1], p[2]};
    vel = {v[0], v[1], v[2]};
}

bool JPLEphemeris::get_raw_state(int component, double jd, Vec3& pos, Vec3& vel) const {
    if (component < 0 || component >= num_components_) return false;
    if (components_[component].num_coeffs == 0) return false;
    eval_chebyshev(component, jd, pos, vel);
    return true;
}

void JPLEphemeris::equatorial_to_ecliptic(Vec3& pos, Vec3& vel) {
    double cos_e = std::cos(constants::OBLIQUITY_J2000);
    double sin_e = std::sin(constants::OBLIQUITY_J2000);

    // Rotate around X-axis by obliquity
    double py = pos.y, pz = pos.z;
    pos.y =  cos_e * py + sin_e * pz;
    pos.z = -sin_e * py + cos_e * pz;

    double vy = vel.y, vz = vel.z;
    vel.y =  cos_e * vy + sin_e * vz;
    vel.z = -sin_e * vy + cos_e * vz;
}

bool JPLEphemeris::get_state(int body_id, double jd, State& state) const {
    if (!loaded_ || !covers(jd)) return false;

    // Component indices (0-based): 0=Mercury..8=Pluto, 9=Moon(geo), 10=Sun
    // body_id (1-based): 1=Mercury..9=Pluto, 10=Moon, 11=Sun

    // Get Sun's barycentric state (needed for heliocentric conversion)
    Vec3 sun_pos, sun_vel;
    if (!get_raw_state(10, jd, sun_pos, sun_vel)) return false; // component 10 = Sun

    if (body_id == 11) {
        // Sun: heliocentric = (0,0,0) by definition
        state.pos = {0, 0, 0};
        state.vel = {0, 0, 0};
        return true;
    }

    if (body_id == 3) {
        // Earth: compute from EMB (component 2) and Moon (component 9)
        Vec3 emb_pos, emb_vel, moon_pos, moon_vel;
        if (!get_raw_state(2, jd, emb_pos, emb_vel)) return false;
        if (!get_raw_state(9, jd, moon_pos, moon_vel)) return false;

        double factor = 1.0 / (1.0 + emrat_);
        // Earth = EMB - Moon_geo * factor (Moon_geo is geocentric)
        Vec3 earth_pos = emb_pos - moon_pos * factor;
        Vec3 earth_vel = emb_vel - moon_vel * factor;

        // Barycentric to heliocentric
        state.pos = earth_pos - sun_pos;
        state.vel = earth_vel - sun_vel;
        equatorial_to_ecliptic(state.pos, state.vel);
        return true;
    }

    if (body_id == 10) {
        // Moon (geocentric): get geocentric moon, then convert to heliocentric
        // by adding Earth's heliocentric position
        Vec3 moon_geo_pos, moon_geo_vel;
        if (!get_raw_state(9, jd, moon_geo_pos, moon_geo_vel)) return false;

        // We need this in ecliptic for the caller
        equatorial_to_ecliptic(moon_geo_pos, moon_geo_vel);
        state.pos = moon_geo_pos;
        state.vel = moon_geo_vel;
        return true;
    }

    // Planets 1-9 (Mercury through Pluto): components 0-8
    int component = body_id - 1;
    Vec3 body_pos, body_vel;
    if (!get_raw_state(component, jd, body_pos, body_vel)) return false;

    // Special case: body 3 (EMB) handled above; bodies 1,2,4-9 are barycentric
    // Convert to heliocentric
    state.pos = body_pos - sun_pos;
    state.vel = body_vel - sun_vel;
    equatorial_to_ecliptic(state.pos, state.vel);
    return true;
}

int body_name_to_de_id(const std::string& name) {
    std::string lower = name;
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    if (lower == "mercury")  return 1;
    if (lower == "venus")    return 2;
    if (lower == "earth")    return 3;
    if (lower == "mars")     return 4;
    if (lower == "jupiter")  return 5;
    if (lower == "saturn")   return 6;
    if (lower == "uranus")   return 7;
    if (lower == "neptune")  return 8;
    if (lower == "pluto")    return 9;
    if (lower == "moon")     return 10;
    if (lower == "sun")      return 11;
    return -1;
}

} // namespace solar
