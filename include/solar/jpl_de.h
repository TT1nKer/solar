#pragma once

#include "solar/body.h"
#include <vector>
#include <string>
#include <unordered_map>

namespace solar {

// Coefficient pointer info for one DE component (from Group 1050)
struct DEComponentInfo {
    int offset;           // 1-based offset into coefficient record
    int num_coeffs;       // Chebyshev coefficients per subinterval per coordinate
    int num_subintervals; // subintervals per record
};

// JPL Development Ephemeris reader (DE430, DE440, DE441, etc.)
// Reads the old JPL ASCII format (header.xxx + ascpXXXX.xxx)
class JPLEphemeris {
public:
    // Load from a single concatenated ASCII file (header + data)
    // or from separate header and data files
    bool load_ascii(const std::string& path);

    // Get heliocentric state of a body at Julian date (ecliptic J2000, km, km/s)
    // body_id: 1=Mercury, 2=Venus, 3=Earth, 4=Mars, 5=Jupiter,
    //          6=Saturn, 7=Uranus, 8=Neptune, 9=Pluto, 10=Moon(geocentric), 11=Sun
    bool get_state(int body_id, double jd, State& state) const;

    bool is_loaded() const { return loaded_; }
    bool covers(double jd) const { return loaded_ && jd >= jd_start_ && jd <= jd_end_; }

    double jd_start() const { return jd_start_; }
    double jd_end() const { return jd_end_; }
    double record_span() const { return record_span_; }
    double emrat() const { return emrat_; }
    double au() const { return au_km_; }
    int num_records() const { return static_cast<int>(records_.size()); }

    // Get a named constant from the file
    double get_constant(const std::string& name) const;

private:
    bool loaded_ = false;
    double jd_start_ = 0, jd_end_ = 0, record_span_ = 0;
    int ncoeff_ = 0;  // coefficients per record
    double emrat_ = 81.3005690699;  // Earth/Moon mass ratio (overwritten from file)
    double au_km_ = 149597870.7;    // AU in km (overwritten from file)

    // 15 component slots (13 used + 2 padding in some DE versions)
    static constexpr int MAX_COMPONENTS = 15;
    DEComponentInfo components_[MAX_COMPONENTS] = {};
    int num_components_ = 13;

    // Constants from the file
    std::vector<std::string> const_names_;
    std::vector<double> const_values_;
    std::unordered_map<std::string, double> constants_;

    // Coefficient records: records_[i] holds ncoeff_ doubles
    std::vector<std::vector<double>> records_;
    double data_jd_start_ = 0; // JD of first data record

    // Parse header groups
    bool parse_header(std::istream& in);
    // Parse data records
    bool parse_data(std::istream& in);

    // Evaluate Chebyshev for a component, returning position and velocity
    void eval_chebyshev(int component, double jd, Vec3& pos, Vec3& vel) const;

    // Get raw barycentric state for a component (ICRF equatorial)
    bool get_raw_state(int component, double jd, Vec3& pos, Vec3& vel) const;

    // Transform ICRF equatorial to ecliptic J2000
    static void equatorial_to_ecliptic(Vec3& pos, Vec3& vel);
};

// Map body name to DE body ID (1-11)
int body_name_to_de_id(const std::string& name);

// Global DE ephemeris instance
JPLEphemeris* global_de_ephemeris();
void set_global_de_ephemeris(JPLEphemeris* eph);

} // namespace solar
