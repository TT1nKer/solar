#include "solar/jpl_de.h"
#include "solar/ephemeris.h"
#include "solar/constants.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace solar;
using namespace solar::constants;

static int pass = 0, fail = 0;

static void check(const char* name, double got, double expected, double tol) {
    double err = std::fabs(got - expected);
    if (err < tol) {
        std::cout << "  PASS: " << name << " = " << got << " (err=" << err << ")\n";
        pass++;
    } else {
        std::cout << "  FAIL: " << name << " = " << got
                  << " (expected " << expected << ", err=" << err << ")\n";
        fail++;
    }
}

int main(int argc, char* argv[]) {
    std::string de_path = "data/de440.asc";
    if (argc > 1) de_path = argv[1];

    JPLEphemeris de;
    if (!de.load_ascii(de_path)) {
        if (argc == 1) {
            std::cout << "SKIP: optional DE440 fixture not found at "
                      << de_path << "\n";
            return 0;
        }
        std::cerr << "Cannot load DE file: " << de_path << "\n";
        return 1;
    }
    set_global_de_ephemeris(&de);

    std::cout << std::scientific << std::setprecision(10);
    std::cout << "DE loaded: " << de.num_records() << " records\n";
    std::cout << "EMRAT = " << de.emrat() << "\n";
    std::cout << "AU = " << de.au() << " km\n\n";

    // ============================================
    // Test 1: Earth at 2000-01-01 00:00 (JD 2451544.5)
    // JPL Horizons reference (ecliptic J2000, heliocentric):
    //   X = -2.521092855899356e+07  Y = 1.449279195838006e+08  Z = -6.164165719002485e+02
    //   VX = -2.983983333677879e+01  VY = -5.207633902410673  VZ = 6.168441184239981e-05
    // ============================================
    std::cout << "=== Earth at JD 2451544.5 (vs JPL Horizons) ===\n";
    {
        double jd = 2451544.5;
        State s;
        de.get_state(3, jd, s);

        check("Earth X", s.pos.x, -2.521092855899356e+07, 10.0);  // 10 km tolerance
        check("Earth Y", s.pos.y,  1.449279195838006e+08, 10.0);
        check("Earth Z", s.pos.z, -6.164165719002485e+02, 1.0);

        check("Earth VX", s.vel.x, -2.983983333677879e+01, 0.01); // 10 m/s tolerance
        check("Earth VY", s.vel.y, -5.207633902410673e+00, 0.01);
        check("Earth VZ", s.vel.z,  6.168441184239981e-05, 0.001);

        double dist_err = std::sqrt(
            (s.pos.x - (-2.521092855899356e+07)) * (s.pos.x - (-2.521092855899356e+07)) +
            (s.pos.y - 1.449279195838006e+08) * (s.pos.y - 1.449279195838006e+08) +
            (s.pos.z - (-6.164165719002485e+02)) * (s.pos.z - (-6.164165719002485e+02)));
        std::cout << "  Total position error: " << std::fixed << std::setprecision(3) << dist_err << " km\n\n";
    }

    // ============================================
    // Test 2: Mars position
    // ============================================
    std::cout << "=== Mars at JD 2451544.5 ===\n";
    {
        State s;
        de.get_state(4, 2451544.5, s);
        double dist = s.pos.norm();
        std::cout << "  Mars distance from Sun: " << std::scientific << dist << " km ("
                  << std::fixed << std::setprecision(4) << dist / AU << " AU)\n";
        // Mars should be ~1.4-1.7 AU
        check("Mars dist (AU)", dist / AU, 1.39, 0.3); // wide tolerance for orbit position
    }

    // ============================================
    // Test 3: Fallback — without DE, use Keplerian
    // ============================================
    std::cout << "\n=== Fallback test (Keplerian) ===\n";
    {
        set_global_de_ephemeris(nullptr);
        auto bodies = load_solar_system();
        const Body* earth = find_body(bodies, "Earth");
        State s_kep = compute_position(*earth, 2451544.5);
        std::cout << "  Keplerian Earth X = " << std::scientific << s_kep.pos.x << "\n";

        set_global_de_ephemeris(&de);
        State s_de = compute_position(*earth, 2451544.5);
        std::cout << "  DE440 Earth X    = " << s_de.pos.x << "\n";

        double diff = std::fabs(s_de.pos.x - s_kep.pos.x);
        std::cout << "  Difference: " << std::fixed << std::setprecision(0) << diff << " km\n";
        // They should be different (DE is more accurate)
        if (diff > 100) {
            std::cout << "  PASS: DE and Keplerian give different results (expected)\n";
            pass++;
        } else {
            std::cout << "  FAIL: DE and Keplerian should differ\n";
            fail++;
        }
    }

    std::cout << "\n=== Results: " << pass << " passed, " << fail << " failed ===\n";
    return fail > 0 ? 1 : 0;
}
