#include "solar/body.h"
#include "solar/constants.h"
#include "solar/dynamics/ltb_collapse.h"
#include "solar/dynamics/pn_collapse.h"
#include "solar/dynamics/pn_gravity.h"
#include "solar/dynamics/turbulent_cloud.h"
#include "solar/nbody.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

// Traceable collapse movie: renders the three verified stages of the
// Newton -> 1PN -> GR chain as PPM frames with per-frame captions of
// the measured physical quantities, plus a JSON sidecar that ties every
// frame to its numbers.
//
//   stage 1 (nebula): turbulent 1e4 M_sun / 1 pc cloud + embedded
//                     coherent core, PnCollapseForce, 0.55 t_sing
//   stage 2 (compact): 200 G M / c^2 ball, field-based 1PN collapse
//   stage 3 (GR): Oppenheimer-Snyder surface, horizon ring, and the
//                     observer-time redshift / luminosity curves
//
// Usage: test_collapse_movie <output-directory>
// Assemble: ffmpeg -y -framerate 30 -i DIR/nebula-%03d.ppm -c:v libx264
// -pix_fmt yuv420p nebula.mp4 (repeat per stage, then concat).

using solar::Body;
using solar::NBodySim;
using solar::Vec3;
using solar::dynamics::LTBCollapse;
using solar::dynamics::PnCollapseForce;
using solar::dynamics::PostNewtonianGravity;
using solar::dynamics::TurbulentCloudGenerator;

namespace {

int failures = 0;

void check(const std::string& name, bool condition) {
    std::cout << (condition ? "  PASS: " : "  FAIL: ") << name << '\n';
    if (!condition) ++failures;
}

constexpr double pi = 3.14159265358979323846;
constexpr double parsec_km = 3.085677581e13;
constexpr double solar_mass_kg = 1.98892e30;

const int width = 960;
const int height = 540;

// Minimal 5x7 bitmap font for the caption strip.
const char* const glyphs =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789 .:|=<>-()/_,~e"
    "abcdfgilnoprstuy+";
const unsigned char font[68][7] = {
    {0x0E,0x11,0x11,0x1F,0x11,0x11,0x11}, // A
    {0x1E,0x11,0x11,0x1E,0x11,0x11,0x1E}, // B
    {0x0E,0x11,0x10,0x10,0x10,0x11,0x0E}, // C
    {0x1E,0x11,0x11,0x11,0x11,0x11,0x1E}, // D
    {0x1F,0x10,0x10,0x1E,0x10,0x10,0x1F}, // E
    {0x1F,0x10,0x10,0x1E,0x10,0x10,0x10}, // F
    {0x0E,0x11,0x10,0x17,0x11,0x11,0x0F}, // G
    {0x11,0x11,0x11,0x1F,0x11,0x11,0x11}, // H
    {0x1F,0x04,0x04,0x04,0x04,0x04,0x1F}, // I
    {0x07,0x02,0x02,0x02,0x02,0x12,0x0C}, // J
    {0x11,0x12,0x14,0x18,0x14,0x12,0x11}, // K
    {0x10,0x10,0x10,0x10,0x10,0x10,0x1F}, // L
    {0x11,0x1B,0x15,0x15,0x11,0x11,0x11}, // M
    {0x11,0x19,0x15,0x13,0x11,0x11,0x11}, // N
    {0x0E,0x11,0x11,0x11,0x11,0x11,0x0E}, // O
    {0x1E,0x11,0x11,0x1E,0x10,0x10,0x10}, // P
    {0x0E,0x11,0x11,0x11,0x15,0x12,0x0D}, // Q
    {0x1E,0x11,0x11,0x1E,0x14,0x12,0x11}, // R
    {0x0F,0x10,0x10,0x0E,0x01,0x01,0x1E}, // S
    {0x1F,0x04,0x04,0x04,0x04,0x04,0x04}, // T
    {0x11,0x11,0x11,0x11,0x11,0x11,0x0E}, // U
    {0x11,0x11,0x11,0x11,0x11,0x0A,0x04}, // V
    {0x11,0x11,0x11,0x15,0x15,0x15,0x0A}, // W
    {0x11,0x11,0x0A,0x04,0x0A,0x11,0x11}, // X
    {0x11,0x11,0x0A,0x04,0x04,0x04,0x04}, // Y
    {0x1F,0x01,0x02,0x04,0x08,0x10,0x1F}, // Z
    {0x0E,0x11,0x13,0x15,0x19,0x11,0x0E}, // 0
    {0x04,0x0C,0x04,0x04,0x04,0x04,0x0E}, // 1
    {0x0E,0x11,0x01,0x02,0x04,0x08,0x1F}, // 2
    {0x1F,0x02,0x04,0x02,0x01,0x11,0x0E}, // 3
    {0x02,0x06,0x0A,0x12,0x1F,0x02,0x02}, // 4
    {0x1F,0x10,0x1E,0x01,0x01,0x11,0x0E}, // 5
    {0x06,0x08,0x10,0x1E,0x11,0x11,0x0E}, // 6
    {0x1F,0x01,0x02,0x04,0x08,0x08,0x08}, // 7
    {0x0E,0x11,0x11,0x0E,0x11,0x11,0x0E}, // 8
    {0x0E,0x11,0x11,0x0F,0x01,0x02,0x0C}, // 9
    {0x00,0x00,0x00,0x00,0x00,0x00,0x00}, // space
    {0x00,0x00,0x00,0x00,0x00,0x04,0x04}, // .
    {0x00,0x04,0x04,0x1F,0x04,0x04,0x00}, // :
    {0x04,0x04,0x04,0x04,0x04,0x04,0x04}, // |
    {0x00,0x00,0x00,0x1F,0x00,0x00,0x00}, // =
    {0x00,0x04,0x02,0x01,0x02,0x04,0x00}, // <
    {0x00,0x01,0x02,0x04,0x02,0x01,0x00}, // >
    {0x00,0x00,0x00,0x00,0x00,0x00,0x0F}, // -
    {0x00,0x00,0x00,0x1F,0x0A,0x04,0x00}, // (
    {0x00,0x00,0x00,0x00,0x00,0x00,0x00}, // )
    {0x01,0x01,0x02,0x04,0x08,0x08,0x08}, // /
    {0x00,0x00,0x00,0x00,0x00,0x00,0x00}, // ,
    {0x00,0x00,0x00,0x00,0x00,0x00,0x1F}, // _
    {0x00,0x0A,0x15,0x0A,0x00,0x00,0x00}, // ~
    {0x0E,0x11,0x1F,0x10,0x10,0x11,0x0E}, // e
    {0x00,0x00,0x0E,0x01,0x0F,0x11,0x0F}, // a
    {0x10,0x10,0x1E,0x11,0x11,0x11,0x1E}, // b
    {0x00,0x00,0x0E,0x10,0x10,0x10,0x0E}, // c
    {0x01,0x01,0x0F,0x11,0x11,0x11,0x0F}, // d
    {0x00,0x00,0x0E,0x11,0x1E,0x10,0x0E}, // f
    {0x00,0x00,0x0E,0x11,0x11,0x0F,0x01}, // g
    {0x10,0x10,0x1E,0x11,0x11,0x11,0x11}, // i
    {0x10,0x10,0x10,0x10,0x10,0x10,0x1C}, // l
    {0x00,0x00,0x16,0x19,0x11,0x11,0x11}, // n
    {0x00,0x00,0x0E,0x11,0x11,0x11,0x0E}, // o
    {0x00,0x00,0x1E,0x11,0x1E,0x10,0x10}, // p
    {0x00,0x00,0x16,0x1A,0x0C,0x04,0x04}, // r
    {0x00,0x00,0x0E,0x10,0x0E,0x01,0x1E}, // s
    {0x00,0x04,0x1F,0x04,0x04,0x04,0x04}, // t
    {0x00,0x00,0x11,0x11,0x11,0x11,0x0E}, // u
    {0x00,0x00,0x11,0x11,0x0A,0x04,0x04}, // y
    {0x00,0x04,0x04,0x1F,0x04,0x04,0x00}, // +
};

int glyph_index(char ch) {
    const char* found = std::strchr(glyphs, ch);
    const int index = found ? static_cast<int>(found - glyphs) : 27;
    return (index >= 0 && index < 68) ? index : 27;
}

struct Frame {
    std::vector<unsigned char> rgb;
    Frame() : rgb(width * height * 3, 0) {}
    void blend(int x, int y, double r, double g, double b, double a) {
        if (x < 0 || y < 0 || x >= width || y >= height) return;
        const std::size_t i = (static_cast<std::size_t>(y) * width + x) * 3;
        rgb[i] = static_cast<unsigned char>(
            std::min(255.0, rgb[i] + r * a));
        rgb[i + 1] = static_cast<unsigned char>(
            std::min(255.0, rgb[i + 1] + g * a));
        rgb[i + 2] = static_cast<unsigned char>(
            std::min(255.0, rgb[i + 2] + b * a));
    }
    // Gaussian splat for a particle: soft, additive, brighter than a
    // single pixel so sparse clouds read as nebulosity.
    void splat(double cx, double cy, int r, int g, int b, double a) {
        const int x0 = static_cast<int>(cx) - 2;
        const int y0 = static_cast<int>(cy) - 2;
        for (int dy = 0; dy <= 4; ++dy) {
            for (int dx = 0; dx <= 4; ++dx) {
                const double d2 = (dx - 2) * (dx - 2) +
                                  (dy - 2) * (dy - 2);
                blend(x0 + dx, y0 + dy, r, g, b,
                      a * std::exp(-0.5 * d2));
            }
        }
    }
    // Deep-space gradient background with fixed stars.
    void background(const std::vector<std::array<int, 3>>& stars) {
        for (int y = 0; y < height; ++y) {
            const double shade = 0.55 + 0.45 * y / height;
            const std::size_t row = static_cast<std::size_t>(y) * width * 3;
            for (int x = 0; x < width; ++x) {
                rgb[row + static_cast<std::size_t>(x) * 3] =
                    static_cast<unsigned char>(8 * shade);
                rgb[row + static_cast<std::size_t>(x) * 3 + 1] =
                    static_cast<unsigned char>(12 * shade);
                rgb[row + static_cast<std::size_t>(x) * 3 + 2] =
                    static_cast<unsigned char>(26 * shade);
            }
        }
        for (const auto& star : stars) {
            blend(star[0], star[1], star[2], star[2], star[2] * 1.2,
                  0.8);
        }
    }
    void ring(double cx, double cy, double radius_px, int r, int g, int b) {
        const int r0 = static_cast<int>(radius_px);
        for (int dy = -r0; dy <= r0; ++dy) {
            for (int dx = -r0; dx <= r0; ++dx) {
                const double d = std::sqrt(dx * dx + dy * dy);
                if (std::abs(d - radius_px) < 1.5) {
                    blend(static_cast<int>(cx) + dx,
                          static_cast<int>(cy) + dy, r, g, b, 1.0);
                }
            }
        }
    }
    // Soft glowing ring: wide faint halo + bright core line.
    void glow_ring(double cx, double cy, double radius_px, int r, int g,
                   int b) {
        const int r0 = static_cast<int>(radius_px);
        for (int dy = -r0 - 4; dy <= r0 + 4; ++dy) {
            for (int dx = -r0 - 4; dx <= r0 + 4; ++dx) {
                const double d = std::sqrt(dx * dx + dy * dy);
                const double gap = std::abs(d - radius_px);
                if (gap < 4.0) {
                    blend(static_cast<int>(cx) + dx,
                          static_cast<int>(cy) + dy, r, g, b,
                          std::exp(-0.35 * gap * gap));
                }
            }
        }
    }
    void fill_rect(int x0, int y0, int x1, int y1, int r, int g, int b,
                   double a) {
        for (int y = y0; y <= y1; ++y) {
            for (int x = x0; x <= x1; ++x) {
                blend(x, y, r, g, b, a);
            }
        }
    }
    void text(int x, int y, const std::string& caption, int scale = 1) {
        for (const char ch : caption) {
            const int glyph = glyph_index(ch);
            for (int row = 0; row < 7; ++row) {
                for (int col = 0; col < 5; ++col) {
                    if (font[glyph][row] & (1 << (4 - col))) {
                        fill_rect(x + col * scale, y + row * scale,
                                  x + col * scale + scale - 1,
                                  y + row * scale + scale - 1,
                                  255, 255, 255, 0.95);
                    }
                }
            }
            x += 6 * scale;
        }
    }
};

void write_ppm(const std::string& path, const Frame& frame) {
    std::ofstream out(path, std::ios::binary);
    out << "P6\n" << width << " " << height << "\n255\n";
    out.write(reinterpret_cast<const char*>(frame.rgb.data()),
              static_cast<std::streamsize>(frame.rgb.size()));
}

struct Sidecar {
    std::vector<std::string> entries;
    void add(const std::string& stage, int index, double t, double r,
             double eps, double jm, double energy,
             const std::string& caption) {
        char buffer[600];
        std::snprintf(buffer, sizeof(buffer),
            "{\"stage\":\"%s\",\"frame\":%d,\"t\":%.6e,"
            "\"r_surface\":%.6e,\"eps\":%.6e,\"jm\":%.6e,"
            "\"energy\":%.6e,\"caption\":\"%s\"}",
            stage.c_str(), index, t, r, eps, jm, energy,
            caption.c_str());
        entries.push_back(buffer);
    }
    void write(const std::string& path) const {
        std::ofstream out(path);
        out << "[\n";
        for (std::size_t i = 0; i < entries.size(); ++i) {
            out << entries[i]
                << (i + 1 < entries.size() ? ",\n" : "\n");
        }
        out << "]\n";
    }
};

void project(const Vec3& pos, double box, int& x, int& y) {
    x = static_cast<int>((pos.x / box + 0.5) * width);
    y = static_cast<int>((0.5 - pos.z / box) * height);
}

std::vector<Body> nebula_ic(double mass_solar, std::size_t n_ambient,
                            std::size_t n_clump, double r_c_pc,
                            double clump_mass_solar) {
    TurbulentCloudGenerator::Config config;
    config.radius_pc = 1.0;
    config.mass_solar = mass_solar;
    config.particle_count = n_ambient;
    config.seed = 7;
    config.rotation_omega_per_myr = 1.913;
    config.rotation_profile_power = 2.0;
    config.rotation_core_cutoff_pc = 0.5;
    TurbulentCloudGenerator generator(config);
    auto realization = generator.generate();
    std::vector<Body> bodies = std::move(realization.particles);

    const double buffer_radius = 0.5 * parsec_km;
    std::vector<Body> outside;
    double buffer_mass = 0.0;
    std::size_t buffer_count = 0;
    for (const Body& body : bodies) {
        if (body.state.pos.norm() < buffer_radius) {
            buffer_mass += body.mass;
            ++buffer_count;
        } else {
            outside.push_back(body);
        }
    }
    bodies = std::move(outside);
    const double buffer_particle_mass =
        buffer_mass / static_cast<double>(buffer_count);
    for (std::size_t k = 0; k < buffer_count; ++k) {
        const double u = (k + 0.5) / buffer_count;
        const double r = buffer_radius * std::cbrt(u);
        const double z = 1.0 - 2.0 * u;
        const double st = std::sqrt(1.0 - z * z);
        const double phi = 2.399963229728653 * k;
        Body body;
        body.name = "buffer";
        body.mass = buffer_particle_mass;
        body.mu = solar::constants::G * buffer_particle_mass;
        body.state.pos = {r * st * std::cos(phi), r * st * std::sin(phi),
                          r * z};
        bodies.push_back(body);
    }

    const double r_c = r_c_pc * parsec_km;
    const std::size_t shell_count = n_clump / 40;
    const std::size_t per_shell = 40;
    const std::size_t grid = 4001;
    std::vector<double> cumulative(grid, 0.0);
    const double du = 3.0 / (grid - 1);
    double total = 0.0;
    for (std::size_t g = 1; g < grid; ++g) {
        const double x = g * du;
        total += x * x * std::exp(-0.5 * x * x) * du;
        cumulative[g] = total;
    }
    const double clump_particle_mass =
        clump_mass_solar * solar_mass_kg /
        static_cast<double>(shell_count * per_shell);
    for (std::size_t shell = 0; shell < shell_count; ++shell) {
        const double quantile = (shell + 0.5) / shell_count;
        std::size_t g = 0;
        while (g + 1 < grid && cumulative[g + 1] / total < quantile) ++g;
        const double f0 = cumulative[g] / total;
        const double f1 = cumulative[g + 1] / total;
        const double x =
            (g + (quantile - f0) / std::max(f1 - f0, 1.0e-30)) * du;
        const double r = r_c * x;
        for (std::size_t p = 0; p < per_shell; ++p) {
            const double u = (p + 0.5) / per_shell;
            const double z = 1.0 - 2.0 * u;
            const double st = std::sqrt(1.0 - z * z);
            const double phi = 2.399963229728653 * (shell + p) + 1.0;
            Body body;
            body.name = "clump";
            body.mass = clump_particle_mass;
            body.mu = solar::constants::G * clump_particle_mass;
            body.state.pos = {r * st * std::cos(phi),
                              r * st * std::sin(phi), r * z};
            bodies.push_back(body);
        }
    }
    return bodies;
}

std::vector<Body> compact_ball(double mass, double radius0,
                               std::size_t count) {
    std::mt19937 generator(7);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::vector<Body> bodies;
    bodies.reserve(count);
    const double particle_mass = mass / count;
    for (std::size_t i = 0; i < count; ++i) {
        const double u = unit(generator);
        const double v = unit(generator);
        const double w = unit(generator);
        const double r = radius0 * std::cbrt(u);
        const double theta = std::acos(2.0 * v - 1.0);
        const double phi = 2.0 * pi * w;
        Body body;
        body.name = "dust";
        body.mass = particle_mass;
        body.mu = solar::constants::G * particle_mass;
        body.state.pos = {r * std::sin(theta) * std::cos(phi),
                          r * std::sin(theta) * std::sin(phi),
                          r * std::cos(theta)};
        bodies.push_back(body);
    }
    return bodies;
}

} // namespace

int main(int argc, char** argv) {
    const std::string out_dir = argc > 1 ? argv[1] : ".";
    const std::string sep = out_dir.empty() || out_dir.back() == '/'
        ? "" : "/";
    const double c = solar::constants::C_LIGHT;
    const int frames = 120;
    Sidecar sidecar;
    char path[512];
    std::mt19937 star_rng(1234);
    std::vector<std::array<int, 3>> stars;
    for (int i = 0; i < 500; ++i) {
        const int x = static_cast<int>(star_rng() % width);
        const int y = static_cast<int>(star_rng() % height);
        const int b = 15 + static_cast<int>(star_rng() % 90);
        stars.push_back({x, y, b});
    }

    // ---------------- Stage 1: nebula-scale collapse -------------------
    {
        const std::size_t n_clump = 800;
        const double r_c_pc = 0.04;
        std::vector<Body> bodies =
            nebula_ic(1.0e4, 8192, n_clump, r_c_pc, 100.0);
        const std::size_t clump_start = bodies.size() - n_clump;
        const double r_c = r_c_pc * parsec_km;
        std::vector<std::size_t> core;
        for (std::size_t i = clump_start; i < bodies.size(); ++i) {
            if (bodies[i].state.pos.norm() < r_c) core.push_back(i);
        }
        double core_mass = 0.0;
        for (const std::size_t i : core) core_mass += bodies[i].mass;
        const double gm_core = solar::constants::G * core_mass;
        const double t_sing =
            pi * std::sqrt(r_c * r_c * r_c / (8.0 * gm_core));
        const double duration = 0.55 * t_sing;

        NBodySim sim;
        sim.init(std::move(bodies));
        sim.clear_forces();
        sim.add_force(std::make_unique<PnCollapseForce>(
            PnCollapseForce::Config{0.5, 1.0e-3 * parsec_km,
                                    1.0e-4, 0.05, 0.1},
            0.0));
        int frame_index = 0;
        sim.run(duration, t_sing / 4000.0, duration / frames,
                [&](double t, const std::vector<Body>& current) {
                    Frame frame;
                    frame.background(stars);
                    // Exponential zoom: 1.3 pc down to 0.12 pc, so the
                    // core collapse is actually visible.
                    const double box_start = 1.3 * parsec_km;
                    const double box_end = 0.12 * parsec_km;
                    const double box = box_start * std::pow(
                        box_end / box_start, t / duration);
                    for (const Body& body : current) {
                        int x, y;
                        project(body.state.pos, box, x, y);
                        const double v = body.state.vel.norm();
                        const double heat = std::min(1.0, v / 1.5);
                        if (body.name == "clump") {
                            frame.splat(x, y, 225 + 30 * heat,
                                        195 - 95 * heat, 85, 0.8);
                        } else if (body.name == "buffer") {
                            frame.splat(x, y, 95, 120, 150, 0.5);
                        } else {
                            frame.splat(x, y, 60 + 60 * heat,
                                        65 + 50 * heat, 80, 0.28);
                        }
                    }
                    double r_core = 0.0;
                    for (const std::size_t i : core) {
                        r_core = std::max(
                            r_core, current[i].state.pos.norm());
                    }
                    const double eps = gm_core / (r_core * c * c);
                    char caption[200];
                    std::snprintf(caption, sizeof(caption),
                        "STAGE 1  nebula collapse | t=%.2f t_sing | "
                        "R_core=%.2f r_c | FOV=%.2f pc | eps=%.2e",
                        t / t_sing, r_core / r_c, box / parsec_km,
                        eps);
                    frame.text(12, 14, caption);
                    std::snprintf(path, sizeof(path),
                                  "%s%snebula-%03d.ppm",
                                  out_dir.c_str(), sep.c_str(),
                                  frame_index);
                    write_ppm(path, frame);
                    sidecar.add("nebula", frame_index, t, r_core, eps,
                                0.0, sim.total_energy(), caption);
                    ++frame_index;
                });
        check("stage 1 rendered >= 120 frames",
              frame_index >= frames);
    }

    // ---------------- Stage 2: compact 1PN collapse --------------------
    {
        const double mass = 10.0 * 1.98892e30;
        const double gm = solar::constants::G * mass;
        const double rg = gm / (c * c);
        const double radius0 = 200.0 * rg;
        const double t_ff = pi * std::sqrt(radius0 * radius0 * radius0 /
                                           (8.0 * gm));
        const double duration = 0.35 * t_ff;
        std::vector<Body> bodies = compact_ball(mass, radius0, 2048);
        NBodySim sim;
        sim.init(std::move(bodies));
        sim.clear_forces();
        sim.add_force(std::make_unique<PostNewtonianGravity>(
            PostNewtonianGravity::Config{0.5, 1.0e-3 * radius0, 0.0},
            0.0));
        int frame_index = 0;
        sim.run(duration, t_ff / 8000.0, duration / frames,
                [&](double t, const std::vector<Body>& current) {
                    Frame frame;
                    frame.background(stars);
                    // Zoom 3 R0 -> 1.1 R0 so the shrinking ball stays
                    // framed; the horizon ring stays visible.
                    const double box = (3.0 - 1.9 * t / duration) *
                                       radius0;
                    for (const Body& body : current) {
                        int x, y;
                        project(body.state.pos, box, x, y);
                        const double v = body.state.vel.norm();
                        const double heat = std::min(1.0, v / 0.15);
                        frame.splat(x, y, 215 + 40 * heat,
                                    175 - 105 * heat, 75, 0.75);
                    }
                    frame.ring(width / 2.0, height / 2.0,
                               (2.0 * rg / box) * width, 225, 65, 65);
                    double r_max = 0.0;
                    for (const Body& body : current) {
                        r_max = std::max(r_max, body.state.pos.norm());
                    }
                    const double eps = gm / (r_max * c * c);
                    char caption[200];
                    std::snprintf(caption, sizeof(caption),
                        "STAGE 2  1PN compact collapse | t=%.2f t_ff | "
                        "R=%.3f R0 | FOV=%.2f R0 | eps=%.3e",
                        t / t_ff, r_max / radius0, box / radius0, eps);
                    frame.text(12, 14, caption);
                    std::snprintf(path, sizeof(path),
                                  "%s%scompact-%03d.ppm",
                                  out_dir.c_str(), sep.c_str(),
                                  frame_index);
                    write_ppm(path, frame);
                    sidecar.add("compact", frame_index, t, r_max, eps,
                                0.0, sim.total_energy(), caption);
                    ++frame_index;
                });
        check("stage 2 rendered >= 120 frames",
              frame_index >= frames);
    }

    // ---------------- Stage 3: OS horizon and frozen star --------------
    {
        const double mass = 10.0 * 1.98892e30;
        const double gm = solar::constants::G * mass;
        const double rg = gm / (c * c);
        const double radius0 = 200.0 * rg;
        const double t_ff = solar::dynamics::os_collapse_time(radius0, mass);
        const double r_obs = 1.0e6 * radius0;
        int frame_index = 0;
        for (; frame_index < frames; ++frame_index) {
            // Sample log-uniformly in delta = R / (2 G M / c^2) - 1 so
            // the freeze-out (the last ~1% of the collapse) is resolved:
            // the redshift and t_obs diverge exactly there.
            const double delta = 99.0 * std::pow(
                1.0e-4 / 99.0,
                static_cast<double>(frame_index) /
                    static_cast<double>(frames - 1));
            const double r_e = 2.0 * rg * (1.0 + delta);
            const double cos_theta = 2.0 * r_e / radius0 - 1.0;
            const double theta_angle = std::acos(cos_theta);
            const double tau = t_ff * (theta_angle +
                                       std::sin(theta_angle)) / pi;
            Frame frame;
            frame.background(stars);
            const double box = 2.5 * radius0;
            // The collapsing surface with a glow whose color tracks the
            // redshift: blue -> red as 1 + z explodes.
            const double r_surface =
                solar::dynamics::os_surface_radius(radius0, mass, tau);
            const double one_plus_z = std::min(
                1.0e4, solar::dynamics::os_surface_redshift(
                           radius0, mass, tau, r_obs));
            const double z_norm = std::log10(one_plus_z) / 4.0;  // 0..1
            frame.glow_ring(width / 2.0, height / 2.0,
                            (r_surface / box) * width,
                            static_cast<int>(90 + 165 * z_norm),
                            static_cast<int>(150 - 110 * z_norm),
                            static_cast<int>(220 - 190 * z_norm));
            // The horizon ring (fixed).
            frame.glow_ring(width / 2.0, height / 2.0,
                            (2.0 * rg / box) * width, 235, 90, 90);
            // Scrolling observer charts: t_obs diverges (blue), the
            // luminosity decays (orange) as the surface freezes.
            const double t_obs =
                solar::dynamics::os_observed_time(radius0, mass, tau,
                                                  r_obs);
            const double lum = solar::dynamics::os_luminosity(
                radius0, mass, tau, r_obs);
            const int chart_x0 = 60;
            const int chart_y = 500;
            frame.fill_rect(56, 386, 336, 518, 20, 28, 48, 0.55);
            const double t_norm = std::min(1.0, t_obs / (1000.0 * t_ff));
            const double l_norm = std::max(
                0.0, std::min(1.0, 1.0 + std::log10(std::max(lum, 1e-12)) /
                                      12.0));
            const int cx = chart_x0 + static_cast<int>(
                240.0 * frame_index / (frames - 1));
            frame.fill_rect(cx, chart_y - static_cast<int>(100 * t_norm),
                            cx + 2, chart_y, 120, 200, 255, 0.95);
            frame.fill_rect(cx + 6,
                            chart_y - static_cast<int>(100 * l_norm),
                            cx + 8, chart_y, 255, 170, 90, 0.95);
            frame.text(64, 394, "t_obs (blue) and L (orange)");
            char caption[200];
            std::snprintf(caption, sizeof(caption),
                "STAGE 3  GR horizon | tau=%.3f t_ff | R=%.3f R0 | "
                "1+z=%.2e | frozen star",
                tau / t_ff, r_surface / radius0, one_plus_z);
            frame.text(12, 14, caption);
            std::snprintf(path, sizeof(path), "%s%shorizon-%03d.ppm",
                          out_dir.c_str(), sep.c_str(), frame_index);
            write_ppm(path, frame);
            sidecar.add("horizon", frame_index, tau, r_surface,
                        gm / (r_surface * c * c), 0.0, 0.0, caption);
        }
        check("stage 3 rendered 120 frames", frame_index == frames);
    }

    sidecar.write(out_dir + sep + "collapse-movie.json");
    std::cout << "  frames and sidecar written to " << out_dir << '\n';
    std::cout << (failures == 0 ? "PASS: collapse movie frames"
                                : "FAILURES PRESENT")
              << '\n';
    return failures == 0 ? 0 : 1;
}