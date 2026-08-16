#include "solar/constants.h"
#include "solar/dynamics/ltb_collapse.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

// The outside-observer view of gravitational collapse, rendered with the
// exact Schwarzschild null-geodesic physics:
//
//  - a pixel at apparent angle theta corresponds to impact parameter
//    b = D tan(theta) ~ D theta; rays with b below the star's image
//    radius b_star(R) = R / sqrt(1 - 2 G M / (R c^2)) hit the surface
//    (for R above the photon sphere 3 G M / c^2) and rays with
//    b < 3 sqrt(3) G M / c^2 are captured by the hole;
//  - the image radius therefore shrinks from ~R0 to the shadow size
//    3 sqrt(3) G M / c^2 and FREEZES there once the surface passes the
//    photon sphere - the frozen star IS the shadow;
//  - the surface's observed temperature is T_em / (1+z) and its
//    surface brightness (1+z)^-4, so it reddens and fades while the
//    observer time diverges logarithmically;
//  - rays with b > b_crit escape; each pixel then samples the
//    background starfield at the lensed sky angle beta = theta +
//    alpha(b), with alpha(b) the exact Schwarzschild bending angle
//    (precomputed quadrature table) - the Einstein-ring lensing that is
//    all that remains once the star has faded.
//
// Usage: test_obsview_movie <output-directory>
// ffmpeg -y -framerate 30 -i DIR/obs-%03d.ppm -c:v libx264 -pix_fmt
// yuv420p obsview.mp4

namespace {

int failures = 0;
void check(const std::string& name, bool condition) {
    std::cout << (condition ? "  PASS: " : "  FAIL: ") << name << '\n';
    if (!condition) ++failures;
}
constexpr double pi = 3.14159265358979323846;
const int width = 480;
const int height = 270;

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



// Blackbody temperature (K) -> RGB (Tanner Helland approximation).
void temperature_rgb(double t_kelvin, double& r, double& g, double& b) {
    const double t = std::max(1000.0, t_kelvin) / 100.0;
    r = t <= 66 ? 255.0
                : 329.698727446 * std::pow(t - 60.0, -0.1332047592);
    g = t <= 66 ? 99.4708025861 * std::log(t) - 161.1195681661
                : 288.1221695283 * std::pow(t - 60.0, -0.0755148492);
    b = t >= 66 ? 255.0
        : (t <= 19 ? 0.0
                   : 138.5177312231 * std::log(t - 10.0) - 305.0447927307);
    r = std::min(255.0, std::max(0.0, r));
    g = std::min(255.0, std::max(0.0, g));
    b = std::min(255.0, std::max(0.0, b));
}

struct Star {
    double angle_x;
    double angle_y;
    int brightness;
};

// Exact Schwarzschild bending angle alpha(b): the total deflection of a
// photon with impact parameter b, integrated over u = 1/r.
double bending_angle(double b, double gm_over_c2) {
    const double b_crit = 3.0 * std::sqrt(3.0) * gm_over_c2;
    if (b <= b_crit) return 1.0e30;
    double u0 = 1.0 / b;
    for (int i = 0; i < 20; ++i) {
        const double f = u0 * u0 * (1.0 - 2.0 * gm_over_c2 * u0) -
                         1.0 / (b * b);
        const double fp = 2.0 * u0 - 6.0 * gm_over_c2 * u0 * u0;
        if (fp == 0.0) break;
        const double step = f / fp;
        u0 -= step;
        if (std::fabs(step) < 1e-15 * u0) break;
    }
    const int n = 512;
    double sum = 0.0;
    for (int i = 0; i <= n; ++i) {
        const double u = u0 * static_cast<double>(i) / n;
        if (i == n) continue;  // root endpoint contributes zero
        const double denom = 1.0 / (b * b) - u * u +
                             2.0 * gm_over_c2 * u * u * u;
        double f = 1.0 / std::sqrt(std::max(denom, 1.0e-300));
        if (i == 0) f = b;  // u = 0 endpoint: integrand = b exactly
        if (i % 2 == 1) f *= 4.0;
        else if (i != 0) f *= 2.0;
        sum += f;
    }
    return 2.0 * sum * u0 / (3.0 * n) - pi;
}

} // namespace

int main(int argc, char** argv) {
    const std::string out_dir = argc > 1 ? argv[1] : ".";
    const std::string sep = out_dir.empty() || out_dir.back() == '/'
        ? "" : "/";
    const double c = solar::constants::C_LIGHT;
    const double mass = 10.0 * 1.98892e30;
    const double gm = solar::constants::G * mass;
    const double gm_c2 = gm / (c * c);
    const double b_crit = 3.0 * std::sqrt(3.0) * gm_c2;
    const double radius0 = 200.0 * gm_c2;
    const double t_ff = solar::dynamics::os_collapse_time(radius0, mass);
    const double t_em = 6000.0;
    const double observer_distance = 2.0e4 * b_crit;

    // Deflection table on a log grid of b.
    const int table_n = 2048;
    std::vector<double> table_b(table_n), table_alpha(table_n);
    const double b_lo = b_crit * 1.0005;
    const double b_hi = b_crit * 40.0;
    for (int i = 0; i < table_n; ++i) {
        table_b[i] = b_lo * std::pow(b_hi / b_lo,
                                     static_cast<double>(i) / (table_n - 1));
        table_alpha[i] = bending_angle(table_b[i], gm_c2);
    }
    auto alpha_at = [&](double b) {
        if (b <= b_crit) return 1.0e30;
        const double x = std::log(b / b_lo) / std::log(b_hi / b_lo);
        const double idx = std::min(
            static_cast<double>(table_n - 2),
            std::max(0.0, x * (table_n - 1)));
        const int i = static_cast<int>(idx);
        const double frac = idx - i;
        return table_alpha[i] * (1.0 - frac) + table_alpha[i + 1] * frac;
    };

    // Background starfield on the sky (lensed by the hole); the sky must
    // span the deflected source angles: beta = theta + alpha(b), with
    // alpha ~ 0.25 rad at the frame edge.
    std::mt19937 rng(77);
    std::vector<Star> stars;
    const double sky_angle = 0.32;
    for (int i = 0; i < 1600; ++i) {
        const double r = sky_angle * std::sqrt(
            static_cast<double>(rng() % 100000) / 100000.0);
        const double phi = 2.0 * pi *
            static_cast<double>(rng() % 100000) / 100000.0;
        stars.push_back({r * std::cos(phi), r * std::sin(phi),
                         40 + static_cast<int>(rng() % 180)});
    }
    const int grid_n = 512;
    const double cell = 2.0 * sky_angle / grid_n;
    std::vector<std::vector<int>> star_grid(grid_n * grid_n);
    for (int i = 0; i < static_cast<int>(stars.size()); ++i) {
        const int gx = static_cast<int>((stars[i].angle_x + sky_angle) /
                                        cell);
        const int gy = static_cast<int>((stars[i].angle_y + sky_angle) /
                                        cell);
        if (gx >= 0 && gx < grid_n && gy >= 0 && gy < grid_n) {
            star_grid[gy * grid_n + gx].push_back(i);
        }
    }
    auto star_at = [&](double bx, double by) {
        if (std::abs(bx) > sky_angle || std::abs(by) > sky_angle) {
            return 0.0;
        }
        const int gx = static_cast<int>((bx + sky_angle) / cell);
        const int gy = static_cast<int>((by + sky_angle) / cell);
        for (const int index : star_grid[gy * grid_n + gx]) {
            const Star& star = stars[index];
            const double dx = star.angle_x - bx;
            const double dy = star.angle_y - by;
            if (dx * dx + dy * dy < 3.0e-9) {
                return static_cast<double>(star.brightness);
            }
        }
        return 0.0;
    };

    const int collapse_frames = 120;
    const int settled_frames = 40;
    const int total = collapse_frames + settled_frames;
    char path[512];
    int bright_pixels_frame0 = 0;
    int bright_pixels_last = 0;
    bool image_froze_at_shadow = true;

    for (int frame_index = 0; frame_index < total; ++frame_index) {
        double r_surface;
        double zp1 = 1.0;
        double t_obs = 0.0;
        if (frame_index < collapse_frames) {
            const double log_r = std::log(radius0) +
                (std::log(2.0 * gm_c2 * 1.000002) - std::log(radius0)) *
                    static_cast<double>(frame_index) /
                    (collapse_frames - 1);
            r_surface = std::exp(log_r);
            const double cos_theta = 2.0 * r_surface / radius0 - 1.0;
            const double theta_angle = std::acos(cos_theta);
            const double tau = t_ff * (theta_angle +
                                       std::sin(theta_angle)) / pi;
            zp1 = solar::dynamics::os_surface_redshift(
                radius0, mass, tau, observer_distance);
            t_obs = solar::dynamics::os_observed_time(
                radius0, mass, tau, observer_distance);
        } else {
            r_surface = 2.0 * gm_c2;
        }
        double b_star;
        if (r_surface >= 3.0 * gm_c2) {
            b_star = r_surface / std::sqrt(1.0 - 2.0 * gm_c2 / r_surface);
        } else {
            b_star = b_crit;
        }
        if (frame_index == collapse_frames - 1 &&
            std::abs(b_star - b_crit) > 1.0e-9 * b_crit) {
            image_froze_at_shadow = false;
        }

        Frame frame;
        frame.background({});
        // Collapse phase: tight camera on the shadow. Settled phase:
        // wide camera so the Einstein ring (~124 shadow radii out) is
        // visible around the tiny shadow.
        const bool wide = frame_index >= collapse_frames;
        const double wide_distance = 1.0e5 * b_crit;
        const double max_b = wide ? 300.0 * b_crit : 3.2 * b_crit;
        const double cam_distance =
            wide ? wide_distance : observer_distance;
        const double theta_max = max_b / cam_distance;
        const double zp4 = 1.0 / (zp1 * zp1 * zp1 * zp1);
        double tr, tg, tb;
        temperature_rgb(t_em / zp1, tr, tg, tb);
        const double surface_alpha = std::min(1.0, 2.0 * zp4);
        int bright = 0;
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                const double thx = (x - width / 2.0) / (width / 2.0) *
                                   theta_max;
                const double thy = (y - height / 2.0) / (height / 2.0) *
                                   theta_max * width / height;
                const double theta = std::sqrt(thx * thx + thy * thy);
                const double b = cam_distance * theta;
                if (b <= b_star && !wide) {
                    frame.blend(x, y, tr * surface_alpha,
                                tg * surface_alpha, tb * surface_alpha,
                                1.0);
                    if (tr * surface_alpha > 5.0) ++bright;  // visible
                } else if (b <= b_crit) {
                    // captured: the shadow
                } else {
                    const double alpha = alpha_at(b);
                    const double scale = 1.0 - alpha / theta;
                    const double light = star_at(thx * scale, thy * scale);
                    if (light > 0.0) {
                        frame.blend(x, y, light, light, light * 1.05, 1.0);
                    }
                }
            }
        }
        if (frame_index == 0) bright_pixels_frame0 = bright;
        if (frame_index == collapse_frames - 1) bright_pixels_last = bright;
        char caption[200];
        if (frame_index < collapse_frames) {
            std::snprintf(caption, sizeof(caption),
                "OUTSIDE VIEW | T_OBS=%.3g S | R=%.4g KM | 1+Z=%.3g | "
                "IMAGE=%.2f X SHADOW",
                t_obs, r_surface, zp1, b_star / b_crit);
        } else {
            std::snprintf(caption, sizeof(caption),
                "SETTLED | SCHWARZSCHILD SHADOW | LENSED SKY | "
                "AFTERGLOW E-FOLD ~ GM/C3 = 49 US");
        }
        frame.text(8, 8, caption);
        std::snprintf(path, sizeof(path), "%s%sobs-%03d.ppm",
                      out_dir.c_str(), sep.c_str(), frame_index);
        write_ppm(path, frame);
    }

    check("image froze at the shadow size when R < 3 G M / c^2",
          image_froze_at_shadow);
    std::cout << "  bright pixels: frame 0 = " << bright_pixels_frame0
              << ", last collapse frame = " << bright_pixels_last << '\n';
    check("surface fades to black by the last collapse frame",
          bright_pixels_frame0 > 100000 && bright_pixels_last == 0);

    std::cout << (failures == 0 ? "PASS: outside-observer view frames"
                                : "FAILURES PRESENT")
              << '\n';
    return failures == 0 ? 0 : 1;
}