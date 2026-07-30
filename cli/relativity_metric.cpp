#include "relativity_metric.h"

#include "solar/relativity/kerr_bl_metric.h"
#include "solar/relativity/minkowski_metric.h"
#include "solar/relativity/schwarzschild_metric.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

namespace solar::cli {
namespace {

using relativity::Contravariant4;
using relativity::Mat4;
using relativity::Metric;
using relativity::Vec4;

struct MetricCommandOptions {
    std::string metric_name;
    double mass_M = 0.0;
    double spin_chi = 0.0;
    Contravariant4 point{};
    bool has_metric = false;
    bool has_mass = false;
    bool has_spin = false;
    bool has_point = false;
    bool json = false;
};

double parse_finite_double(
    const std::string& text, const char* option_name) {
    std::size_t consumed = 0;
    double value = 0.0;
    try {
        value = std::stod(text, &consumed);
    } catch (const std::exception&) {
        throw std::invalid_argument(
            std::string(option_name) + " requires a number");
    }
    if (consumed != text.size() || !std::isfinite(value)) {
        throw std::invalid_argument(
            std::string(option_name) + " requires a finite number");
    }
    return value;
}

Contravariant4 parse_point(const std::string& text) {
    Vec4 coordinates{};
    std::size_t token_start = 0;
    for (std::size_t coordinate = 0; coordinate < 4; ++coordinate) {
        const std::size_t comma = text.find(',', token_start);
        if ((coordinate < 3 && comma == std::string::npos) ||
            (coordinate == 3 && comma != std::string::npos)) {
            throw std::invalid_argument(
                "--x requires exactly four comma-separated numbers");
        }
        const std::size_t token_end =
            comma == std::string::npos ? text.size() : comma;
        if (token_end == token_start) {
            throw std::invalid_argument(
                "--x coordinate values cannot be empty");
        }
        coordinates[coordinate] = parse_finite_double(
            text.substr(token_start, token_end - token_start), "--x");
        token_start = token_end + 1;
    }
    return Contravariant4{coordinates};
}

const char* require_option_value(
    int argc, char* argv[], int& index, const char* option_name) {
    if (index + 1 >= argc) {
        throw std::invalid_argument(
            std::string(option_name) + " requires a value");
    }
    return argv[++index];
}

void reject_duplicate(bool already_present, const char* option_name) {
    if (already_present) {
        throw std::invalid_argument(
            std::string(option_name) + " may be specified only once");
    }
}

MetricCommandOptions parse_metric_options(int argc, char* argv[]) {
    MetricCommandOptions options;
    for (int index = 3; index < argc; ++index) {
        const std::string option = argv[index];
        if (option == "--metric") {
            reject_duplicate(options.has_metric, "--metric");
            options.metric_name = require_option_value(
                argc, argv, index, "--metric");
            options.has_metric = true;
        } else if (option == "--M") {
            reject_duplicate(options.has_mass, "--M");
            options.mass_M = parse_finite_double(
                require_option_value(argc, argv, index, "--M"), "--M");
            options.has_mass = true;
        } else if (option == "--spin") {
            reject_duplicate(options.has_spin, "--spin");
            options.spin_chi = parse_finite_double(
                require_option_value(argc, argv, index, "--spin"),
                "--spin");
            options.has_spin = true;
        } else if (option == "--x") {
            reject_duplicate(options.has_point, "--x");
            options.point = parse_point(
                require_option_value(argc, argv, index, "--x"));
            options.has_point = true;
        } else if (option == "--json") {
            reject_duplicate(options.json, "--json");
            options.json = true;
        } else {
            throw std::invalid_argument("unknown metric option: " + option);
        }
    }

    if (!options.has_metric) {
        throw std::invalid_argument("--metric is required");
    }
    if (!options.has_point) {
        throw std::invalid_argument("--x is required");
    }
    return options;
}

std::unique_ptr<Metric> create_metric(
    const MetricCommandOptions& options) {
    if (options.metric_name == "minkowski") {
        if (options.has_mass || options.has_spin) {
            throw std::invalid_argument(
                "Minkowski metric does not accept --M or --spin");
        }
        return std::make_unique<relativity::MinkowskiMetric>();
    }
    if (options.metric_name == "schwarzschild") {
        if (!options.has_mass) {
            throw std::invalid_argument(
                "Schwarzschild metric requires --M");
        }
        if (options.has_spin) {
            throw std::invalid_argument(
                "Schwarzschild metric does not accept --spin");
        }
        return std::make_unique<
            relativity::SchwarzschildBoyerLindquistMetric>(
                options.mass_M);
    }
    if (options.metric_name == "kerr-bl") {
        if (!options.has_mass || !options.has_spin) {
            throw std::invalid_argument(
                "Kerr BL metric requires --M and --spin");
        }
        return std::make_unique<
            relativity::KerrBoyerLindquistMetric>(
                options.mass_M, options.spin_chi);
    }
    throw std::invalid_argument(
        "unknown metric: " + options.metric_name);
}

double inverse_identity_error(
    const Mat4& covariant, const Mat4& contravariant) {
    const Mat4 product = relativity::multiply(covariant, contravariant);
    double maximum = 0.0;
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            const double expected = row == column ? 1.0 : 0.0;
            maximum = std::max(
                maximum, std::fabs(product[row][column] - expected));
        }
    }
    return maximum;
}

void print_json_vector(const Vec4& vector) {
    std::cout << '[';
    for (std::size_t index = 0; index < 4; ++index) {
        if (index != 0) {
            std::cout << ',';
        }
        std::cout << vector[index];
    }
    std::cout << ']';
}

void print_json_matrix(const Mat4& matrix) {
    std::cout << '[';
    for (std::size_t row = 0; row < 4; ++row) {
        if (row != 0) {
            std::cout << ',';
        }
        std::cout << '[';
        for (std::size_t column = 0; column < 4; ++column) {
            if (column != 0) {
                std::cout << ',';
            }
            std::cout << matrix[row][column];
        }
        std::cout << ']';
    }
    std::cout << ']';
}

void print_json(
    const Metric& metric,
    const Contravariant4& point,
    const Mat4& covariant,
    const Mat4& contravariant,
    double inverse_error) {
    std::cout << std::setprecision(
        std::numeric_limits<double>::max_digits10);
    std::cout << "{\"metric\":\"" << metric.name()
              << "\",\"chart\":\""
              << relativity::chart_name(metric.chart())
              << "\",\"x\":";
    print_json_vector(point.v);
    std::cout << ",\"covariant\":";
    print_json_matrix(covariant);
    std::cout << ",\"contravariant\":";
    print_json_matrix(contravariant);
    std::cout << ",\"inverse_error\":" << inverse_error << "}\n";
}

void print_human_matrix(const char* label, const Mat4& matrix) {
    std::cout << label << ":\n";
    for (const auto& row : matrix) {
        std::cout << "  ";
        for (double value : row) {
            std::cout << std::setw(24) << value;
        }
        std::cout << '\n';
    }
}

void print_human(
    const Metric& metric,
    const Contravariant4& point,
    const Mat4& covariant,
    const Mat4& contravariant,
    double inverse_error) {
    std::cout << std::setprecision(
        std::numeric_limits<double>::max_digits10);
    std::cout << "Metric: " << metric.name() << '\n'
              << "Chart: " << relativity::chart_name(metric.chart())
              << "\nCoordinates: ";
    for (std::size_t index = 0; index < 4; ++index) {
        if (index != 0) {
            std::cout << ", ";
        }
        std::cout << point.v[index];
    }
    std::cout << '\n';
    print_human_matrix("Covariant g_mu_nu", covariant);
    print_human_matrix("Contravariant g^mu_nu", contravariant);
    std::cout << "Inverse identity max error: "
              << inverse_error << '\n';
}

int inspect_metric(int argc, char* argv[]) {
    const MetricCommandOptions options =
        parse_metric_options(argc, argv);
    const std::unique_ptr<Metric> metric = create_metric(options);
    if (!metric->valid_point(options.point)) {
        throw std::domain_error(
            "coordinate is outside the metric's valid domain");
    }

    const Mat4 covariant = metric->covariant(options.point);
    const Mat4 contravariant = metric->contravariant(options.point);
    const double inverse_error =
        inverse_identity_error(covariant, contravariant);
    if (!std::isfinite(inverse_error)) {
        throw std::runtime_error(
            "metric inverse validation produced a non-finite error");
    }

    if (options.json) {
        print_json(
            *metric, options.point, covariant, contravariant,
            inverse_error);
    } else {
        print_human(
            *metric, options.point, covariant, contravariant,
            inverse_error);
    }
    return 0;
}

} // namespace

int dispatch_relativity(int argc, char* argv[]) {
    if (argc < 3) {
        throw std::invalid_argument(
            "usage: solar relativity metric [options]");
    }
    const std::string subcommand = argv[2];
    if (subcommand == "metric") {
        return inspect_metric(argc, argv);
    }
    throw std::invalid_argument(
        "unknown relativity subcommand: " + subcommand);
}

} // namespace solar::cli
