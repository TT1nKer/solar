#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <unistd.h>

namespace {

int passed = 0;
int failed = 0;
int command_sequence = 0;

struct CommandResult {
    int status;
    std::string output;
};

void check(const char* name, bool condition) {
    std::cout << (condition ? "  PASS: " : "  FAIL: ") << name << '\n';
    condition ? ++passed : ++failed;
}

CommandResult run_command(const std::string& command) {
    const std::string output_path =
        "/tmp/solar-relativity-metric-cli-" +
        std::to_string(static_cast<long long>(getpid())) + "-" +
        std::to_string(command_sequence++) + ".out";
    const int status = std::system(
        (command + " > " + output_path + " 2>&1").c_str());

    std::ifstream stream(output_path);
    std::ostringstream output;
    output << stream.rdbuf();
    std::remove(output_path.c_str());
    return {status, output.str()};
}

bool contains(const std::string& text, const std::string& expected) {
    return text.find(expected) != std::string::npos;
}

void check_rejected(const char* name, const std::string& command) {
    check(name, run_command(command).status != 0);
}

} // namespace

int main() {
    const CommandResult kerr = run_command(
        "./solar relativity metric --metric kerr-bl --M 1 --spin 0.9 "
        "--x 0,10,1.5707963267948966,0 --json");
    check("Kerr metric command succeeds", kerr.status == 0);
    check("Kerr JSON identifies the metric",
          contains(kerr.output, "\"metric\":\"kerr-bl\""));
    check("Kerr JSON reports inverse error",
          contains(kerr.output, "\"inverse_error\":"));

    const CommandResult schwarzschild = run_command(
        "./solar relativity metric --metric schwarzschild --M 1 "
        "--x 0,10,1.5707963267948966,0 --json");
    check("Schwarzschild metric command succeeds",
          schwarzschild.status == 0);
    check("Schwarzschild JSON identifies the metric",
          contains(
              schwarzschild.output,
              "\"metric\":\"schwarzschild\""));

    const CommandResult minkowski = run_command(
        "./solar relativity metric --metric minkowski "
        "--x 0,1,2,3 --json");
    check("Minkowski metric command succeeds", minkowski.status == 0);
    check("Minkowski JSON identifies the metric",
          contains(minkowski.output, "\"metric\":\"minkowski\""));

    check_rejected(
        "missing coordinate is rejected",
        "./solar relativity metric --metric kerr-bl "
        "--M 1 --spin 0.9 --json");
    check_rejected(
        "missing mass is rejected",
        "./solar relativity metric --metric schwarzschild "
        "--x 0,10,1.2,0 --json");
    check_rejected(
        "missing spin is rejected",
        "./solar relativity metric --metric kerr-bl "
        "--M 1 --x 0,10,1.2,0 --json");
    check_rejected(
        "extremal spin is rejected",
        "./solar relativity metric --metric kerr-bl "
        "--M 1 --spin 1 --x 0,10,1.2,0 --json");
    check_rejected(
        "unknown metric is rejected",
        "./solar relativity metric --metric nope "
        "--x 0,10,1.2,0 --json");
    check_rejected(
        "malformed coordinate is rejected",
        "./solar relativity metric --metric minkowski "
        "--x 0,1,2 --json");
    check_rejected(
        "non-finite coordinate is rejected",
        "./solar relativity metric --metric minkowski "
        "--x 0,1,nan,3 --json");
    check_rejected(
        "numeric trailing characters are rejected",
        "./solar relativity metric --metric schwarzschild "
        "--M 1oops --x 0,10,1.2,0 --json");
    check_rejected(
        "duplicate option is rejected",
        "./solar relativity metric --metric minkowski "
        "--x 0,1,2,3 --x 0,4,5,6 --json");
    check_rejected(
        "Kerr horizon point is rejected",
        "./solar relativity metric --metric kerr-bl "
        "--M 1 --spin 0.9 --x 0,1.4358898943540672,1.2,0 --json");
    check_rejected(
        "unknown option is rejected",
        "./solar relativity metric --metric minkowski "
        "--x 0,1,2,3 --wat");

    std::cout << "\n=== Results: " << passed << " passed, " << failed
              << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
