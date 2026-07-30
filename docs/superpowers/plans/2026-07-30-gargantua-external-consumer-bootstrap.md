# Gargantua External Consumer Bootstrap Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create an independent Gargantua Studio repository that builds one real Kerr shadow probe through Solar's public Phase 2 API, adding only the minimum CMake package boundary Solar needs for that consumer.

**Architecture:** Solar remains the CPU physics authority and exports `Solar::Core` and `Solar::Relativity` CMake targets without changing the existing Makefile workflow. Gargantua is a sibling Git repository with a locked Solar commit, a local-source override for development, and a probe that prints versioned Phase 2 Bardeen shadow evidence. No CUDA, scene system, radiation model, or film output enters this bootstrap.

**Tech Stack:** C++17, CMake 3.20+, existing Solar Makefile, CTest, Git.

## Global Constraints

- Solar Phase 2 is the first usable baseline; this plan does not claim Solar is feature-complete or API-stable.
- Dependency direction is `Gargantua -> Solar`; Solar must not reference Gargantua.
- Gargantua may include only installed/public `include/solar/...` headers.
- Preserve `make`, `make test`, all existing public APIs, and the zero-mandatory-dependency Solar CPU core.
- The initial contract identifiers are version `0.2.0-alpha.1` and physics contract `relativity-v3-phase2`.
- Gargantua must reject a local Solar checkout whose Git commit differs from its lock.
- No CUDA, OpenEXR, OCIO, FFmpeg, renderer, scene parser, or empty future module is added in this plan.
- Use TDD, focused commits, and local verification before any push.

---

## File map

### Solar repository

- `CMakeLists.txt`: build, alias, install, and export the two public library targets.
- `cmake/SolarConfig.cmake.in`: installed `find_package(Solar CONFIG)` entry point.
- `include/solar/version.h`: public version and physics-contract constants.
- `tests/external_consumer/CMakeLists.txt`: a package consumer with no access to Solar source internals.
- `tests/external_consumer/probe.cpp`: verifies the installed package and public Kerr shadow API.
- `tests/test_external_consumer.sh`: isolated configure/build/install/consume acceptance flow.
- `Makefile`: adds a discoverable `test-external-consumer` target without changing `make test`.
- `.github/workflows/ci.yml`: runs the external-consumer gate once per candidate.
- `docs/relativity/EXTERNAL_CONSUMERS.md`: records the development-package contract and limitations.

### Gargantua repository

- `CMakeLists.txt`: defines the first executable and CTest contract.
- `cmake/SolarDependency.cmake`: loads the lock, verifies a local checkout, or fetches the locked commit.
- `cmake/WriteSolarLock.cmake`: records an exact local Solar commit without a mutable branch reference.
- `cmake/solar-lock.cmake`: immutable Solar repository, commit, version, and physics-contract values.
- `src/probe_main.cpp`: calls Solar's public Bardeen curve API and prints compact JSON.
- `.gitignore`: excludes local build products.
- `README.md`: build and ownership boundary.
- `docs/validation/00_solar_consumer.md`: reproducible evidence for the first cross-repository call.

---

### Task 1: Solar public version contract

**Files:**
- Create: `include/solar/version.h`
- Create: `tests/external_consumer/probe.cpp`
- Create: `tests/external_consumer/CMakeLists.txt`

**Interfaces:**
- Consumes: `solar::relativity::KerrBoyerLindquistMetric` and `bardeen_shadow_curve`.
- Produces: `solar::version`, `solar::physics_contract`, and an external executable named `solar_external_consumer_probe`.

- [ ] **Step 1: Add the failing external consumer source**

Create `tests/external_consumer/probe.cpp`:

```cpp
#include "solar/relativity/kerr_shadow.h"
#include "solar/version.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string_view>
#include <vector>

int main() {
    using solar::relativity::KerrBoyerLindquistMetric;
    using solar::relativity::ShadowCriticalPoint;
    using solar::relativity::bardeen_shadow_curve;

    constexpr double half_pi = 1.5707963267948966;
    const KerrBoyerLindquistMetric metric(1.0, 0.5);
    const std::vector<ShadowCriticalPoint> curve =
        bardeen_shadow_curve(metric, half_pi, 65);
    const auto edges = std::minmax_element(
        curve.begin(),
        curve.end(),
        [](const ShadowCriticalPoint& left,
           const ShadowCriticalPoint& right) {
            return left.alpha < right.alpha;
        });

    const bool valid =
        solar::version == std::string_view{"0.2.0-alpha.1"} &&
        solar::physics_contract ==
            std::string_view{"relativity-v3-phase2"} &&
        curve.size() == 128 &&
        edges.first != curve.end() &&
        std::fabs(
            edges.first->alpha -
            (-4.096266658713869)) < 1.0e-13 &&
        std::fabs(
            edges.second->alpha -
            6.138155724715452) < 1.0e-13;
    if (!valid) {
        std::cerr << "installed Solar package contract mismatch\n";
        return 1;
    }

    std::cout
        << "{\"solar_version\":\"" << solar::version
        << "\",\"physics_contract\":\""
        << solar::physics_contract
        << "\",\"samples\":" << curve.size()
        << ",\"left\":" << edges.first->alpha
        << ",\"right\":" << edges.second->alpha
        << "}\n";
    return 0;
}
```

- [ ] **Step 2: Add the consumer-only CMake project**

Create `tests/external_consumer/CMakeLists.txt`:

```cmake
cmake_minimum_required(VERSION 3.20)
project(SolarExternalConsumer LANGUAGES CXX)

find_package(Solar CONFIG REQUIRED)

add_executable(
  solar_external_consumer_probe
  probe.cpp
)
target_compile_features(
  solar_external_consumer_probe
  PRIVATE cxx_std_17
)
target_link_libraries(
  solar_external_consumer_probe
  PRIVATE Solar::Relativity
)
```

- [ ] **Step 3: Run the consumer configure to verify RED**

Run:

```bash
cmake -S tests/external_consumer \
  -B /tmp/solar-external-consumer-red
```

Expected: configure fails because no installed `SolarConfig.cmake` exists.

- [ ] **Step 4: Add the public version header**

Create `include/solar/version.h`:

```cpp
#pragma once

#include <string_view>

namespace solar {

inline constexpr std::string_view version{
    "0.2.0-alpha.1"};
inline constexpr std::string_view physics_contract{
    "relativity-v3-phase2"};

} // namespace solar
```

- [ ] **Step 5: Compile the existing Makefile suite**

Run:

```bash
make clean
make
make tests/relativity/test_kerr_shadow
./tests/relativity/test_kerr_shadow
```

Expected: existing build succeeds and the Kerr shadow executable reports zero failures.

- [ ] **Step 6: Commit the version contract and RED consumer**

```bash
git add \
  include/solar/version.h \
  tests/external_consumer/CMakeLists.txt \
  tests/external_consumer/probe.cpp
git commit -m "test: define Solar external consumer contract"
```

---

### Task 2: Installable Solar CMake package

**Files:**
- Create: `CMakeLists.txt`
- Create: `cmake/SolarConfig.cmake.in`

**Interfaces:**
- Consumes: all `src/*.cpp`, `src/relativity/*.cpp`, and `include/solar/...`.
- Produces: build-tree aliases and installed targets `Solar::Core` and `Solar::Relativity`; options `SOLAR_BUILD_CLI`.

- [ ] **Step 1: Add the package configuration template**

Create `cmake/SolarConfig.cmake.in`:

```cmake
@PACKAGE_INIT@

include("${CMAKE_CURRENT_LIST_DIR}/SolarTargets.cmake")
check_required_components(Solar)
```

- [ ] **Step 2: Add the top-level CMake build**

Create `CMakeLists.txt`:

```cmake
cmake_minimum_required(VERSION 3.20)
project(Solar VERSION 0.2.0 LANGUAGES CXX)

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

option(SOLAR_BUILD_CLI "Build the Solar command-line application" ON)

file(
  GLOB SOLAR_CORE_SOURCES
  CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)
file(
  GLOB SOLAR_RELATIVITY_SOURCES
  CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/src/relativity/*.cpp"
)

add_library(solar_core STATIC ${SOLAR_CORE_SOURCES})
add_library(Solar::Core ALIAS solar_core)
set_target_properties(
  solar_core
  PROPERTIES
    EXPORT_NAME Core
    OUTPUT_NAME solar_core
)
target_compile_features(solar_core PUBLIC cxx_std_17)
target_include_directories(
  solar_core
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

add_library(
  solar_relativity
  STATIC
  ${SOLAR_RELATIVITY_SOURCES}
)
add_library(Solar::Relativity ALIAS solar_relativity)
set_target_properties(
  solar_relativity
  PROPERTIES
    EXPORT_NAME Relativity
    OUTPUT_NAME solar_relativity
)
target_compile_features(solar_relativity PUBLIC cxx_std_17)
target_include_directories(
  solar_relativity
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_link_libraries(solar_relativity PUBLIC Solar::Core)

if(SOLAR_BUILD_CLI)
  file(
    GLOB SOLAR_CLI_SOURCES
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/cli/*.cpp"
  )
  add_executable(solar ${SOLAR_CLI_SOURCES})
  target_compile_features(solar PRIVATE cxx_std_17)
  target_link_libraries(
    solar
    PRIVATE Solar::Core Solar::Relativity
  )
endif()

install(
  TARGETS solar_core solar_relativity
  EXPORT SolarTargets
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
)
install(
  DIRECTORY include/solar
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)
install(
  EXPORT SolarTargets
  FILE SolarTargets.cmake
  NAMESPACE Solar::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Solar"
)

configure_package_config_file(
  cmake/SolarConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/SolarConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Solar"
)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/SolarConfigVersion.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY SameMajorVersion
)
install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/SolarConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/SolarConfigVersion.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Solar"
)
```

- [ ] **Step 3: Configure and build the two Solar libraries**

Run:

```bash
cmake -S . \
  -B /tmp/solar-package-build \
  -DSOLAR_BUILD_CLI=OFF \
  -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/solar-package-build --parallel
```

Expected: `solar_core` and `solar_relativity` build without warnings introduced by the new configuration.

- [ ] **Step 4: Install to an isolated prefix**

Run:

```bash
cmake --install /tmp/solar-package-build \
  --prefix /tmp/solar-package-prefix
```

Expected: the prefix contains public headers, both static libraries, `SolarConfig.cmake`, `SolarConfigVersion.cmake`, and `SolarTargets.cmake`.

- [ ] **Step 5: Configure and run the external consumer**

Run:

```bash
cmake -S tests/external_consumer \
  -B /tmp/solar-external-consumer-green \
  -DCMAKE_PREFIX_PATH=/tmp/solar-package-prefix
cmake --build /tmp/solar-external-consumer-green
/tmp/solar-external-consumer-green/solar_external_consumer_probe
```

Expected output contains:

```json
{"solar_version":"0.2.0-alpha.1","physics_contract":"relativity-v3-phase2","samples":128
```

and the process exits zero.

- [ ] **Step 6: Run the existing Makefile build**

Run:

```bash
make clean
make
make test
```

Expected: all existing relativity and fixture-independent legacy assertions pass; the optional DE440 fixture may remain a visible skip.

- [ ] **Step 7: Commit the CMake package**

```bash
git add CMakeLists.txt cmake/SolarConfig.cmake.in
git commit -m "build: export Solar CMake package"
```

---

### Task 3: Reproducible Solar package acceptance gate

**Files:**
- Create: `tests/test_external_consumer.sh`
- Modify: `Makefile`
- Modify: `.github/workflows/ci.yml`
- Create: `docs/relativity/EXTERNAL_CONSUMERS.md`

**Interfaces:**
- Consumes: the installed package from Task 2.
- Produces: `make test-external-consumer` and a documented `0.2.0-alpha.1` development contract.

- [ ] **Step 1: Add the isolated package test script**

Create `tests/test_external_consumer.sh`:

```bash
#!/usr/bin/env bash
set -euo pipefail

solar_root="$(
  cd "$(dirname "${BASH_SOURCE[0]}")/.." &&
  pwd
)"
task_root="$(mktemp -d "${TMPDIR:-/tmp}/solar-consumer.XXXXXX")"
trap 'rm -rf "$task_root"' EXIT

cmake \
  -S "$solar_root" \
  -B "$task_root/solar-build" \
  -DSOLAR_BUILD_CLI=OFF \
  -DCMAKE_BUILD_TYPE=Release
cmake --build "$task_root/solar-build" --parallel
cmake \
  --install "$task_root/solar-build" \
  --prefix "$task_root/prefix"

cmake \
  -S "$solar_root/tests/external_consumer" \
  -B "$task_root/consumer-build" \
  -DCMAKE_PREFIX_PATH="$task_root/prefix"
cmake --build "$task_root/consumer-build" --parallel
"$task_root/consumer-build/solar_external_consumer_probe"
```

- [ ] **Step 2: Make the script executable and run it**

Run:

```bash
chmod +x tests/test_external_consumer.sh
./tests/test_external_consumer.sh
```

Expected: configure, install, consumer build, and probe all succeed in a new temporary directory.

- [ ] **Step 3: Add a Makefile entry without changing `make test`**

Modify `.PHONY` and append:

```make
.PHONY: all clean test test-external-consumer

test-external-consumer:
	./tests/test_external_consumer.sh
```

Keep the existing `test` recipe unchanged.

- [ ] **Step 4: Add the package gate to CI**

Append after the existing `Run test suite` step in `.github/workflows/ci.yml`:

```yaml
      - name: Verify external CMake consumer
        run: make test-external-consumer
```

- [ ] **Step 5: Document the development contract**

Create `docs/relativity/EXTERNAL_CONSUMERS.md` with:

```markdown
# Solar external consumers

Solar relativity Phase 2 is available as a development library contract. It
is not a stable 1.0 API and does not imply completion of the v3 roadmap.

## Contract

- Package version: `0.2.0-alpha.1`
- Physics contract: `relativity-v3-phase2`
- CMake targets: `Solar::Core`, `Solar::Relativity`
- Language level: C++17
- Mandatory runtime dependencies: none

An external project consumes an installed package with:

```cmake
find_package(Solar CONFIG REQUIRED)
target_link_libraries(app PRIVATE Solar::Relativity)
```

Run the isolated acceptance check with:

```bash
make test-external-consumer
```

Consumers must pin an immutable Solar commit during development and record
both `solar::version` and `solar::physics_contract` in generated artifacts.
No compatibility is promised between `0.x-alpha` contracts.
```

- [ ] **Step 6: Run focused and full verification**

Run:

```bash
make test-external-consumer
make clean
make
make test
git diff --check
```

Expected: package probe passes; all existing tests pass; diff check is clean.

- [ ] **Step 7: Commit the package gate and documentation**

```bash
git add \
  .github/workflows/ci.yml \
  Makefile \
  docs/relativity/EXTERNAL_CONSUMERS.md \
  tests/test_external_consumer.sh
git commit -m "test: verify Solar external consumers"
```

---

### Task 4: Independent Gargantua repository and locked Solar dependency

**Files:**
- Create repository: `/Users/hostsjim/project/gargantua-studio`
- Create: `/Users/hostsjim/project/gargantua-studio/CMakeLists.txt`
- Create: `/Users/hostsjim/project/gargantua-studio/cmake/SolarDependency.cmake`
- Create: `/Users/hostsjim/project/gargantua-studio/cmake/WriteSolarLock.cmake`
- Create: `/Users/hostsjim/project/gargantua-studio/cmake/solar-lock.cmake`
- Create: `/Users/hostsjim/project/gargantua-studio/src/probe_main.cpp`
- Create: `/Users/hostsjim/project/gargantua-studio/.gitignore`
- Create: `/Users/hostsjim/project/gargantua-studio/README.md`

**Interfaces:**
- Consumes: `Solar::Relativity`, `solar::version`, `solar::physics_contract`, and a locked Solar Git commit.
- Produces: executable `gargantua-probe` and CTest `gargantua.solar_probe`.

- [ ] **Step 1: Create the independent repository**

Run:

```bash
mkdir /Users/hostsjim/project/gargantua-studio
git -C /Users/hostsjim/project/gargantua-studio init -b main
```

Expected: Gargantua has its own `.git` directory and is not nested inside the Solar worktree.

- [ ] **Step 2: Add the deterministic Solar lock writer**

Create `cmake/WriteSolarLock.cmake`:

```cmake
if(
  NOT DEFINED SOLAR_SOURCE_DIR OR
  NOT DEFINED LOCK_OUTPUT
)
  message(
    FATAL_ERROR
    "SOLAR_SOURCE_DIR and LOCK_OUTPUT are required"
  )
endif()

execute_process(
  COMMAND git -C "${SOLAR_SOURCE_DIR}" rev-parse HEAD
  RESULT_VARIABLE solar_git_result
  OUTPUT_VARIABLE solar_git_commit
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
string(LENGTH "${solar_git_commit}" solar_commit_length)
if(
  NOT solar_git_result EQUAL 0 OR
  NOT solar_commit_length EQUAL 40 OR
  NOT solar_git_commit MATCHES "^[0-9a-f]+$"
)
  message(FATAL_ERROR "Solar source is not a valid Git commit")
endif()

file(
  WRITE "${LOCK_OUTPUT}"
  "set(GARGANTUA_SOLAR_REPOSITORY \"https://github.com/TT1nKer/solar.git\")\n"
  "set(GARGANTUA_SOLAR_COMMIT \"${solar_git_commit}\")\n"
  "set(GARGANTUA_SOLAR_VERSION \"0.2.0-alpha.1\")\n"
  "set(GARGANTUA_SOLAR_PHYSICS_CONTRACT \"relativity-v3-phase2\")\n"
)
```

- [ ] **Step 3: Generate and inspect the exact Solar lock**

Run:

```bash
cmake \
  -DSOLAR_SOURCE_DIR=/Users/hostsjim/project/solar \
  -DLOCK_OUTPUT=/Users/hostsjim/project/gargantua-studio/cmake/solar-lock.cmake \
  -P cmake/WriteSolarLock.cmake
sed -n '1,20p' cmake/solar-lock.cmake
```

Expected: `GARGANTUA_SOLAR_COMMIT` contains the exact 40-character Solar HEAD
after Tasks 1-3, with no branch name or symbolic ref.

- [ ] **Step 4: Add strict local-or-fetch dependency resolution**

Create `cmake/SolarDependency.cmake`:

```cmake
include(FetchContent)
include("${CMAKE_CURRENT_LIST_DIR}/solar-lock.cmake")

set(SOLAR_BUILD_CLI OFF CACHE BOOL "" FORCE)
set(
  GARGANTUA_SOLAR_SOURCE_DIR
  ""
  CACHE PATH
  "Optional local Solar checkout matching the locked commit"
)

if(GARGANTUA_SOLAR_SOURCE_DIR)
  execute_process(
    COMMAND
      git -C "${GARGANTUA_SOLAR_SOURCE_DIR}" rev-parse HEAD
    RESULT_VARIABLE solar_git_result
    OUTPUT_VARIABLE solar_git_commit
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if(
    NOT solar_git_result EQUAL 0 OR
    NOT solar_git_commit STREQUAL GARGANTUA_SOLAR_COMMIT
  )
    message(
      FATAL_ERROR
      "Local Solar checkout does not match cmake/solar-lock.cmake"
    )
  endif()
  add_subdirectory(
    "${GARGANTUA_SOLAR_SOURCE_DIR}"
    "${CMAKE_BINARY_DIR}/_deps/solar-build"
    EXCLUDE_FROM_ALL
  )
else()
  FetchContent_Declare(
    solar
    GIT_REPOSITORY "${GARGANTUA_SOLAR_REPOSITORY}"
    GIT_TAG "${GARGANTUA_SOLAR_COMMIT}"
    GIT_SHALLOW FALSE
  )
  FetchContent_MakeAvailable(solar)
endif()
```

- [ ] **Step 5: Add the first Gargantua CPU probe**

Create `src/probe_main.cpp`:

```cpp
#include "solar/relativity/kerr_shadow.h"
#include "solar/version.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string_view>
#include <vector>

int main() {
    using solar::relativity::KerrBoyerLindquistMetric;
    using solar::relativity::ShadowCriticalPoint;
    using solar::relativity::bardeen_shadow_curve;

    constexpr std::string_view required_version{
        "0.2.0-alpha.1"};
    constexpr std::string_view required_contract{
        "relativity-v3-phase2"};
    if (solar::version != required_version ||
        solar::physics_contract != required_contract) {
        std::cerr << "Solar package contract mismatch\n";
        return 1;
    }

    constexpr double half_pi = 1.5707963267948966;
    const KerrBoyerLindquistMetric metric(1.0, 0.5);
    const std::vector<ShadowCriticalPoint> curve =
        bardeen_shadow_curve(metric, half_pi, 65);
    const auto edges = std::minmax_element(
        curve.begin(),
        curve.end(),
        [](const ShadowCriticalPoint& left,
           const ShadowCriticalPoint& right) {
            return left.alpha < right.alpha;
        });
    if (curve.size() != 128 ||
        edges.first == curve.end() ||
        !std::isfinite(edges.first->alpha) ||
        !std::isfinite(edges.second->alpha)) {
        std::cerr << "Solar shadow probe failed\n";
        return 1;
    }

    std::cout << std::setprecision(17)
              << "{\"engine\":\"solar\","
              << "\"solar_version\":\"" << solar::version
              << "\",\"physics_contract\":\""
              << solar::physics_contract
              << "\",\"metric\":\"kerr-bl\","
              << "\"mass\":1,\"spin\":0.5,"
              << "\"samples\":" << curve.size()
              << ",\"left\":" << edges.first->alpha
              << ",\"right\":" << edges.second->alpha
              << "}\n";
    return 0;
}
```

- [ ] **Step 6: Add the Gargantua build and test**

Create `CMakeLists.txt`:

```cmake
cmake_minimum_required(VERSION 3.20)
project(GargantuaStudio VERSION 0.1.0 LANGUAGES CXX)

include(CTest)
include(cmake/SolarDependency.cmake)

add_executable(gargantua-probe src/probe_main.cpp)
target_compile_features(gargantua-probe PRIVATE cxx_std_17)
target_link_libraries(
  gargantua-probe
  PRIVATE Solar::Relativity
)

add_test(
  NAME gargantua.solar_probe
  COMMAND gargantua-probe
)
set_tests_properties(
  gargantua.solar_probe
  PROPERTIES
    PASS_REGULAR_EXPRESSION
      "\"physics_contract\":\"relativity-v3-phase2\""
)
```

- [ ] **Step 7: Add the repository README**

Create `.gitignore`:

```gitignore
/build/
```

Create `README.md`:

```markdown
# Gargantua Studio

Gargantua Studio is an offline Kerr black-hole film renderer powered and
validated by the independent Solar physics library.

This repository currently contains only the first cross-repository CPU probe.
It does not yet contain a renderer, CUDA backend, radiation model, or claim of
cinema-quality output.

## Build against the locked local Solar checkout

```bash
cmake -S . -B build \
  -DGARGANTUA_SOLAR_SOURCE_DIR=/Users/hostsjim/project/solar
cmake --build build
ctest --test-dir build --output-on-failure
./build/gargantua-probe
```

The local checkout must match `cmake/solar-lock.cmake`. Generated artifacts
must record the Solar version and physics contract.
```

- [ ] **Step 8: Verify wrong-checkout rejection**

Point the dependency option at the Gargantua repository itself:

```bash
cmake -S . \
  -B /tmp/gargantua-wrong-solar \
  -DGARGANTUA_SOLAR_SOURCE_DIR=/Users/hostsjim/project/gargantua-studio
```

Expected: configure fails with `Local Solar checkout does not match cmake/solar-lock.cmake`.

- [ ] **Step 9: Build and test the real local dependency**

Run:

```bash
cmake -S . \
  -B build \
  -DGARGANTUA_SOLAR_SOURCE_DIR=/Users/hostsjim/project/solar
cmake --build build --parallel
ctest --test-dir build --output-on-failure
./build/gargantua-probe
```

Expected: one CTest passes and JSON contains the required Solar version,
physics contract, 128 samples, and finite left/right shadow edges.

- [ ] **Step 10: Commit the independent repository**

```bash
git add \
  .gitignore \
  CMakeLists.txt \
  README.md \
  cmake/SolarDependency.cmake \
  cmake/WriteSolarLock.cmake \
  cmake/solar-lock.cmake \
  src/probe_main.cpp
git commit -m "feat: bootstrap Gargantua with Solar probe"
```

---

### Task 5: Cross-repository validation evidence

**Files:**
- Create: `/Users/hostsjim/project/gargantua-studio/docs/validation/00_solar_consumer.md`
- Modify: `/Users/hostsjim/project/gargantua-studio/README.md`

**Interfaces:**
- Consumes: successful package and probe results from both repositories.
- Produces: a reproducible validation report and an honest next boundary.

- [ ] **Step 1: Run fresh verification in Solar**

Run:

```bash
make test-external-consumer
make clean
make
make test
git diff --check
```

Expected: package acceptance and all existing Solar tests pass; the optional DE440 fixture may remain a visible skip.

- [ ] **Step 2: Run fresh verification in Gargantua**

Run:

```bash
cmake -S . \
  -B build \
  -DGARGANTUA_SOLAR_SOURCE_DIR=/Users/hostsjim/project/solar
cmake --build build --parallel
ctest --test-dir build --output-on-failure
./build/gargantua-probe
git diff --check
```

Expected: build succeeds, one CTest passes, and the probe prints Phase 2 Kerr edge evidence.

- [ ] **Step 3: Write the validation report from actual output**

Create `docs/validation/00_solar_consumer.md` with these exact headings:

```markdown
# Validation: Gargantua consumes Solar Phase 2

## Claim

## Model boundary

## Solar dependency

## Command

## Inputs

## Expected

## Actual

## Error

## Result

## Limitations

## Fastest falsification
```

Populate `Solar dependency` with the exact 40-character lock, package version,
and physics contract. Populate `Actual` only from the fresh Task 5 outputs.
State explicitly that this proves an external CPU library call, not CUDA,
radiative transfer, a renderer, or cinema-quality output.

- [ ] **Step 4: Link the evidence from the README**

Add under the build instructions:

```markdown
## Validation

The first external-library call is documented in
[`docs/validation/00_solar_consumer.md`](docs/validation/00_solar_consumer.md).
```

- [ ] **Step 5: Commit the evidence**

```bash
git add README.md docs/validation/00_solar_consumer.md
git commit -m "docs: validate Gargantua Solar dependency"
```

- [ ] **Step 6: Final responsibility and cleanliness audit**

In Solar, verify:

```bash
git status --short
git diff --check HEAD~3..HEAD
```

In Gargantua, verify:

```bash
git status --short
git diff --check HEAD~2..HEAD
```

Expected: no uncommitted project changes; Solar contains no Gargantua source or film dependency; Gargantua includes no Solar private header or copied implementation.
