# Solar external consumer contract

Solar `0.2.0-alpha.1` is the first package contract intended for an
independent downstream project. It exposes two CMake targets:

- `Solar::Core` for the existing astrodynamics library.
- `Solar::Relativity` for metrics, observers, geodesics, Kerr invariants, and
  shadow geometry. It transitively links `Solar::Core`.

The current physics contract identifier is `relativity-v3-phase2`. Downstream
projects should record both `solar::version` and `solar::physics_contract` in
validation artifacts. A changed physics contract means numerical results must
be revalidated even when source compatibility is preserved.

## Install and consume

```sh
cmake -S /path/to/solar -B /tmp/solar-build \
  -DSOLAR_BUILD_CLI=OFF \
  -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/solar-build --parallel
cmake --install /tmp/solar-build --prefix /tmp/solar-install
```

An independent CMake project can then use:

```cmake
find_package(Solar CONFIG REQUIRED)

add_executable(my_renderer main.cpp)
target_compile_features(my_renderer PRIVATE cxx_std_17)
target_link_libraries(my_renderer PRIVATE Solar::Relativity)
```

Configure that project with:

```sh
cmake -S . -B build -DCMAKE_PREFIX_PATH=/tmp/solar-install
cmake --build build --parallel
```

The full prerelease identifier is defined in `solar/version.h`. The generated
CMake compatibility file uses the numeric base version `0.2.0`, because CMake
package version comparison does not model SemVer prerelease precedence.

## Acceptance boundary

`make test-external-consumer` builds and installs Solar into an isolated
temporary prefix, configures a separate project through `find_package`, and
runs a Kerr shadow probe. The probe locks:

- the public version and physics-contract identifiers;
- 128 closed-curve samples for an equatorial Kerr black hole with
  geometrized mass `M=1` and dimensionless spin `chi=0.5`;
- the left and right critical-curve edges to a tolerance of `1e-13`.

This is a packaging and public-physics smoke test. It complements, rather than
replaces, Solar's focused numerical and invariant test suites.
