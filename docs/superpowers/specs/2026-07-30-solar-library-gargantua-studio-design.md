# Solar Library and Gargantua Studio Program Design

**Status:** Approved direction; written design pending user review
**Date:** 2026-07-30
**Source contract:** `SOLAR_RELATIVITY_KERR_完整实施主提示词_v3.md`

## 1. Program goal

This program has two independent repositories with a one-way dependency:

```text
Solar Library
    |
    | versioned public API, validation fixtures, CPU reference results
    v
Gargantua Studio
    |
    | offline rendering, film pipeline, shot assets
    v
4K master and, later, a reduced real-time web experience
```

Solar becomes a reusable orbital, relativistic, and numerical library.
Gargantua Studio is a real downstream consumer whose first product is a
cinema-quality Kerr black-hole sequence. The film is both a portfolio artifact
and a forcing function for better Solar APIs. Solar must never depend on
Gargantua, CUDA, film assets, or a DCC application.

The initial public brand is:

```text
GARGANTUA — A Kerr Black Hole Powered by Solar
```

Repository and product names may change without changing this architecture.

## 2. Why two repositories

The separation creates an honest product boundary:

- Solar owns reusable physics, numerical contracts, CPU reference behavior,
  public headers, validation reports, and releases.
- Gargantua owns CUDA rendering, scenes, cameras, artistic emission,
  production assets, EXR frames, encoding, and shot-specific choices.
- Gargantua consumes a pinned Solar release or commit. It never includes
  Solar private headers or reaches into `src/`.
- Solar may advertise Gargantua as a downstream project; Gargantua links to
  the exact Solar version used for every render.

Copying Solar into the film repository, using filesystem-relative source
paths, or letting a render preset alter Solar physics are prohibited.

## 3. Solar as a consumable library

Solar retains the existing C++17 library and zero-dependency CPU reference
core. The existing Makefile remains supported while an installable CMake
package is added incrementally.

The first external-consumption contract is:

```cmake
find_package(Solar CONFIG REQUIRED)
target_link_libraries(gargantua PRIVATE Solar::Relativity)
```

During pre-release development, Gargantua may use CMake `FetchContent` pinned
to an immutable Solar commit. Public film manifests must use a Solar release
tag. Solar will provide:

- exported `Solar::Core` and `Solar::Relativity` targets;
- installable public headers and package configuration;
- semantic pre-release versions;
- a machine-readable version and physics-contract identifier;
- deterministic validation fixtures and reference outputs;
- an explicit compatibility policy for public types and serialized data.

C ABI and WASM packages belong to the later web track. CUDA is not added as a
mandatory Solar dependency. The first GPU implementation remains a
Gargantua backend and is accepted only through comparison with Solar's CPU
reference.

## 4. Gargantua Studio product boundary

Gargantua Studio is an independent C++/CUDA application and focused film
toolkit, not a second general-purpose physics engine and not a replacement
for Blender or Houdini.

Its first modules are:

```text
scene-contract       validated scene, camera, timeline, and render settings
solar-adapter        conversion to versioned Solar public types
cpu-reference        Solar-backed reference rays and diagnostic probes
cuda-renderer        RTX 3080 high-throughput image backend
radiation            disk/fluid/emission models and invariant transfer
render-scheduler     tiles, frames, checkpoints, retries, and manifests
film-output          OpenEXR AOVs, OCIO/ACES transforms, and encode handoff
gargantua-cli        validate, probe, render, resume, inspect, and encode
```

The dependency direction is:

```text
CLI -> shot flow -> scheduler -> renderer interface
                              -> CPU Solar adapter
                              -> CUDA backend
renderer -> radiation models -> scene-contract
film-output -> completed frame data
```

Scene and output modules do not call CUDA. Solar-facing code does not know
about film presets. CUDA kernels do not parse files or own the render queue.

## 5. Build versus reuse

Gargantua implements only the parts that create scientific and product value:

- Kerr/Schwarzschild ray propagation for image workloads;
- invariant radiative transfer and redshift;
- disk crossings and volumetric integration;
- critical-region adaptive sampling;
- temporal stability and deterministic sampling;
- CPU/GPU comparison and per-pixel physics diagnostics;
- resumable shot rendering.

It reuses established film infrastructure:

- OpenEXR for scene-linear HDR, tiled frames, metadata, and AOV channels;
- OpenColorIO with an ACES working/output configuration;
- FFmpeg or a DCC application for final mezzanine and delivery encoding;
- later, USD or a DCC bridge for ordinary production geometry.

OptiX or RT cores may later accelerate conventional mesh intersection. They
are not treated as a curved-spacetime geodesic integrator. Version numbers are
pinned only after probing the actual RTX 3080 driver/toolkit combination.

## 6. Physics and artistic modes

Geometry and ray propagation always use validated physics. Radiation is
explicitly classified:

- `SCIENCE`: physical model, physical parameters, full diagnostic output;
- `ART_DIRECTED_RADIATION`: physical geodesics with documented artistic
  emission, flare, exposure, or frequency compression;
- `ENGINE_DEBUG`: no bloom, flare, grain, vignette, or hidden repair.

Every frame manifest records the mode, Solar version, Gargantua version,
metric, chart, mass, spin, observer, integrator/backend, tolerances, scene
hash, sample counts, unconverged pixels, and output checksums. The film may be
art directed, but no artistic choice may be represented as a physical result.

## 7. RTX 3080 execution model

The RTX 3080 is the first production render worker. The local Apple Silicon
machine remains the development host and CPU-reference environment.

The renderer is designed for bounded GPU memory:

- tile-based execution with configurable tile size;
- no whole-sequence or full AOV set resident on the GPU;
- atomic tile/frame completion records;
- deterministic restart from a manifest;
- bounded queues and explicit allocation failures;
- per-tile non-finite and unconverged masks;
- CPU recheck of selected critical pixels.

GPU precision is a measured policy, not an assumption. Ordinary rays may use
mixed precision only after CPU comparison passes. Critical, near-turning, and
failed rays are retried with a stricter path or classified as unconverged.
They are never silently painted black.

The production-time target is at most 48 hours for the first 288-frame
sequence on the RTX 3080. This is a performance target and may not weaken
physics or convergence gates.

## 8. First film acceptance target

The first release is one continuous exterior Kerr hero shot:

- 3840 by 2160;
- 24 frames per second;
- approximately 12 seconds, 288 frames;
- black hole, physical lensing, star environment, and an emitting disk;
- no spacecraft, dialogue, editor UI, or multi-shot asset system;
- scene-linear multi-channel OpenEXR frame sequence;
- ACES-managed 4K master plus a web delivery encode;
- deterministic render manifest and reproducible command;
- a validation report linking the rendered shot to Solar CPU fixtures.

Beauty channels may use half precision where the film pipeline supports it.
Physics diagnostics such as constraint error and conserved quantities remain
32-bit channels or sidecar data. Required AOVs are:

```text
beauty.RGBA
emission.RGB
redshift
min_radius
winding
disk_crossing_count
constraint_error
sample_count
classification
unconverged
```

The final sequence contains zero silently unconverged pixels. Any unresolved
pixel remains visible in diagnostics and blocks the master until rerendered or
explicitly waived in the validation report.

## 9. Frame data flow

```text
versioned scene file
  -> schema and resource validation
  -> immutable shot manifest
  -> camera state for frame and temporal sample
  -> tile plan
  -> CPU reference probes or CUDA ray batches
  -> geodesic events and radiation integration
  -> beauty plus physics AOVs
  -> atomic tiled EXR
  -> frame validation and checksums
  -> OCIO/ACES review transform
  -> master and web encode
```

The scene contract uses explicit units and rejects unknown enum values,
non-finite numbers, invalid paths, impossible camera states, unsupported
Solar contracts, and resource limits before allocating a render.

## 10. Failure and recovery contract

Every failure is one of:

```text
invalid_scene
incompatible_solar_contract
resource_exhausted
gpu_unavailable
ray_unconverged
physics_constraint_violation
tile_write_failed
frame_incomplete
encode_failed
```

Completed tiles are immutable and checksummed. Partial files use a temporary
name and become visible only after a successful close and atomic rename.
Resume verifies the scene hash, software versions, device policy, and every
completed tile checksum. A changed scene or physics version starts a new
render generation instead of mixing frames.

## 11. Verification gates

Solar remains the authority for CPU equations and validation. Gargantua adds:

- CPU/GPU comparisons for ordinary, near-critical, spin-sign, observer, and
  disk rays;
- Schwarzschild and Bardeen shadow boundary checks from the v3 contract;
- Hamiltonian and Carter diagnostics on selected production rays;
- deterministic same-device rerender checks;
- tile-order and resume equivalence;
- EXR channel/metadata round trips;
- malformed scene and resource-limit tests;
- static-camera temporal stability and sample-convergence reports;
- one higher-quality audit render for selected production frames.

The exact v3 error gates remain hard requirements. Image quality settings may
increase work but may not relabel unconverged physics as success.

## 12. Delivery sequence

This program is too large for one implementation plan. Work proceeds through
separate reviewed specs and plans:

1. **Solar package boundary:** CMake export/install, version contract, and a
   minimal external consumer test.
2. **Solar physics continuation:** Phase 3 separated Kerr solver, followed by
   the remaining v3 CPU validation needed by rendering.
3. **Gargantua repository foundation:** independent repository, pinned Solar
   dependency, scene contract, CPU still-image probe, and render manifest.
4. **CUDA backend:** RTX 3080 worker, tiled scheduler, CPU/GPU comparison,
   checkpointing, and diagnostic buffers.
5. **Radiation and film output:** invariant transfer, emitting disk, OpenEXR,
   OCIO/ACES, and convergence reports.
6. **Hero sequence:** camera, look development, 288-frame render, validation,
   master, and website delivery encode.
7. **Web track:** reduced interactive renderer derived from validated assets
   and contracts, not a replacement physics implementation.

The next implementation plan covers only item 1, the Solar package boundary.
The Gargantua repository is not scaffolded until that external-consumer
contract passes.

## 13. Publication and licensing boundary

No license is changed by this design. The new repository remains private
until the owner explicitly selects publication and licensing terms. Before a
public release:

- Solar must have a clear library license and versioned release;
- Gargantua must declare its own code and asset licenses;
- third-party binary, texture, sky, font, and audio rights must be recorded;
- every public film must identify its Solar and Gargantua versions.

This prevents the desire for propagation from silently making legal or asset
licensing decisions.

## 14. Authoritative references

- Project v3 contract, especially sections 13, 20.4, 24, 25, 27, and 28.
- CUDA C++ Programming Guide:
  <https://docs.nvidia.com/cuda/cuda-c-programming-guide/>.
- OpenEXR documentation:
  <https://openexr.com/en/latest/index.html>.
- OpenColorIO documentation:
  <https://opencolorio.org/>.

