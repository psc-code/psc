# PSC

The Plasma Simulation Code: a 3D fully electromagnetic particle-in-cell code for high-performance kinetic plasma simulations, used in academic research. Production runs are generally on Linux supercomputer clusters, so that is the platform to target. It began in Fortran, mainly for laser-plasma interactions, and was ported to C++ to get GPU support. It uses CMake and gtensor, with optional CUDA/HIP. `psc` itself only requires `cxx_std_11`, but gtensor pulls in C++17, so in practice it builds as C++17.

## Priorities

When these conflict, the earlier one wins:

0. **Physical accuracy.** A technique must be correct to its stated order, within the expected error. A lower-order method is fine; an incorrect one is worthless.
1. **Performance.** Scientists need ever-larger runs, and performance is the most important metric by which codes are compared.
2. **Extensibility.** Scientists need to adapt the code to test new physics.
3. **Usability.** User-facing code (cases, parameters, setup) should be well documented and avoid sophisticated C++ such as metaprogramming. Installing and running from scratch should be simple.
4. **Consistency.** See [Naming](#naming).

## Commands

```bash
# Build (assumes an already-configured ./build)
cmake --build build -j                                  # everything
cmake --build build -t test_neumann_boundary_injector   # one test binary
./build/src/libpsc/tests/test_neumann_boundary_injector --gtest_filter='*.InwardsLo'
cd build && ctest -R NeumannBoundaryInjectorTest        # ctest names = gtest names, not binary names

# What CI runs (.github/workflows/psc-cpu.yml, GCC 13 in ghcr.io/psc-code/psc-cpp-ubuntu-20.04)
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo -G Ninja && cmake --build build -t test
```

GPU backend is chosen with `-DPSC_GPU=host|cuda|hip` (README's `USE_CUDA` is legacy).
Formatting: clang-format 14, `.clang-format` (Mozilla-based); CI checks `src/` on PRs, except paths in `.clang-format-ignore`.

## Layout

- `src/psc_*.cxx`: cases ("input decks"); setup is hardcoded or read from a parameter file, with almost no CLI options. Registered with `add_psc_executable` in `src/CMakeLists.txt` (some only when `NOT USE_CUDA`, e.g. `psc_shock`).
- `src/include/`: core headers (fields, particles, grid, push/bnd interfaces).
- `src/libpsc/<component>/`: implementations (`psc_push_particles`, `psc_bnd_*`, `psc_particle_injectors`, `cuda/`, ...).
- `src/libpsc/tests/`: gtest tests, one `test_*.cxx` per binary, registered with `add_psc_test(name)` in that dir's `CMakeLists.txt`.
- `src/kg/`, `src/libmrc/`: support libraries (libmrc is legacy C, not clang-formatted).
- `python/`: reading output / analysis.

## State of the codebase

Inconsistencies are everywhere, so verify assumptions against the code instead of generalizing from one place. Template metaprogramming is common.

Migrations in progress (there may be others). New code should follow the target direction:

- `src/include/` → `src/libpsc/`, with namespacing.
- The integrator in `psc.hxx` is moving from compile-time templating on components to runtime polymorphism. The underlying data types stay templated for performance; virtual dispatch between components is cheap, and the templating is inflexible.
- Short names → longer descriptive ones, especially in interpolation code.
- `libmrc` ("mrc" = magnetic reconnection) is being phased out. It was an external library absorbed for domain decomposition. It is C, with macros that emulate classes.
- The conversion functions between container types are mostly legacy. Some are still needed, e.g. GPU particles ↔ CPU particles.
- Cases are moving toward parameter files (see `psc_bgk.cxx`, `psc_shock.cxx`) instead of fully hardcoded setups.
- Parts of VPIC (another PIC code) used to be included; they have mostly been removed.

Planned, longer-term: finish the migrations above; make `Mparticles` struct-of-arrays instead of array-of-structs; simplify `MfieldsState`/`Mfields`; formally move to C++17 or newer; improve the test framework, including integration tests for parameter-file cases; GPU support via Kokkos; general cleanup.

## Naming

These short names are standard throughout. Elsewhere, prefer longer descriptive names. Add to this list as more conventions are identified.

- `p`: patch index (in loops). `d`: dimension index. `m`: component index of a vector field.
- `ib` / `ie` / `im`: index begin / end / extent (`ie - ib`).
- `ldims` / `gdims`: local (per-patch) / global grid dimensions.
- `CC` / `NC` / `FC` / `EC`: cell-, node-, face-, edge-centered relative to the Yee cell.
