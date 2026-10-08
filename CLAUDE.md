# PSC

3D electromagnetic particle-in-cell code (C++17, CMake, gtensor; optional CUDA/HIP).

## Commands

```bash
# Local build (existing ./build: Release, Unix Makefiles, PSC_GPU=host, Apple Clang)
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

- `src/psc_*.cxx`: cases ("input decks"); everything is hardcoded, there are almost no CLI options. Registered with `add_psc_executable` in `src/CMakeLists.txt` (some only when `NOT USE_CUDA`, e.g. `psc_shock`).
- `src/include/`: core headers (fields, particles, grid, push/bnd interfaces).
- `src/libpsc/<component>/`: implementations (`psc_push_particles`, `psc_bnd_*`, `psc_particle_injectors`, `cuda/`, ...).
- `src/libpsc/tests/`: gtest tests, one `test_*.cxx` per binary, registered with `add_psc_test(name)` in that dir's `CMakeLists.txt`.
- `src/kg/`, `src/libmrc/`: support libraries (libmrc is legacy C, not clang-formatted).
- `python/`: reading output / analysis.
