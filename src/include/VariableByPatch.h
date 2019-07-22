
#ifndef VARIABLE_BY_PATCH_H
#define VARIABLE_BY_PATCH_H

#include <kg/io.h>

// ======================================================================
// VariableByPatch

template <typename T>
struct VariableByPatch;

template <typename T>
struct VariableByPatch<std::vector<Vec3<T>>>
{
  using value_type = std::vector<Vec3<T>>;

  void put(kg::io::Engine& writer, const value_type& vec, const Grid_t& grid,
           const kg::io::Mode launch = kg::io::Mode::NonBlocking)
  {
    kg::io::Dims shape = {static_cast<size_t>(grid.nGlobalPatches()), 3};
    kg::io::Dims start = {
      static_cast<size_t>(grid.localPatchInfo(0).global_patch), 0};
    kg::io::Dims count = {static_cast<size_t>(grid.n_patches()), 3};
    writer.putVariable(vec[0].data(), launch, shape, {start, count});
  }

  void get(kg::io::Engine& reader, value_type& vec, const Grid_t& grid,
           const kg::io::Mode launch = kg::io::Mode::NonBlocking)
  {
    kg::io::Dims shape = {static_cast<size_t>(grid.nGlobalPatches()), 3};
    kg::io::Dims start = {
      static_cast<size_t>(grid.localPatchInfo(0).global_patch), 0};
    kg::io::Dims count = {static_cast<size_t>(grid.n_patches()), 3};
    assert(reader.variableShape<T>() == shape);
    vec.resize(count[0]);
    reader.getVariable(vec[0].data(), launch, {start, count});
  }
};

template <typename T>
struct VariableByPatch<std::vector<T>>
{
  using value_type = std::vector<T>;

  void put(kg::io::Engine& writer, const value_type& vec, const Grid_t& grid,
           const kg::io::Mode launch = kg::io::Mode::NonBlocking)
  {
    kg::io::Dims shape = {static_cast<size_t>(grid.nGlobalPatches())};
    kg::io::Dims start = {
      static_cast<size_t>(grid.localPatchInfo(0).global_patch)};
    kg::io::Dims count = {static_cast<size_t>(grid.n_patches())};
    writer.putVariable(vec.data(), launch, shape, {start, count});
  }

  void get(kg::io::Engine& reader, value_type& vec, const Grid_t& grid,
           const kg::io::Mode launch = kg::io::Mode::NonBlocking)
  {
    kg::io::Dims shape = {static_cast<size_t>(grid.nGlobalPatches())};
    kg::io::Dims start = {
      static_cast<size_t>(grid.localPatchInfo(0).global_patch)};
    kg::io::Dims count = {static_cast<size_t>(grid.n_patches())};
    assert(reader.variableShape<T>() == shape);
    vec.resize(count[0]);
    reader.getVariable(vec.data(), launch, {start, count});
  }
};

#endif
