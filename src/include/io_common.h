
#pragma once

#include <kg/io.h>

namespace
{

inline kg::io::Dims makeDims(int m, const Int3& dims)
{
  return kg::io::Dims{static_cast<size_t>(m), static_cast<size_t>(dims[2]),
                      static_cast<size_t>(dims[1]),
                      static_cast<size_t>(dims[0])};
}

inline kg::io::Dims makeDims(const Int3& dims)
{
  return kg::io::Dims{static_cast<size_t>(dims[2]),
                      static_cast<size_t>(dims[1]),
                      static_cast<size_t>(dims[0])};
}

template <typename E>
inline void write_4d(kg::io::Engine& file, const Int3& ldims, const Int3& gdims,
                     const std::vector<Int3>& patch_off, const E& h_expr)
{
  auto launch = kg::io::Mode::Blocking;

  int n_comps = h_expr.shape(3);
  int n_patches = h_expr.shape(4);
  Int3 im = {h_expr.shape(0), h_expr.shape(1), h_expr.shape(2)};
  Int3 ib = -(im - ldims) / 2;

  file.put("ib", ib, launch);
  file.put("im", im, launch);

  auto shape = makeDims(n_comps, gdims);
  for (int p = 0; p < n_patches; p++) {
    auto start = makeDims(0, patch_off[p]);
    auto count = makeDims(n_comps, ldims);
    auto _ib = makeDims(0, -ib);
    auto _im = makeDims(n_comps, im);
    file.putVariable(&h_expr(ib[0], ib[1], ib[2], 0, p), launch, shape,
                     {start, count}, {_ib, _im});
  }
  file.performPuts();
}

template <typename E>
inline void write_3d(kg::io::Engine& file, const Int3& ldims, const Int3& gdims,
                     const std::vector<Int3>& patch_off, const E& h_expr,
                     const std::string& name)
{
  // FIXME, should not pass name, but have prefixes set first, but gotta resolve
  // separator inconsistency first
  auto launch = kg::io::Mode::Blocking;

  int n_patches = h_expr.shape(3);
  Int3 im = {h_expr.shape(0), h_expr.shape(1), h_expr.shape(2)};
  Int3 ib = -(im - ldims) / 2;

  file.prefixes_.push_back(name);
  file.put("ib", ib, launch);
  file.put("im", im, launch);

  auto shape = makeDims(gdims);
  for (int p = 0; p < n_patches; p++) {
    auto start = makeDims(patch_off[p]);
    auto count = makeDims(ldims);
    auto _ib = makeDims(-ib);
    auto _im = makeDims(im);
    file.putVariable(&h_expr(ib[0], ib[1], ib[2], p), launch, shape,
                     {start, count}, {_ib, _im});
  }
  file.prefixes_.pop_back();
  // FIXME, it'd be better to write within the prefix, but xarray-adios2
  // expects a "/" separator instead of "::"
  file.put(std::string(name) + "/dimensions", std::string("z y x"));
  file.performPuts();
}

} // namespace
