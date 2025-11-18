#pragma once

#include <mrc_common.h>
#include <iostream>
#include "../libpsc/vpic/psc_vpic_bits.h"

namespace psc
{
namespace grid
{

// ======================================================================
// Domain

template <class R>
struct Domain
{
  using Real = R;
  using Real3 = Vec3<R>;

  Domain() {}

  Domain(Int3 gdims, Real3 length, Real3 corner = {0., 0., 0.},
         Int3 np = {1, 1, 1})
    : gdims(gdims), length(length), corner(corner), np(np)
  {
    for (int d = 0; d < 3; d++) {
      if (gdims[d] % np[d] != 0) {
        LOG_ERROR("in dimension %d, number of patches (%d) doesn't divide "
                  "number of cells (%d)\n",
                  d, np[d], gdims[d]);
      }

      if (length[d] <= 0.0) {
        LOG_ERROR("dimension %d has non-positive length (%f)\n", d, length[d]);
      }

      if (gdims[d] <= 0) {
        LOG_ERROR("dimension %d has non-positive number of cells (%d)\n", d,
                  gdims[d]);
      }
    }

    ldims = gdims / np;
    dx = length / Real3(gdims);
  }

  void view() const
  {
    mprintf("Grid_::Domain: gdims %d x %d x %d\n", gdims[0], gdims[1],
            gdims[2]);
  }

  bool isInvar(int d) const { return gdims[d] == 1; }

  /// @brief Total number of grid cells in each dimension.
  Int3 gdims;
  /// @brief Side lengths of the domain, in physical units.
  Real3 length;
  /// @brief Location of lower-right corner of the domain, in physical units.
  Real3 corner;
  /// @brief Number of patches in each dimension.
  Int3 np;
  /// @brief Number of grid cells per patch in each dimension.
  Int3 ldims;
  /// @brief Side lengths of each grid cell, in physical units.
  Real3 dx;
};

template <typename R>
inline std::ostream& operator<<(std::ostream& os, const Domain<R>& domain)
{
  os << "Domain{gdims=" << domain.gdims;
  os << ", length=" << domain.length;
  os << ", corner=" << domain.corner;
  os << ", np=" << domain.np;
  os << ", ldims=" << domain.ldims;
  os << ", dx=" << domain.dx;
  os << "}";
  return os;
}

} // namespace grid
} // namespace psc
