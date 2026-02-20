
#pragma once

#include <type_traits>

#include "kg/Vec3.h"

template <bool IX = false, bool IY = false, bool IZ = false>
struct Invar
{
  using InvarX = std::integral_constant<bool, IX>;
  using InvarY = std::integral_constant<bool, IY>;
  using InvarZ = std::integral_constant<bool, IZ>;

  /// @brief Gets a mask of noninvariant dimensions (i.e., 0 for invariant
  /// dimensions, 1 for noninvariant dimensions). This is useful for multiplying
  /// by a constant to, say, get the number of ghost cells in each dimension
  /// (where invariant dimensions have no ghosts).
  /// @return the mask
  static Int3 get_noninvariant_mask()
  {
    return {!InvarX::value, !InvarY::value, !InvarZ::value};
  }

  static bool is_invar(int dim)
  {
    switch (dim) {
      case 0: return InvarX::value;
      case 1: return InvarY::value;
      case 2: return InvarZ::value;
      default: return false;
    }
  }
};

using dim_xyz = Invar<false, false, false>;
using dim_xy = Invar<false, false, true>;
using dim_xz = Invar<false, true, false>;
using dim_yz = Invar<true, false, false>;
using dim_x = Invar<false, true, true>;
using dim_y = Invar<true, false, true>;
using dim_z = Invar<true, true, false>;
using dim_1 = Invar<true, true, true>;
