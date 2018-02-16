
#pragma once

#include <type_traits>

template<bool IX = false, bool IY = false, bool IZ = false>
struct Invar
{
  using InvarX = std::integral_constant<bool, IX>;
  using InvarY = std::integral_constant<bool, IY>;
  using InvarZ = std::integral_constant<bool, IZ>;
};

using dim_xyz = Invar<false, false, false>;
using dim_xy  = Invar<false, false, true >;
using dim_xz  = Invar<false, true , false>;
using dim_yz  = Invar<true , false, false>;
using dim_x   = Invar<false, true , true >;
using dim_y   = Invar<true , false, true >;
using dim_z   = Invar<true , true , false>;
using dim_1   = Invar<true , true , true >;

