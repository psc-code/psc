
#pragma once

#include "kg/Vec3.h"

namespace centering
{

enum CenterStyle
{
  NC, // node-centered
  EC, // edge-centered
  FC, // face-centered
  CC, // cell-centered
};

enum Component
{
  X,
  Y,
  Z
};

/*  Centering Guide
         X   Y   Z
       -------------
NC  |   nnn nnn nnn
EC  |   cnn ncn nnc
FC  |   ncc cnc ccn
CC  |   ccc ccc ccc
*/

template <typename PATCH>
Double3 getPos(PATCH patch, Int3 index, int style, int comp = X)
{
  Double3 pos;
  for (int a = 0; a < 3; a++) {
    if (style == CC || (style == FC && a != comp) ||
        (style == EC && a == comp)) {
      pos[a] = patch.get_cc(index[a], a);
    } else {
      pos[a] = patch.get_nc(index[a], a);
    }
  }
  return pos;
}

struct Centerer
{
  CenterStyle style;
  Centerer(CenterStyle style) : style(style) {}

  template <typename PATCH>
  inline Double3 getPos(PATCH patch, Int3 index, int comp = X) const
  {
    return centering::getPos(patch, index, style, comp);
  }
};

} // namespace centering