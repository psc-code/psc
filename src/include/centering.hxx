
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

/*  Centering Guide
         X   Y   Z
       -------------
NC  |   nnn nnn nnn
EC  |   cnn ncn nnc
FC  |   ncc cnc ccn
CC  |   ccc ccc ccc
*/

template <typename PATCH>
Double3 get_pos(PATCH patch, Int3 index, CenterStyle style, int component = -1)
{
  if (component == -1) {
    assert(style == NC || style == CC);
  }

  Double3 pos;
  for (int d = 0; d < 3; d++) {
    if (style == CC || (style == FC && d != component) ||
        (style == EC && d == component)) {
      pos[d] = patch.get_cc(index[d], d);
    } else {
      pos[d] = patch.get_nc(index[d], d);
    }
  }
  return pos;
}

struct Centerer
{
  CenterStyle style;
  Centerer(CenterStyle style) : style(style) {}

  template <typename PATCH>
  inline Double3 get_pos(PATCH patch, Int3 index, int component = -1) const
  {
    return centering::get_pos(patch, index, style, component);
  }
};

} // namespace centering