
#pragma once

#include <cmath>

#include "pushp.hxx"

template <typename real_t, typename Real3 = Vec3<real_t>>
void find_idx_pos_1st_rel(Real3 final_pos, Real3 dxi, Int3 final_index,
                          Real3 final_pos_normalized, real_t shift)
{
  for (int d = 0; d < 3; d++) {
    final_pos_normalized[d] = final_pos[d] * dxi[d] + shift;
    final_index[d] = fint(final_pos_normalized[d]);
  }
}
