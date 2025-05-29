
#pragma once

#include <cmath>

#include "pushp.hxx"

template <typename real_t>
class PI
{
public:
  using Real3 = Vec3<real_t>;

  PI(const Grid_t& grid) : dxi_{Real3{1., 1., 1.} / Real3(grid.domain.dx)} {}

  // ----------------------------------------------------------------------
  // find_idx_off_pos_1st_rel

  void find_idx_off_pos_1st_rel(Real3 final_pos, Int3 final_index,
                                Real3 final_offset, Real3 final_pos_normalized,
                                real_t shift)
  {
    for (int d = 0; d < 3; d++) {
      final_pos_normalized[d] = final_pos[d] * dxi_[d] + shift;
      final_index[d] = fint(final_pos_normalized[d]);
      final_offset[d] = final_pos_normalized[d] - final_index[d];
    }
  }

private:
  Real3 dxi_;
};
