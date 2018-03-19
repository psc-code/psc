
#pragma once

#include <cmath>

#include "pushp.hxx"

template<typename real_t>
class PI
{
public:
  using Real3 = Vec3<real_t>;
  
  PI(const Grid_t& grid)
    : dxi_{ Real3{ 1., 1., 1. } / Real3(grid.domain.dx) }
  {}
  
  // ----------------------------------------------------------------------
  // find_idx_off_1st_rel

  void find_idx_off_1st_rel(real_t xi[3], int lg[3], real_t og[3], real_t shift)
  {
    for (int d = 0; d < 3; d++) {
      real_t pos = xi[d] * dxi_[d] + shift;
      lg[d] = fint(pos);
      og[d] = pos - lg[d];
    }
  }

  // ----------------------------------------------------------------------
  // find_idx_off_pos_1st_rel

  void find_idx_off_pos_1st_rel(real_t xi[3], int lg[3], real_t og[3], real_t pos[3], real_t shift)
  {
    for (int d = 0; d < 3; d++) {
      pos[d] = xi[d] * dxi_[d] + shift;
      lg[d] = fint(pos[d]);
      og[d] = pos[d] - lg[d];
    }
  }

private:
  Real3 dxi_;
};

