#pragma once

#include "../kg/include/kg/Vec3.h"
#include "grid.hxx"

template <typename FE>
void add_ghosts_reflecting_lo(const Int3& ldims, FE& mres_gt, const Int3& ib,
                              int p, int d, int mb, int me)
{
  int bx = ldims[0] == 1 ? 0 : 1;
  if (d == 1) {
    for (int iz = -1; iz < ldims[2] + 1; iz++) {
      for (int ix = -bx; ix < ldims[0] + bx; ix++) {
        int iy = 0;
        {
          for (int m = mb; m < me; m++) {
            mres_gt(ix - ib[0], iy - ib[1], iz - ib[2], m, p) +=
              mres_gt(ix - ib[0], iy - 1 - ib[1], iz - ib[2], m, p);
          }
        }
      }
    }
  } else if (d == 2) {
    for (int iy = 0 * -1; iy < ldims[1] + 0 * 1; iy++) {
      for (int ix = -bx; ix < ldims[0] + bx; ix++) {
        int iz = 0;
        {
          for (int m = mb; m < me; m++) {
            mres_gt(ix - ib[0], iy - ib[1], iz - ib[2], m, p) +=
              mres_gt(ix - ib[0], iy - ib[1], iz - 1 - ib[2], m, p);
          }
        }
      }
    }
  } else {
    assert(0);
  }
}

template <typename FE>
void add_ghosts_reflecting_hi(const Int3& ldims, FE& mres_gt, const Int3& ib,
                              int p, int d, int mb, int me)
{
  // FIXME only need to scan 1 cell into the ghost region for 1st-order moments

  Int3 idx_end = ldims - ib;
  Int3 idx_begin = ib;

  idx_begin[d] = ldims[d] + ib[d];
  idx_end[d] = ldims[d];

  Int3 idx;

  for (int m = mb; m < me; m++) {
    for (idx[2] = idx_begin[2]; idx[2] < idx_end[2]; idx[2]++) {
      for (idx[1] = idx_begin[1]; idx[1] < idx_end[1]; idx[1]++) {
        for (idx[0] = idx_begin[0]; idx[0] < idx_end[0]; idx[0]++) {
          Int3 idx_reflected = idx;
          idx_reflected[d] = 2 * ldims[d] - idx[d] - 1;

          auto val_reflected =
            mres_gt(idx_reflected[0] - ib[0], idx_reflected[1] - ib[1],
                    idx_reflected[2] - ib[2], m, p);

          mres_gt(idx[0] - ib[0], idx[1] - ib[1], idx[2] - ib[2], m, p) +=
            val_reflected;
        }
      }
    }
  }
}