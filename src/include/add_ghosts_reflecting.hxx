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
  int bx = ldims[0] == 1 ? 0 : 1;
  if (d == 1) {
    for (int iz = -1; iz < ldims[2] + 1; iz++) {
      for (int ix = -bx; ix < ldims[0] + bx; ix++) {
        int iy = ldims[1] - 1;
        {
          for (int m = mb; m < me; m++) {
            mres_gt(ix - ib[0], iy - ib[1], iz - ib[2], m, p) +=
              mres_gt(ix - ib[0], iy - 1 + ib[1], iz - ib[2], m, p);
          }
        }
      }
    }
  } else if (d == 2) {
    for (int iy = 0 * -1; iy < ldims[1] + 0 * 1; iy++) {
      for (int ix = -bx; ix < ldims[0] + bx; ix++) {
        int iz = ldims[2] - 1;
        {
          for (int m = mb; m < me; m++) {
            mres_gt(ix - ib[0], iy - ib[1], iz - ib[2], m, p) +=
              mres_gt(ix - ib[0], iy - ib[1], iz + 1 - ib[2], m, p);
          }
        }
      }
    }
  } else {
    assert(0);
  }
}