#pragma once

#include "../kg/include/kg/Vec3.h"
#include "grid.hxx"

template <typename FE>
void add_ghosts_reflecting_lo_cc(const Int3& ldims, FE& mres_gt, const Int3& ib,
                                 int p, int d, int mb, int me)
{
  // FIXME only need to scan 1 cell into the ghost region for 1st-order moments

  Int3 idx_begin = ib;
  Int3 idx_end = ldims - ib;

  idx_begin[d] = 0;
  idx_end[d] = -ib[d];

  Int3 idx;

  for (int m = mb; m < me; m++) {
    for (idx[2] = idx_begin[2]; idx[2] < idx_end[2]; idx[2]++) {
      for (idx[1] = idx_begin[1]; idx[1] < idx_end[1]; idx[1]++) {
        for (idx[0] = idx_begin[0]; idx[0] < idx_end[0]; idx[0]++) {
          Int3 idx_reflected = idx;
          idx_reflected[d] = -idx[d] - 1;

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

template <typename FE>
void add_ghosts_reflecting_hi_cc(const Int3& ldims, FE& mres_gt, const Int3& ib,
                                 int p, int d, int mb, int me)
{
  // FIXME only need to scan 1 cell into the ghost region for 1st-order moments

  Int3 idx_begin = ib;
  Int3 idx_end = ldims - ib;

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

template <typename FE>
void add_ghosts_reflecting_lo_nc(const Int3& ldims, FE& mres_gt, const Int3& ib,
                                 int p, int d, int mb, int me)
{
  // When lower boundary has 2 ghosts, upper boundary of nc grid only has 1
  // because nc has +1 point compared to cc, and takes it from upper ghost
  // region. Thus, ignore a ghost from lower ghost region to balance it out.

  Int3 unused_ghost = {!!ib[0], !!ib[1], !!ib[2]};

  Int3 idx_begin = ib + unused_ghost;
  Int3 idx_end = ldims - ib;

  // Also note that the boundary nc point is the point about which reflections
  // occur, and isn't itself reflected.

  idx_begin[d] = 1;
  idx_end[d] = 1 - ib[d] - unused_ghost[d];

  Int3 idx;

  for (int m = mb; m < me; m++) {
    for (idx[2] = idx_begin[2]; idx[2] < idx_end[2]; idx[2]++) {
      for (idx[1] = idx_begin[1]; idx[1] < idx_end[1]; idx[1]++) {
        for (idx[0] = idx_begin[0]; idx[0] < idx_end[0]; idx[0]++) {
          Int3 idx_reflected = idx;
          idx_reflected[d] = -idx[d];

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

template <typename FE>
void add_ghosts_reflecting_hi_nc(const Int3& ldims, FE& mres_gt, const Int3& ib,
                                 int p, int d, int mb, int me)
{
  // When lower boundary has 2 ghosts, upper boundary of nc grid only has 1
  // because nc has +1 point compared to cc, and takes it from upper ghost
  // region. Thus, ignore a ghost from lower ghost region to balance it out.

  Int3 unused_ghost = {!!ib[0], !!ib[1], !!ib[2]};

  Int3 idx_begin = ib + unused_ghost;
  Int3 idx_end = ldims - ib;

  // Also note that the extra nc point is the point about which reflections
  // occur, and isn't itself reflected.

  idx_begin[d] = ldims[d] + ib[d] + unused_ghost[d];
  idx_end[d] = ldims[d];

  Int3 idx;

  for (int m = mb; m < me; m++) {
    for (idx[2] = idx_begin[2]; idx[2] < idx_end[2]; idx[2]++) {
      for (idx[1] = idx_begin[1]; idx[1] < idx_end[1]; idx[1]++) {
        for (idx[0] = idx_begin[0]; idx[0] < idx_end[0]; idx[0]++) {
          Int3 idx_reflected = idx;
          idx_reflected[d] = 2 * ldims[d] - idx[d];

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
