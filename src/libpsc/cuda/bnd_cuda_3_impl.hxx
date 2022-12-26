
#pragma once

#include "bnd.hxx"
#include "psc_fields_cuda.h"

// ======================================================================
// BndCuda3
//
// just wrapping CudaBnd doing the actual work

struct CudaBnd;

template <typename MF>
struct BndCuda3 : BndBase
{
  BndCuda3(const Grid_t& grid, const int ibn[3]);
  ~BndCuda3();

  static void clear();

  template <typename MF2>
  void add_ghosts(MF2& mflds, int mb, int me);
  template <typename MF2>
  void fill_ghosts(MF2& mflds, int mb, int me);
  template <typename S>
  void add_ghosts(const Grid_t& grid, S& mflds_gt, const Int3& mflds_ib, int mb,
                  int me);
  template <typename S>
  void fill_ghosts(const Grid_t& grid, S& mflds_gt, const Int3& mflds_ib,
                   int mb, int me);

private:
  static CudaBnd* cbnd_;
};
