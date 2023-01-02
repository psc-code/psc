
#pragma once

#include "bnd.hxx"
#include "psc_fields_cuda.h"

// ======================================================================
// BndCuda3
//
// just wrapping CudaBnd doing the actual work

struct CudaBnd;

struct BndCuda3 : BndBase
{
  BndCuda3();
  ~BndCuda3();

  static void clear();

  template <typename MF>
  void add_ghosts(MF& mflds, int mb, int me);
  template <typename MF>
  void fill_ghosts(MF& mflds, int mb, int me);
  template <typename S>
  void add_ghosts(const Grid_t& grid, S& mflds_gt, const Int3& mflds_ib, int mb,
                  int me);
  template <typename S>
  void fill_ghosts(const Grid_t& grid, S& mflds_gt, const Int3& mflds_ib,
                   int mb, int me);

private:
  static CudaBnd* cbnd_;
};
