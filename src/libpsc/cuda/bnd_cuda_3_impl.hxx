
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
  using Mfields = MF;
  using storage_type = typename Mfields::Storage;

  BndCuda3(const Grid_t& grid, const int ibn[3]);
  ~BndCuda3();

  static void clear();

  void add_ghosts(Mfields& mflds, int mb, int me);
  void fill_ghosts(Mfields& mflds, int mb, int me);
  void add_ghosts(const Grid_t& grid, storage_type& mflds_gt,
                  const Int3& mflds_ib, int mb, int me);
  void fill_ghosts(const Grid_t& grid, storage_type& mflds_gt,
                   const Int3& mflds_ib, int mb, int me);

private:
  static CudaBnd* cbnd_;
};
