
#pragma once

#include "bnd.hxx"

// ======================================================================
// BndCuda3
//
// just wrapping CudaBnd doing the actual work

struct CudaBnd;

template<typename MF>
struct BndCuda3 : BndBase
{
  using Mfields = MF;

  BndCuda3(const Grid_t& grid, int ibn[3]);
  ~BndCuda3();
  
  void reset(const Grid_t& grid);
  void add_ghosts(Mfields& mflds, int mb, int me);
  void fill_ghosts(Mfields& mflds, int mb, int me);

  void add_ghosts(MfieldsBase& mflds_base, int mb, int me) override
  {
    assert(0);
  }

  void fill_ghosts(MfieldsBase& mflds_base, int mb, int me) override
  {
    assert(0);
  }

private:
  CudaBnd* cbnd_;
  int balance_generation_cnt_;
};

