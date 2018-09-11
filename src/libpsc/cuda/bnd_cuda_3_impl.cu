
#include "bnd_cuda_3_impl.hxx"
#include "cuda_bnd.cuh"
#include "balance.hxx"

// ----------------------------------------------------------------------
// ctor

template<typename MF>
BndCuda3<MF>::BndCuda3(const Grid_t& grid, int ibn[3])
  : cbnd_{new CudaBnd{grid, ibn}},
    balance_generation_cnt_{psc_balance_generation_cnt}
{}

// ----------------------------------------------------------------------
// dtor

template<typename MF>
BndCuda3<MF>::~BndCuda3()
{
  delete cbnd_;
}
  
// ----------------------------------------------------------------------
// reset

template<typename MF>
void BndCuda3<MF>::reset(const Grid_t& grid)
{
  // FIXME, not really a pretty way of doing this
  delete cbnd_;
  cbnd_ = new CudaBnd{grid, grid.ibn};
}
  
// ----------------------------------------------------------------------
// add_ghosts

template<typename MF>
void BndCuda3<MF>::add_ghosts(Mfields& mflds, int mb, int me)
{
  if (psc_balance_generation_cnt != balance_generation_cnt_) {
    reset(mflds.grid());
    balance_generation_cnt_ = psc_balance_generation_cnt;
  }
  cbnd_->add_ghosts(*mflds.cmflds(), mb, me);
}

// ----------------------------------------------------------------------
// fill_ghosts

template<typename MF>
void BndCuda3<MF>::fill_ghosts(Mfields& mflds, int mb, int me)
{
  if (psc_balance_generation_cnt != balance_generation_cnt_) {
    reset(mflds.grid());
    balance_generation_cnt_ = psc_balance_generation_cnt;
  }
  cbnd_->fill_ghosts(*mflds.cmflds(), mb, me);
}

template struct BndCuda3<MfieldsCuda>;
template struct BndCuda3<MfieldsStateCuda>;
