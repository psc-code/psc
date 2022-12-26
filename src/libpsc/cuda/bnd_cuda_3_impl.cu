
#include "bnd_cuda_3_impl.hxx"
#include "cuda_bnd.cuh"
#include "balance.hxx"

// ----------------------------------------------------------------------
// ctor

template <typename MF>
BndCuda3<MF>::BndCuda3(const Grid_t& grid, const int ibn[3])
{
  if (!cbnd_) {
    cbnd_ = new CudaBnd{grid, Int3::fromPointer(ibn)};
  }
}

// ----------------------------------------------------------------------
// dtor

template <typename MF>
BndCuda3<MF>::~BndCuda3()
{
  // FIXME, if we're the last user, we should clean up cbnd_?
}

// ----------------------------------------------------------------------
// add_ghosts

template <typename MF>
void BndCuda3<MF>::add_ghosts(const Grid_t& grid, storage_type& mflds_gt,
                              const Int3& mflds_ib, int mb, int me)
{
  cbnd_->add_ghosts(grid, mflds_gt, mflds_ib, mb, me);
}

template <typename MF>
void BndCuda3<MF>::add_ghosts(Mfields& mflds, int mb, int me)
{
  add_ghosts(mflds.grid(), mflds.storage(), mflds.ib(), mb, me);
}

// ----------------------------------------------------------------------
// fill_ghosts

template <typename MF>
void BndCuda3<MF>::fill_ghosts(const Grid_t& grid, storage_type& mflds_gt,
                               const Int3& mflds_ib, int mb, int me)
{
  cbnd_->fill_ghosts(grid, mflds_gt, mflds_ib, mb, me);
}

template <typename MF>
void BndCuda3<MF>::fill_ghosts(Mfields& mflds, int mb, int me)
{
  fill_ghosts(mflds.grid(), mflds.storage(), mflds.ib(), mb, me);
}

template <typename MF>
void BndCuda3<MF>::clear()
{
  if (cbnd_) {
    cbnd_->clear();
  }
}

template <typename MF>
CudaBnd* BndCuda3<MF>::cbnd_;

template struct BndCuda3<MfieldsCuda>;
template struct BndCuda3<MfieldsStateCuda>;
