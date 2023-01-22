
#include "bnd_cuda_3_impl.hxx"
#include "cuda_bnd.hxx"
#include "balance.hxx"

// ----------------------------------------------------------------------
// ctor

BndCuda3::BndCuda3()
{
  if (!cbnd_) {
    cbnd_ = new CudaBnd{};
  }
}

// ----------------------------------------------------------------------
// dtor

BndCuda3::~BndCuda3()
{
  // FIXME, if we're the last user, we should clean up cbnd_?
}

// ----------------------------------------------------------------------
// add_ghosts

template <typename S>
void BndCuda3::add_ghosts(const Grid_t& grid, S& mflds_gt, const Int3& mflds_ib,
                          int mb, int me)
{
  cbnd_->add_ghosts(grid, mflds_gt, mflds_ib, mb, me);
}

template <typename MF>
void BndCuda3::add_ghosts(MF& mflds, int mb, int me)
{
  add_ghosts(mflds.grid(), mflds.storage(), mflds.ib(), mb, me);
}

// ----------------------------------------------------------------------
// fill_ghosts

template <typename S>
void BndCuda3::fill_ghosts(const Grid_t& grid, S& mflds_gt,
                           const Int3& mflds_ib, int mb, int me)
{
  cbnd_->fill_ghosts(grid, mflds_gt, mflds_ib, mb, me);
}

template <typename MF>
void BndCuda3::fill_ghosts(MF& mflds, int mb, int me)
{
  fill_ghosts(mflds.grid(), mflds.storage(), mflds.ib(), mb, me);
}

void BndCuda3::clear()
{
  if (cbnd_) {
    cbnd_->clear();
  }
}

CudaBnd* BndCuda3::cbnd_;

template void BndCuda3::add_ghosts(MfieldsCuda&, int, int);
template void BndCuda3::add_ghosts(MfieldsStateCuda&, int, int);
template void BndCuda3::add_ghosts(const Grid_t&,
                                   psc::gtensor_device<float, 5>&, const Int3&,
                                   int, int);
template void BndCuda3::fill_ghosts(MfieldsCuda&, int, int);
template void BndCuda3::fill_ghosts(MfieldsStateCuda&, int, int);
template void BndCuda3::fill_ghosts(const Grid_t&,
                                    psc::gtensor_device<float, 5>&, const Int3&,
                                    int, int);
