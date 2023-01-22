
#include "cuda_mparticles.hxx"

#include "psc_bits.h"
#include "cuda_bits.h"
#include "bs.hxx"

// ======================================================================
// cuda_mparticles_base

// ----------------------------------------------------------------------
// ctor

template <typename BS, typename S>
cuda_mparticles_base<BS, S>::cuda_mparticles_base(const Grid_t& grid)
  : cuda_mparticles_indexer<BS>(grid), grid_(grid), by_block_(this->n_blocks)
{}

// ----------------------------------------------------------------------
// reserve_all

template <typename BS, typename S>
void cuda_mparticles_base<BS, S>::resize(uint size)
{
  storage.resize(size);
}

// ----------------------------------------------------------------------
// sizeByPatch

template <typename BS, typename S>
std::vector<uint> cuda_mparticles_base<BS, S>::sizeByPatch() const
{
  std::vector<uint> n_prts_by_patch(this->n_patches());
  thrust::host_vector<uint> h_off(by_block_.d_off);

  for (int p = 0; p < this->n_patches(); p++) {
    n_prts_by_patch[p] = h_off[(p + 1) * this->n_blocks_per_patch] -
                         h_off[p * this->n_blocks_per_patch];
    // printf("p %d n_prts_by_patch %d\n", p, n_prts_by_patch[p]);
  }
  return n_prts_by_patch;
}

template struct cuda_mparticles_base<BS144, MparticlesCudaStorage>;
template struct cuda_mparticles_base<BS444, MparticlesCudaStorage>;
