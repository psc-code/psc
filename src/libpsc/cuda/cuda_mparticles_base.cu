
#include "cuda_mparticles.cuh"

#include "psc_bits.h"
#include "cuda_bits.h"
#include "bs.hxx"

// ======================================================================
// cuda_mparticles_base

// ----------------------------------------------------------------------
// ctor

template<typename BS>
cuda_mparticles_base<BS>::cuda_mparticles_base(const Grid_t& grid)
  : cuda_mparticles_indexer<BS>(grid),
    grid_(grid),
    by_block_(this->n_blocks)
{}

// ----------------------------------------------------------------------
// reserve_all
  
template<typename BS>
void cuda_mparticles_base<BS>::resize(uint size)
{
  storage.xi4.resize(size);
  storage.pxi4.resize(size);
}

// ----------------------------------------------------------------------
// get_size_all

template<typename BS>
void cuda_mparticles_base<BS>::get_size_all(uint *n_prts_by_patch)
{
  thrust::host_vector<uint> h_off(by_block_.d_off);

  for (int p = 0; p < this->n_patches(); p++) {
    n_prts_by_patch[p] = h_off[(p+1) * this->n_blocks_per_patch] - h_off[p * this->n_blocks_per_patch];
    //printf("p %d n_prts_by_patch %d\n", p, n_prts_by_patch[p]);
  }
}

template struct cuda_mparticles_base<BS144>;
template struct cuda_mparticles_base<BS444>;